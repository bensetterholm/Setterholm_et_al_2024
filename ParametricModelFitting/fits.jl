include("data.jl")
include("models.jl")

using OrderedCollections    # v1.4.1
using StatsBase             # v0.33.21
using Dates                 # Julia Standard Library
using UltraNest             # v0.1.0 (with python ultranest v3.3.3)
using Measurements          # v2.7.2
using Distributions         # v0.25.66
using PyCall                # v1.93.1 (with anaconda python v3.7.3)
using FileIO                # v1.15.0 (with JLD2 v0.4.23)
using Printf                # Julia Standard Library
using CMPFit                # v0.3.3
using ProgressMeter
stepsampler = pyimport("ultranest.stepsampler")

v2 = Dict{String, Vis2{Float64}}()
cp = Dict{String, CP{Float64}}()

function get_params(res)
    distfit = [Distributions.fit(Normal, res["samples"][:,s]) for s in 1:length(res["paramnames"])]
    return [f.μ ± f.σ for f in distfit]
end

cvis = Dict{String, Function}()

cvis["gaussian"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc,  PA,  fwhm,  g0,  gX,  h0 = params
    hX = 0

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fgauss = @. g0 * (λ / λ0) ^ gX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        g = fgauss .* gauss(_u, _v, fwhm)
        return (fstar .+ g) ./ (fstar .+ fhalo .+ fgauss)
    end
end

cvis["skew_ring_fix"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc, PA, ringr, ringw, ringa, h0, r0, rX = params
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(besselj0.(π * d .* _B), -ringa * besselj1.(π * d .* _B) .* cosd.(_PA)) .* gauss(_u, _v, ringw)
        return (fstar .+ ring) ./ (fstar .+ fhalo .+ fring)
    end
end

cvis["skew_ring"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc, PA, ringr, ringw, ringa, ringϕ, h0, r0, rX = params
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(besselj0.(π * d .* _B), -ringa * besselj1.(π * d .* _B) .* cosd.(_PA .+ (ringϕ - 90))) .* gauss(_u, _v, ringw)
        return (fstar .+ ring) ./ (fstar .+ fhalo .+ fring)
    end
end

cvis["skew_ring2"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc, PA, ringr, ringw, ringa, ringϕ, ringa2, ringϕ2, h0, r0, rX = params
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(besselj0.(π * d .* _B) .- (ringa2 * besselj.(2, π * d .* _B) .* cosd.(2 * (_PA .+ (ringϕ2 - 90)))), -ringa * besselj1.(π * d .* _B) .* cosd.(_PA .+ (ringϕ - 90))) .* gauss(_u, _v, ringw)
        return (fstar .+ ring) ./ (fstar .+ fhalo .+ fring)
    end
end

cvis["skew_ring_fix_plus_blob"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc, PA, ringr, ringw, ringa, blobw, blobα, blobδ, h0, r0, rX, b0, bX = params
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX
    fblob = @. b0 * (λ / λ0) ^ bX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(besselj0.(π * d .* _B), -ringa * besselj1.(π * d .* _B) .* cosd.(_PA)) .* gauss(_u, _v, ringw)
        blob = shift(u, v, fblob .* gauss(u, v, blobw), blobα, blobδ)
        return (fstar .+ ring .+ blob) ./ (fstar .+ fhalo .+ fring .+ fblob)
    end
end

cvis["skew_ring_extra_fix_plus_blob"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    ringr, ringw, ringa, blobw, blobα, blobδ, h0, r0, rX, b0, bX = params
    inc = 48
    PA = 136
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX
    fblob = @. b0 * (λ / λ0) ^ bX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(besselj0.(π * d .* _B), -ringa * besselj1.(π * d .* _B) .* cosd.(_PA)) .* gauss(_u, _v, ringw)
        blob = shift(u, v, fblob .* gauss(u, v, blobw), blobα, blobδ)
        return (fstar .+ ring .+ blob) ./ (fstar .+ fhalo .+ fring .+ fblob)
    end
end

cvis["centered_skew_ring"] = function (params, λ=1.6e-6; λ0=1.6e-6)
    inc = 46
    PA = 133
    rorder = (length(params) - 5) ÷ 2
    ringr, ringw, r0, rX, h0, ringp... = params
    ring_params = zeros(8 * 2)
    ring_params[1:length(ringp)] = ringp
    ringa, ringϕ, ringa2, ringϕ2, ringa3, ringϕ3, ringa4, ringϕ4, ringa5, ringϕ5, ringa6, ringϕ6, ringa7, ringϕ7, ringa8, ringϕ8 = ring_params
    hX = 0

    d = mas2rad(2*ringr)

    fstar = @. (λ / λ0) ^ -4
    fhalo = @. h0 * (λ / λ0) ^ hX
    fring = @. r0 * (λ / λ0) ^ rX

    function (u, v)
        _u, _v = baseline(u, v, inc, PA)
        _PA = atand.(_v, _u)
        _B = hypot.(_u, _v)
        ring = fring .* Complex.(
            besselj0.(π * d .* _B)
            .- (rorder ≥ 2 && (ringa2 * besselj.(2, π * d .* _B) .* cosd.(2 * (_PA .+ (ringϕ2 - 90)))))
            .+ (rorder ≥ 4 && (ringa4 * besselj.(4, π * d .* _B) .* cosd.(4 * (_PA .+ (ringϕ4 - 90)))))
            .- (rorder ≥ 6 && (ringa6 * besselj.(6, π * d .* _B) .* cosd.(6 * (_PA .+ (ringϕ6 - 90)))))
            .+ (rorder ≥ 8 && (ringa8 * besselj.(8, π * d .* _B) .* cosd.(8 * (_PA .+ (ringϕ8 - 90)))))
            , 0
            .- (rorder ≥ 1 && (ringa * besselj1.(π * d .* _B) .* cosd.(_PA .+ (ringϕ - 90))))
            .+ (rorder ≥ 3 && (ringa3 * besselj.(3, π * d .* _B) .* cosd.(3 * (_PA .+ (ringϕ3 - 90)))))
            .- (rorder ≥ 5 && (ringa5 * besselj.(5, π * d .* _B) .* cosd.(5 * (_PA .+ (ringϕ5 - 90)))))
            .+ (rorder ≥ 7 && (ringa7 * besselj.(7, π * d .* _B) .* cosd.(7 * (_PA .+ (ringϕ7 - 90)))))
            ) .* gauss(_u, _v, ringw)
        return (fstar .+ ring) ./ (fstar .+ fhalo .+ fring)
    end
end

function chi2(visfunc, v2, cp, parameters::Vector{Measurement{T}}; kwargs...) where T<:AbstractFloat
    return chi2(visfunc, v2, cp, Measurements.value.(parameters); kwargs...)
end

function report(visfunc, v2, cp, res::Dict; kwargs...)
    best = get_params(res)
    println(chi2(visfunc, v2, cp, best; kwargs...))
    for i in eachindex(params)
        println(params[i], "\t", best[i])
    end
end

function report(file_str, result, res, params=nothing; kwargs...)
    println(chi2(cvis[file_str], v2[result], cp[result], res; kwargs...))
    _params = isnothing(params) ? load("$(file_str).jld2")[result]["paramnames"] : params
    for i in eachindex(_params)
        println(_params[i], "\t", res[i])
    end
end

function chi2(visfunc, v2, cp, params; total_only = false, l2array = false, showinfo=true, penalty::Union{Function, Nothing} = nothing)
    v2mod = visfunc(params, v2.λ)(v2.ucoord, v2.vcoord) .|> abs2

    cpvis = visfunc(params, cp.λ)

    cpmod = (
        cpvis(cp.u1coord, cp.v1coord) .*
        cpvis(cp.u2coord, cp.v2coord) .*
        conj.(cpvis(cp.u3coord, cp.v3coord))
    ) .|> angle .|> rad2deg

    χv2 = (v2mod .- v2.data) ./ v2.error
    χcp = (cpmod .- cp.data) ./ cp.error
    χtot = vcat(χv2, χcp)

    # impose penalties for unphysical solutions
    penalty_arr = if isnothing(penalty)
        0
    else
        penalty_arr = sign.(χtot) * penalty(params)
    end
    return if total_only
        χtot .|> abs2 |> sum
    elseif l2array
        χtot .+ penalty_arr
    else
        lv2 = length(χv2)
        lcp = length(χcp)
        lp = length(params)
        if showinfo
            @info "length v2: $(lv2), length cp: $(lcp), length params: $(lp)"
        end
        (total=sum(abs2.(χtot .+ penalty_arr)) / (lv2 + lcp - lp), v2=sum(χv2.^2) / (lv2 - lp), cp=sum(χcp.^2) / (lcp - lp))
    end
end

function loglike(visfunc, v2, cp)
    return parameters -> -0.5 * chi2(visfunc, v2, cp, parameters, total_only=true)
end

function unest(setup, visfunc, dataset; min_num_live_points = 400)
    params = String.(setup[:,1])
    low_val = Float64.(setup[:,2])
    high_val = Float64.(setup[:,3])
    wrapped_params = Bool.(setup[:,4])
    nsteps = 2 * length(params)

    function prior_transform(cube)
        return cube .* (high_val .- low_val) .+ low_val
    end

    sampler = ultranest.ReactiveNestedSampler(params, loglike(visfunc, v2[dataset], cp[dataset]), prior_transform; wrapped_params)
    sampler.stepsampler = stepsampler.SliceSampler(; nsteps, generate_direction=stepsampler.generate_mixture_random_direction)
    try
        sampler.run(; min_num_live_points)
    catch
        sampler.results
    end
end

function lmfit(visfunc, v2, cp, params)
    cmpfit_callback(params::Vector{Float64}) = chi2(visfunc, v2, cp, params, l2array=true)
    config = CMPFit.Config()
    config.ftol=1e-7
    config.epsfcn=1e-7
    config.maxfev=5000
    config.maxiter=5000
    parinfo = CMPFit.Parinfo(length(params))
    return cmpfit(cmpfit_callback, params; config, parinfo)
end

lmfit(file_str::String, result::String) = lmfit(cvis[file_str], v2[result], cp[result], load("$(file_str).jld2")[result]["maximum_likelihood"]["point"])

function randfit(setup::Matrix{Any}, visfunc, v2, cp, point::Union{Vector{Float64}, Nothing}=nothing; randsearch::Int=0)
    bounds = Float64.(setup[:,2:3])
    lowbound = bounds[:,1]
    boundspan = diff(bounds, dims=2) |> vec
    _point = if isnothing(point)
        mean(bounds, dims=2) |> vec
    else
        point
    end
    best = lmfit(visfunc, v2, cp, _point)
    for s in 1:randsearch
        test = lmfit(visfunc, v2, cp, lowbound .+ (rand(length(boundspan)) .* boundspan))
        if test.bestnorm < best.bestnorm
            best = test
        end
    end
    return best
end

function randfit(res::Dict{Any, Any}, visfunc, v2, cp; kwargs...)
    setup = hcat(res["paramnames"], minimum(res["samples"], dims=1)', maximum(res["samples"], dims=1)', zeros(Bool, length(res["paramnames"])))
    randfit(setup, visfunc, v2, cp, res["maximum_likelihood"]["point"]; kwargs...)
end

function randfit(file_str::String, result::String; kwargs...)
    res = load(file_str * ".jld2")[result]
    return randfit(res, cvis[file_str], v2[result], cp[result]; kwargs...)
end

function bootfit(setup::Matrix{Any}, visfunc, v2, cp, point::Union{Vector{Float64}, Nothing}=nothing; randsearch::Int=50, bootstraps::Int=1000)
    restree = fill(NaN, size(setup)[1], bootstraps)
    v2l = length(v2.data)
    cpl = length(cp.data)
    rv = Random.Sampler(Random.MersenneTwister, 1:v2l)
    rc = Random.Sampler(Random.MersenneTwister, 1:cpl)
    @showprogress for strap in 1:bootstraps
        restree[:, strap] = randfit(setup, visfunc, v2[rand(rv, v2l)], cp[rand(rc, cpl)], point; randsearch).param
    end
    return restree
end

function bootfit(res::Dict{Any, Any}, visfunc, v2, cp; kwargs...)
    setup = hcat(res["paramnames"], minimum(res["samples"], dims=1)', maximum(res["samples"], dims=1)', zeros(Bool, length(res["paramnames"])))
    bootfit(setup, visfunc, v2, cp, res["maximum_likelihood"]["point"]; kwargs...)
end
