import Base: getindex
using OIFITS    # v1.0.0-a `https://github.com/emmt/OIFITS.jl.git#master`

struct Vis2{T<:AbstractFloat}
    data::Vector{T}
    error::Vector{T}
    ucoord::Vector{T}
    vcoord::Vector{T}
    mjd::Vector{T}
    λ::Vector{T}
    instr::Vector{Symbol}
    function Vis2(data::Vector{T}, error::Vector{T}, ucoord::Vector{T}, vcoord::Vector{T}, mjd::Vector{T}, λ::Vector{T}, instr::Vector{Symbol}) where T<:AbstractFloat
        if all(length(data) .== length.([error, ucoord, vcoord, mjd, λ, instr]))
            new{eltype(data)}(data, error, ucoord, vcoord, mjd, λ, instr)
        else
            error("Failed to create Vis2 object. Arrays are of different lengths")
        end
    end
end

function Vis2(data::T, error::T, ucoord::T, vcoord::T, mjd::T, λ::T, instr::Symbol) where T<:AbstractFloat
    return Vis2([data], [error], [ucoord], [vcoord], [mjd], [λ], [instr])
end

function getindex(a::Vis2, @nospecialize vals)
    return Vis2(a.data[vals], a.error[vals], a.ucoord[vals], a.vcoord[vals], a.mjd[vals], a.λ[vals], a.instr[vals])
end

function set_error!(v2::Vis2; vis2err_pmin=0.066, vis2err_min=2e-4)
    if isnothing(vis2err_pmin) && isnothing(vis2err_min)
        nothing
    elseif isnothing(vis2err_pmin)
        for i in eachindex(v2.error)
            v2.error[i] = max(vis2err_min, v2.error[i])
        end
    elseif isnothing(vis2err_min)
        for i in eachindex(v2.error)
            v2.error[i] = max(v2.data[i] * vis2err_pmin, v2.error[i])
        end
    else
        for i in eachindex(v2.error)
            v2.error[i] = max(v2.data[i] * vis2err_pmin, vis2err_min, v2.error[i])
        end
    end
end

function set_error(v2::Vis2; vis2err_max=0.2, vis2err_pmax=0.2, vis2err_min=2e-4, vis2err_pmin=0.066)
    amask = if isnothing(vis2err_max)
        ones(Bool, length(v2.data))
    else
        v2.error .≤ vis2err_max
    end
    pmask = if isnothing(vis2err_pmax)
        ones(Bool, length(v2.data))
    else
        (v2.error ./ v2.data) .≤ vis2err_pmax
    end
    mask = amask .& pmask
    _v2 = Vis2(v2.data[mask], v2.error[mask], v2.ucoord[mask], v2.vcoord[mask], v2.mjd[mask], v2.λ[mask], v2.instr[mask])
    set_error!(_v2; vis2err_min, vis2err_pmin)
    return _v2
end

struct CP{T<:AbstractFloat}
    ampdata::Vector{T}
    amperror::Vector{T}
    data::Vector{T}
    error::Vector{T}
    u1coord::Vector{T}
    v1coord::Vector{T}
    u2coord::Vector{T}
    v2coord::Vector{T}
    u3coord::Vector{T}
    v3coord::Vector{T}
    mjd::Vector{T}
    λ::Vector{T}
    instr::Vector{Symbol}
    function CP(ampdata::Vector{T}, amperror::Vector{T}, data::Vector{T}, error::Vector{T}, u1coord::Vector{T}, v1coord::Vector{T}, u2coord::Vector{T}, v2coord::Vector{T}, u3coord::Vector{T}, v3coord::Vector{T}, mjd::Vector{T}, λ::Vector{T}, instr::Vector{Symbol}) where T<:AbstractFloat
        if all(length(data) .== length.([ampdata, amperror, error, u1coord, v1coord, u2coord, v2coord, u3coord, v3coord, mjd, λ, instr]))
            new{eltype(data)}(ampdata, amperror, data, error, u1coord, v1coord, u2coord, v2coord, u3coord, v3coord, mjd, λ, instr)
        else
            error("Failed to create CP object. Arrays are of different lengths")
        end
    end
end

function CP(ampdata::T, amperror::T, data::T, error::T, u1coord::T, v1coord::T, u2coord::T, v2coord::T, u3coord::T, v3coord::T, mjd::T, λ::T, instr::Symbol) where T<:AbstractFloat
    return CP([ampdata], [amperror], [data], [error], [u1coord], [v1coord], [u2coord], [v2coord], [u3coord], [v3coord], [mjd], [λ], [instr])
end

function getindex(a::CP, @nospecialize vals)
    return CP(a.ampdata[vals], a.amperror[vals], a.data[vals], a.error[vals], a.u1coord[vals], a.v1coord[vals], a.u2coord[vals], a.v2coord[vals], a.u3coord[vals], a.v3coord[vals], a.mjd[vals], a.λ[vals], a.instr[vals])
end

function set_error!(cp::CP; cperr_min=0.1)
    if !isnothing(cperr_min)
        for i in eachindex(cp.error)
            cp.error[i] = max(cperr_min, cp.error[i])
        end
    end
end

function set_error(cp::CP; cperr_max=10, cperr_min=0.1)
    mask = if isnothing(cperr_max)
        1:length(cp.data)
    else
        cp.error .≤ cperr_max
    end
    _cp = CP(cp.ampdata[mask], cp.amperror[mask], cp.data[mask], cp.error[mask], cp.u1coord[mask], cp.v1coord[mask], cp.u2coord[mask], cp.v2coord[mask], cp.u3coord[mask], cp.v3coord[mask], cp.mjd[mask], cp.λ[mask], cp.instr[mask])
    set_error!(_cp; cperr_min)
    return _cp
end

function read_files(files...; hack_revn=1, exclude_edge=false)
    v2_data = OIVis2[]
    cp_data = OIT3[]

    for file in files
        f = read(OIDataSet, file; hack_revn)
        for v2table in f.vis2
            push!(v2_data, v2table)
        end
        for cptable in f.t3
            push!(cp_data, cptable)
        end
    end

    _vis2data = Float64[]
    _vis2err  = Float64[]
    _ucoord   = Float64[]
    _vcoord   = Float64[]
    _vis2mjd  = Float64[]
    _vis2λ    = Float64[]
    _vis2ins  = Symbol[]

    for v2 in v2_data
        for i in eachindex(v2.vis2err)
            if isnan(v2.vis2err[i]) || isnan(v2.vis2data[i])
                v2.flag[i] = true
            end
        end
        if exclude_edge
            v2.flag[[1,end],:] .= true
        end
        mask = .!vec(v2.flag)
        nλ = length(v2.instr.eff_wave)
        push!(_vis2data, v2.vis2data[mask]...)
        push!(_vis2err, v2.vis2err[mask]...)
        push!(_ucoord, repeat(v2.ucoord, inner=nλ)[mask]...)
        push!(_vcoord, repeat(v2.vcoord, inner=nλ)[mask]...)
        push!(_vis2mjd, repeat(v2.mjd, inner=nλ)[mask]...)
        push!(_vis2λ, repeat(v2.instr.eff_wave, outer=length(v2.mjd))[mask]...)
        push!(_vis2ins, repeat([Symbol(lowercase(prod(split(split(v2.instr.insname, "_")[1], "-"))))], sum(mask))...)
    end

    _cadata   = Float64[]
    _caerr    = Float64[]
    _cpdata   = Float64[]
    _cperr    = Float64[]
    _u1coord  = Float64[]
    _v1coord  = Float64[]
    _u2coord  = Float64[]
    _v2coord  = Float64[]
    _cpmjd    = Float64[]
    _cpλ      = Float64[]
    _cpins    = Symbol[]

    for cp in cp_data
        for i in eachindex(cp.t3phierr)
            if isnan(cp.t3phierr[i]) || isnan(cp.t3phi[i]) || isnan(cp.t3amperr[i]) || isnan(cp.t3amp[i])
                cp.flag[i] = true
            end
        end
        if exclude_edge
            cp.flag[[1,end],:] .= true
        end
        mask = .!vec(cp.flag)
        nλ = length(cp.instr.eff_wave)
        push!(_cadata, cp.t3amp[mask]...)
        push!(_caerr, cp.t3amperr[mask]...)
        push!(_cpdata, cp.t3phi[mask]...)
        push!(_cperr, cp.t3phierr[mask]...)
        push!(_u1coord, repeat(cp.u1coord, inner=nλ)[mask]...)
        push!(_v1coord, repeat(cp.v1coord, inner=nλ)[mask]...)
        push!(_u2coord, repeat(cp.u2coord, inner=nλ)[mask]...)
        push!(_v2coord, repeat(cp.v2coord, inner=nλ)[mask]...)
        push!(_cpmjd, repeat(cp.mjd, inner=nλ)[mask]...)
        push!(_cpλ, repeat(cp.instr.eff_wave, outer=length(cp.mjd))[mask]...)
        push!(_cpins, repeat([Symbol(lowercase(prod(split(split(cp.instr.insname, "_")[1], "-"))))], sum(mask))...)
    end

    return Vis2(_vis2data, _vis2err, _ucoord ./ _vis2λ, _vcoord ./ _vis2λ, _vis2mjd, _vis2λ, _vis2ins), CP(_cadata, _caerr, _cpdata, _cperr, _u1coord ./ _cpλ, _v1coord ./ _cpλ, _u2coord ./ _cpλ, _v2coord ./ _cpλ, (_u1coord .+ _u2coord) ./ _cpλ, (_v1coord .+ _v2coord) ./ _cpλ, _cpmjd, _cpλ, _cpins)
end

function read_files_set_error(files...; vis2err_max=nothing, vis2err_pmax=0.2, vis2err_min=2e-4, vis2err_pmin=0.066, cperr_max=10, cperr_min=0.1, kwargs...)
    v2, cp = read_files(files...; kwargs...)
    return set_error(v2; vis2err_min, vis2err_pmin, vis2err_max, vis2err_pmax), set_error(cp; cperr_min, cperr_max)
end
