include("fits.jl")

filepath = "/path/to/reduced/oifits/files"
savedir = "/path/to/save/output/products/to"

files = Dict{String, Vector{String}}()

# change the following as needed:

files["T1"] = [
    "pionier_july29.oifits",
    "mircx_june5.oifits",
]

files["T2"] = [
    "pionier_july29.oifits",
    "mircx_june6.oifits",
]

files["T3"] = [
    "pionier_july29.oifits",
    "mircx_june8_pol1.oifits",
    "mircx_june8_pol2.oifits",
]

files["T4"] = [
    "pionier_july29.oifits",
    "pionier_june9.oifits",
    "mircx_june9.oifits",
]

files["T5"] = [
    "pionier_july29.oifits",
    "mircx_june10_pol1.oifits",
    "mircx_june10_pol2.oifits",
]

files["T6"] = [
    "pionier_july29.oifits",
    "pionier_july9.oifits",
]

files["T7"] = [
    "pionier_july29.oifits",
    "mircx_july13.oifits",
]

files["T8"] = [
    "pionier_july29.oifits",
    "pionier_july20.oifits",
    "mircx_july20.oifits",
]

for k in keys(files)
    v2[k], cp[k] = read_files_set_error(
        (joinpath.(filepath, files[k]))...;
        vis2err_max=0.2,
        vis2err_pmax=nothing,
        vis2err_min=2e-4,
        vis2err_pmin=0.066,
        cperr_max=20,
        cperr_min=0.1
    )
end

# param_name  min     max     wrap (i.e. whether values can "pacman" from max to/from min)
setup = [
    :inc    0       90      false
    :PA     0       180     true
    :ringr  0       5       false
    :ringw  0       5       false
    :r0     0       6       false
    :rX     -4      10      false
    :h0     0       1       false
    :ringa  0       1       false
    :ringϕ  0       360     true
    :ringa2 0       1       false
    :ringϕ2 0       (360/2) true
    :ringa3 0       1       false
    :ringϕ3 0       (360/3) true
    :ringa4 0       1       false
    :ringϕ4 0       (360/4) true
    :ringa5 0       1       false
    :ringϕ5 0       (360/5) true
    :ringa6 0       1       false
    :ringϕ6 0       (360/6) true
    :ringa7 0       1       false
    :ringϕ7 0       (360/7) true
    :ringa8 0       1       false
    :ringϕ8 0       (360/8) true
]

for s in ["T$ep" for ep in 1:8]
    model = "centered_skew_ring"

    i = 5 # maximum azimuthal order

    file_str = "$(model)_$(i)_$(s)"

    r1 = Dict()

    _setup = setup[3:(7 + 2*i), :] # inc/pa fixed for "centered_skew_ring"

    r1["result"] = unest(_setup, cvis[model], s; min_num_live_points=1000)

    save(joinpath(savedir, "$file_str.jld2"), r1)
end
