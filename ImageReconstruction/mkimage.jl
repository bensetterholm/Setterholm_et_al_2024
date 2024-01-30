ENV["MPLBACKEND"]="qt5agg"
using PyPlot
using OITOOLS
using FITSIO

function mkimage(
    filename,
    pixscale = 1/16,
    nx = 2^8;
    regularizers=[["tvsq", 1e7], ["compactness", 1e9]],
    params_start=[0.36, 0.05, 0.0, 4.0, 1.6e-6],
    x_start = rand(nx*nx)/(nx*nx*4),
    savedir=dirname(filename),
    saveimg=true
)
    data = readoifits(filename, filter_bad_data=true)[1,1];
    ft = setup_nfft(data, nx, pixscale);
    params, x = reconstruct_sparco_gray(x_start, params_start, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=100); #grey environment
    minchi2, params,ret = optimize_sparco_parameters(params, x, ft, data; weights = [1.0,1.0,1.0] )
    params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=100); #grey environment
    minchi2, params,ret = optimize_sparco_parameters(params, x, ft, data; weights = [1.0,1.0,1.0] )
    params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=100); #grey environment

    imdisp(x, pixscale = pixscale)
    splitfile = split(basename(filename), '.')
    length(splitfile) > 1 && any(splitfile[end] |> lowercase .== ["fits", "oifits"]) && pop!(splitfile)
    outfileprefix = joinpath(savedir, join(splitfile, '.'))
    saveimg && savefig("$outfileprefix.png")
    f = FITS("$outfileprefix.fits", "w")
    _x = reshape(x, (nx, nx))
    hdr_key = [
        "PIXSCALE",
        "MINCHI2",
        "FITMSG",
        "STRFRAC",
        "BGFRAC",
        "STRDIAM",
        "ENVSLOPE",
        "REFWL",
    ]
    hdr_val = Any[
        pixscale,
        minchi2,
        string(ret),
        params...
    ]
    hdr_cmt = [
        "mas",
        "",
        "see OITOOLS",
        "flux fraction of star at lamdba_0",
        "flux fraction of background at lambda_0",
        "stellar angular diameter",
        "spectral index of the environment plus 4 (gt 0)",
        "lambda_0 reference wavelength",
    ]
    for (i, reg) in enumerate(regularizers)
        push!(hdr_key, "REG$i")
        push!(hdr_val, reg[2])
        push!(hdr_cmt, reg[1])
    end
    write(f, _x; header=FITSHeader(hdr_key, hdr_val, hdr_cmt))
    close(f)
    return
end
