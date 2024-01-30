using Trapz             # v2.0.3
using SpecialFunctions  # v2.1.7
using Interpolations    # v0.14.4
using FFTW              # v2.0.3

mas2rad(mas) = mas * π / 6.48e8
shift(u, v, Vis, α, β) = @. Vis * cispi(-2*mas2rad(α*u + β*v))

function baseline(u, v, inc, PA)
    _u = @. (cosd(PA) * u - sind(PA) * v) * cosd(inc)
    _v = @. sind(PA) * u + cosd(PA) * v
    return _u, _v
end

function gauss(u, v, fwhm)
    return exp.(-(π * mas2rad.(fwhm) .* hypot.(u, v)).^2 / (4 * log(2)))
end

function disk(u, v, diam)
    Θ = mas2rad.(diam)
    _B = hypot.(u, v)
    Vis = 2 * besselj1.(π * Θ .* _B) ./ (π * Θ .* _B)
    Vis[iszero.(_B)] .= 1
    return Vis
end

function annulus(u, v, diam)
    return besselj0.(π * mas2rad.(diam) .* hypot.(u, v))
end

function pseudolorentzian(u, v, hfdiam)
    return exp.(-π / sqrt(3) * mas2rad.(hfdiam) .* hypot.(u, v))
end

function generalVis(u, v, I, r)
    mr = mas2rad.(r)
    _B = reshape(hypot.(u, v), 1, size(u)...)
    return trapz(mr, I .* mr .* besselj0.(2π * mr .* _B), Val(1))/trapz(mr, I .* mr)
end

function generalVis(u, v, I, r, n::Vector{<:Real}, phi::Vector{<:Real}, amp::Vector{<:Real})
    mr = mas2rad.(r)
    _B = reshape(hypot.(u, v), 1, size(u)...)
    _PA = atand.(v, u)
    
    Vis = trapz(mr, I .* mr .* besselj0.(2π * mr .* _B), Val(1))/trapz(mr, I .* mr) .|> Complex

    Hankel(n) = trapz(mr, I .* mr .* besselj.(n, 2π * mr .* _B), Val(1))/trapz(mr, I .* mr)

    for i in eachindex(n)
        if amp[i] > 0
            Vis .+= amp[i] * (-im)^n[i] .* Hankel(n[i]) .* cosd.(n[i] * (_PA .+ (phi[i] - 90)))
        end
    end

    return Vis
end

function generalVis(u, v, I, r, phi::Real, amp::Real)
    mr = mas2rad.(r)
    _B = reshape(hypot.(u, v), 1, size(u)...)
    _PA = atand.(v, u)
    Vis = trapz(mr, I .* mr .* besselj0.(2π * mr .* _B), Val(1))/trapz(mr, I .* mr)
    hankel = trapz(mr, I .* mr .* besselj1.(2π * mr .* _B), Val(1))/trapz(mr, I .* mr)
    return Complex.(Vis, -amp .* hankel .* cosd.(_PA .+ (phi - 90)))
end

function img2vis(img, pixscale; normalize=true)
    s = size(img)
    npix = s[1]
    @assert npix == s[2]

    uvspan = ((0:(npix-1)) .- (npix ÷ 2)) .* inv(mas2rad(2*pixscale) * (npix÷2))
    cvis = fftshift(fft(ifftshift(img, (1,2)), (1,2)), (1,2))
    if normalize
        cvis ./= sum(img, dims=(1,2))
    end
    return cvis
    return interpolate((uvspan, uvspan), cvis, Gridded(Linear()))
end

function img2vis(img, pixscale, plan; normalize=true)
    s = size(img)
    npix = s[1]
    @assert npix == s[2]

    uvspan = ((0:(npix-1)) .- (npix ÷ 2)) .* inv(mas2rad(2*pixscale) * (npix÷2))
    cvis = fftshift(plan * (ifftshift(img, (1,2))), (1,2))
    if normalize
        cvis ./= sum(img, dims=(1,2))
    end
    return interpolate((uvspan, uvspan), cvis, Gridded(Linear()))
end

function vis2img(vis, uvscale; normalize=true, rmstar=false, itp=true)
    s = size(vis)
    npix = s[1]
    @assert npix == s[2]

    pixspan = ((0:(npix-1)) .- (npix ÷ 2)) .* inv(mas2rad(2*uvscale) * (npix÷2))
    img = ifftshift(ifft(ifftshift(vis, (1,2)), (1,2)), (1,2)) .|> real
    if rmstar
        i = npix÷2 + 1
        img[i, i] = sum([img[i-1, i], img[i+1, i], img[i, i-1], img[i, i+1]]) / 4
    end
    if normalize
        img ./= sum(img, dims=(1,2))
    end
    return if itp
        interpolate((pixspan, pixspan), img, Gridded(Linear()))
    else
        img
    end
end

function intensity_profile(rad_prof::Function, az_prof::Function=θ -> 1; inc=0, PA=0, α0=0, δ0=0)
    function _intensity_profile(α, δ)
        x, y = [cosd(-PA)*secd(inc) sind(-PA)*secd(inc); -sind(-PA) cosd(PA)] * [α - α0, δ - δ0]
        return rad_prof(hypot(x, y)) * (1 + az_prof(-atand(x, y)))
    end
end

function gauss_profile(fwhm)
    σ = fwhm / (2 * sqrt(2 * log(2)))
    return r -> exp(-(r^2)/(2*σ^2)) / (2π * σ^2)
end

function img_extended_emission(p; fov=10, pixscale=0.025)
    inc,  PA,  fwhm,  g0,  gX,  h0,  hX = p

    pixspan = (0:pixscale:fov) .- (fov/2)

    im_disk = intensity_profile(gauss_profile(fwhm); inc, PA).(pixspan, pixspan')
    im_disk[isnan.(im_disk)] .= 0
    im_disk = im_disk ./ sum(im_disk)
    imdisk(λ) = im_disk .* (g0 * (λ/1.6e-6)^gX)
    imhalo(λ) = (h0 * (λ/1.6e-6)^hX) / length(pixspan)^2

    imext(λ) = imdisk(λ) .+ imhalo(λ)

    return imext(1.6e-6)
end
