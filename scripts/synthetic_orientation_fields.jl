# Synthetic channel maps whose structures follow a prescribed orientation field t(x, y), for measuring what VGT does
# to orientation structure (effective resolution, accuracy, noise uncertainty). Included by
# scripts/validate_resolution_uncertainty.jl.
#
# Each channel is a Gaussian random field generated in Fourier space with power k^-beta and anisotropy
# exp(-zeta cos^2(phi_k - t)) (wave vectors along t suppressed, so structures run along t; the same generator as
# scripts/validate_angle_bias.jl). It is generated with the same random phases at `nbins` angles t_m = m pi / nbins
# and, at each pixel, blended between the two bins that bracket t(x, y) (linear weights, variance restored with the
# measured correlation of the two bins). The blend is applied before the beam, so the orientation field is on the
# sky and the telescope beam smooths it like any other structure. Then the beam (Gaussian, FWHM `beam` px) and a
# scale that sets S/N r (rms of the signal minus its 20 px boxcar mean, in units of the noise rms).

using FFTW
using Random
using Statistics

beam_tf(kx, ky, fwhm) = exp(-0.5 * (fwhm / sqrt(8log(2)))^2 * (kx^2 + ky^2))
kgrid(n) = 2pi .* fftfreq(n)

function highpass_tf(kx, ky; L=20)
    d(k) = k == 0 ? 1.0 : sin(L * k / 2) / (L * sin(k / 2))
    return 1 - d(kx) * d(ky)
end

"""
    orientation_signal(rng, tmap; nv=8, beta=8/3, zeta=3.0, beam=3.0, nbins=72, r=2.0)

Signal cube (n, n, nv) following the orientation field `tmap` (radians, structures along t, so the VGT angle should
equal t). Every channel is scaled so its 20 px high-passed rms is `r` (unit noise rms; `r = Inf` scales it to 2).
"""
function orientation_signal(rng::AbstractRNG, tmap::AbstractMatrix; nv=8, beta=8 / 3, zeta=3.0, beam=3.0, nbins=72,
                            r=2.0)
    n = size(tmap, 1)
    size(tmap, 2) == n || throw(ArgumentError("tmap must be square"))
    k = kgrid(n)
    kk = [hypot(a, b) for a in k, b in k]
    phik = [atan(b, a) for a in k, b in k]
    base = map(x -> x == 0 ? 0.0 : x^(-beta / 2), kk)
    btf = [beam_tf(a, b, beam) for a in k, b in k]
    hptf = [highpass_tf(a, b) for a in k, b in k]
    Δ = pi / nbins
    p = mod.(tmap, pi) ./ Δ
    m0 = mod.(floor.(Int, p), nbins)
    w = p .- floor.(p)
    m1 = mod.(m0 .+ 1, nbins)
    pairs = unique(zip(vec(m0), vec(m1)))
    out = zeros(n, n, nv)
    for v in 1:nv
        ph = fft(randn(rng, n, n))
        F = Dict{Int, Matrix{Float64}}()
        field(m) = get!(F, m) do
            f = real.(ifft(ph .* base .* exp.(-0.5 * zeta .* cos.(phik .- m * Δ) .^ 2)))
            f ./ std(f)
        end
        rho = Dict(pr => cor(vec(field(pr[1])), vec(field(pr[2]))) for pr in pairs)
        s = Matrix{Float64}(undef, n, n)
        @inbounds for idx in eachindex(s)
            a, b, ww = field(m0[idx]), field(m1[idx]), w[idx]
            ρ = rho[(m0[idx], m1[idx])]
            s[idx] = ((1 - ww) * a[idx] + ww * b[idx]) / sqrt((1 - ww)^2 + ww^2 + 2ww * (1 - ww) * ρ)
        end
        fs = fft(s) .* btf
        hp = real.(ifft(fs .* hptf))
        out[:, :, v] .= real.(ifft(fs)) .* ((isinf(r) ? 2.0 : r) / std(hp))
    end
    return out
end

"Sinusoidal orientation field t = t0 + a sin(2 pi i / lambda) along axis 1 (i = 1..n), constant along axis 2."
sinusoid_tmap(n, t0, a, lambda) = [t0 + a * sin(2pi * i / lambda) for i in 1:n, j in 1:n]

"""
    random_tmap(rng, n; t0=0.5, delta=0.8, alpha=3.0)

Random orientation field: Q + iU = exp(2i t0) + delta * z, z a complex Gaussian random field with power k^-alpha
and unit rms per component, t = arg(Q + iU) / 2.
"""
function random_tmap(rng::AbstractRNG, n; t0=0.5, delta=0.8, alpha=3.0)
    k = kgrid(n)
    amp = [hypot(a, b) == 0 ? 0.0 : hypot(a, b)^(-alpha / 2) for a in k, b in k]
    zq = real.(ifft(fft(randn(rng, n, n)) .* amp))
    zu = real.(ifft(fft(randn(rng, n, n)) .* amp))
    zq ./= std(zq)
    zu ./= std(zu)
    return 0.5 .* atan.(sin(2t0) .+ delta .* zu, cos(2t0) .+ delta .* zq)
end
