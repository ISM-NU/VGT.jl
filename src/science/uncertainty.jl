# Per-pixel noise uncertainty of the VGT angle by Monte Carlo noise propagation (opt-in; the pipeline does not call
# it, so default outputs are unchanged).
#
# Why: the angle map has no error bar. Cube noise enters every step (gradients, window fits, intensity weights, the
# noise mask, the Stokes smoothing), and the sliding windows and the Gaussian smoothing correlate neighbouring pixels,
# so a closed-form propagation needs the covariance of overlapping windows and breaks down at low S/N, where gradient
# angles become uniform and the mask changes. Instead the observed cube is perturbed with K independent noise
# realisations drawn from the cube's own noise model and rerun through the same VGT steps; the spread of the K angle
# maps about each pixel is sigma_psi, the analogue of Planck's psisig (noise only; it does not include how well
# gradients trace the field). Validation: docs/validation-resolution-uncertainty.md (coverage on synthetic cubes
# against the noise-free VGT map).

"""
    correlated_noise(rng, nx, ny, nv; fwhm_major_px, fwhm_minor_px=fwhm_major_px, angle=0.0, sigma=1.0,
                     spectral_lag1=0.0)

Gaussian noise of shape `(nx, ny, nv)`. Each channel is white noise convolved with an elliptical Gaussian beam
(FWHM `fwhm_major_px` along `angle`, radians counter-clockwise from +axis 1, and `fwhm_minor_px` across it; a FWHM of
0 means white along that axis) and rescaled to unit rms, as noise in a radio map is correlated over the beam.
`spectral_lag1` in [0, 0.5) correlates adjacent channels with a [b, 1, b] kernel giving that lag-1 correlation.
The result is multiplied by `sigma` (a number, or an `nx x ny` per-pixel rms map). The FFT is periodic on a grid
padded by 3 beam FWHM, so opposite edges are not correlated.
"""
function correlated_noise(rng::AbstractRNG, nx::Int, ny::Int, nv::Int; fwhm_major_px::Real, fwhm_minor_px::Real=fwhm_major_px,
                          angle::Real=0.0, sigma=1.0, spectral_lag1::Real=0.0)
    fwhm_major_px >= 0 && fwhm_minor_px >= 0 || throw(ArgumentError("beam FWHM must be >= 0"))
    0 <= spectral_lag1 < 0.5 || throw(ArgumentError("spectral_lag1 must be in [0, 0.5); got $spectral_lag1"))
    sigma isa Number || size(sigma) == (nx, ny) ||
        throw(ArgumentError("sigma map size $(size(sigma)) must be ($nx, $ny)"))
    pad = ceil(Int, 3 * max(fwhm_major_px, fwhm_minor_px))
    mx, my = nx + 2pad, ny + 2pad
    sa, sb = fwhm_major_px / sqrt(8log(2)), fwhm_minor_px / sqrt(8log(2))
    kx = 2pi .* fftfreq(mx)
    ky = 2pi .* fftfreq(my)
    c, s = cos(angle), sin(angle)
    tf = [exp(-0.5 * ((sa * (a * c + b * s))^2 + (sb * (-a * s + b * c))^2)) for a in kx, b in ky]
    gain = sqrt(sum(abs2, tf) / length(tf))  # rms of the filtered unit white noise
    tf ./= gain
    # [b, 1, b] / sqrt(1 + 2 b^2) has lag-1 correlation 2b / (1 + 2b^2); solve for b (smaller root).
    b = spectral_lag1 == 0 ? 0.0 : (1 - sqrt(1 - 2 * spectral_lag1^2)) / (2 * spectral_lag1)
    white = Array{Float64}(undef, nx, ny, nv + 2)
    for v in 1:(nv + 2)
        w = real.(ifft(fft(randn(rng, mx, my)) .* tf))
        white[:, :, v] .= @view w[pad + 1:pad + nx, pad + 1:pad + ny]
    end
    out = Array{Float64}(undef, nx, ny, nv)
    norm = 1 / sqrt(1 + 2b^2)
    @inbounds for v in 1:nv, j in 1:ny, i in 1:nx
        z = (b * white[i, j, v] + white[i, j, v + 1] + b * white[i, j, v + 2]) * norm
        out[i, j, v] = z * (sigma isa Number ? sigma : sigma[i, j])
    end
    return out
end

"""
    circular_spread(psis) -> (; sigma, mean_angle, n)

Per-pixel spread of a vector of axial angle maps (radians): `sigma = sqrt(-2 ln R) / 2`, the circular standard
deviation of the doubled angles halved (equal to the ordinary standard deviation for small spreads and finite for
uniform angles), with R the mean resultant length of exp(2i psi_k) over the finite maps; `mean_angle` is the
circular mean angle; `n` the number of finite maps. Pixels with fewer than two finite maps are NaN.
"""
function circular_spread(psis::AbstractVector{<:AbstractMatrix})
    sz = size(first(psis))
    sigma = fill(NaN, sz)
    mean_angle = fill(NaN, sz)
    n = zeros(Int, sz)
    @inbounds for idx in CartesianIndices(sz)
        c = 0.0
        s = 0.0
        m = 0
        for p in psis
            a = p[idx]
            isfinite(a) || continue
            c += cos(2a)
            s += sin(2a)
            m += 1
        end
        n[idx] = m
        m >= 2 || continue
        r = min(hypot(c, s) / m, 1.0)
        sigma[idx] = r > 0 ? sqrt(-2log(r)) / 2 : Inf
        mean_angle[idx] = 0.5 * atan(s, c)
    end
    return (; sigma, mean_angle, n)
end

"""
    psi_noise_mc(cube, dn, ker, noise_threshold; realizations=16, noise_sigma, fwhm_major_px,
                 fwhm_minor_px=fwhm_major_px, noise_angle=0.0, spectral_lag1=0.0, scale=1.0, rng=Xoshiro(1),
                 min_valid=100, gradient_mode=:central_checked, kwargs...)

Monte Carlo noise uncertainty of the VGT angle of an in-memory cube `(nx, ny, nv)`. Runs VGT (`VGT_QU_mv_2d` and
`finalize_vgt_products` with `dn`, `ker`, the mask threshold `noise_threshold` and `kwargs`) on the cube and on
`realizations` copies with added noise from [`correlated_noise`](@ref) (`noise_sigma` scalar or per-pixel map,
times `scale`), and returns

- `psi`: the angle map of the cube itself (as the pipeline's `psi`, gradient + 90 deg, pixel frame);
- `sigma_psi`: per-pixel spread of the perturbed maps about their circular mean, times sqrt(K / (K - 1)) (unbiased
  for the estimated mean), divided by `scale`
  (radians; linear error propagation makes the spread proportional to the added noise, so `scale < 1` probes the
  local derivative with less non-linearity, and `scale = 1` adds noise as large as the cube's);
- `offset`: circular mean of the perturbed maps minus `psi` (radians), the shift that added noise causes;
- `sigma_reliability`: standard deviation of the reliability `sqrt(Q^2 + U^2) / I` over the perturbed maps / scale;
- `n`: number of perturbed maps finite at each pixel.

The noise mask threshold stays at `noise_threshold`, as when the pipeline runs on the observed cube.
Validation and coverage: docs/validation-resolution-uncertainty.md.
"""
function psi_noise_mc(cube::AbstractArray{<:Real, 3}, dn::Int, ker::Int, noise_threshold;
                      realizations::Int=16, noise_sigma, fwhm_major_px::Real, fwhm_minor_px::Real=fwhm_major_px,
                      noise_angle::Real=0.0, spectral_lag1::Real=0.0, scale::Real=1.0,
                      rng::AbstractRNG=Xoshiro(1), min_valid::Int=100, gradient_mode::Symbol=:central_checked,
                      kwargs...)
    realizations >= 2 || throw(ArgumentError("realizations must be >= 2"))
    scale > 0 || throw(ArgumentError("scale must be > 0"))
    nx, ny, nv = size(cube)
    intensity = dropdims(sum(cube; dims=3); dims=3)
    vgt_run(c, I) = begin
        q, u = VGT_QU_mv_2d(c, dn, 1, nv, noise_threshold; min_valid=min_valid, gradient_mode=gradient_mode, kwargs...)
        finalize_vgt_products(q, u, I, ker)
    end
    base = vgt_run(cube, copy(intensity))
    psis = Matrix{Float64}[]
    rels = Matrix{Float64}[]
    for _ in 1:realizations
        noisy = cube .+ correlated_noise(rng, nx, ny, nv; fwhm_major_px, fwhm_minor_px, angle=noise_angle,
                                         sigma=noise_sigma isa Number ? scale * noise_sigma : scale .* noise_sigma,
                                         spectral_lag1)
        r = vgt_run(noisy, dropdims(sum(noisy; dims=3); dims=3))
        push!(psis, r.psi)
        push!(rels, r.reliability)
    end
    spread = circular_spread(psis)
    offset = headless_angle_difference(spread.mean_angle, base.psi)
    sigma_rel = fill(NaN, nx, ny)
    @inbounds for idx in eachindex(sigma_rel)
        vals = [r[idx] for r in rels if isfinite(r[idx])]
        length(vals) >= 2 && (sigma_rel[idx] = std(vals) / scale)
    end
    # circular_spread is the spread about the sample mean of the K maps (divides by K); sqrt(K / (K - 1)) makes the
    # variance unbiased, so error / sigma_psi follows Student's t with K - 1 degrees of freedom when calibrated.
    return (; psi=base.psi, reliability=base.reliability,
            sigma_psi=spread.sigma .* sqrt(realizations / (realizations - 1)) ./ scale, offset,
            sigma_reliability=sigma_rel, n=spread.n)
end
