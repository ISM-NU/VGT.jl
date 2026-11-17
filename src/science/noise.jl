# Noise modeling for VGT thresholding.
#
# VGT gradients from low-SNR emission are unreliable. This file estimates or
# builds noise models so the science kernel can mask pixels below the configured
# noise threshold.

const NEGATIVE_MAD_NORMAL = 0.6744897501960817

struct NoiseModel
    # mode describes which noise variation is represented:
    # :scalar uses one value, :channel varies with velocity, :spatial varies
    # over the sky, and :separable combines spatial and channel factors.
    mode::Symbol
    scalar::Float64
    channel::Union{Nothing, Vector{Float64}}
    spatial::Union{Nothing, Matrix{Float64}}
end

function robust_negative_noise(values; min_count::Int=100)
    # In emission-line cubes, negative intensities are usually noise-dominated.
    # Their median absolute value gives a robust noise estimate without being
    # driven by bright positive emission.
    negative_abs = Float64[]
    for v in values
        if isfinite(v) && v < 0.0
            push!(negative_abs, abs(Float64(v)))
        end
    end
    length(negative_abs) < min_count && return NaN
    return median(negative_abs) / NEGATIVE_MAD_NORMAL
end

function _safe_spread_pct(values)
    clean = sort!(collect(filter(isfinite, values)))
    isempty(clean) && return NaN
    med = median(clean)
    med == 0.0 && return NaN
    return 100.0 * (quantile(clean, 0.84) - quantile(clean, 0.16)) / med
end

function _finite_positive_or(value, fallback::Real)
    return (isfinite(value) && value > 0.0) ? Float64(value) : Float64(fallback)
end

function _estimated_scalar_noise(cube::AbstractArray{<:Any, 3}, fallback::Real)
    return _finite_positive_or(robust_negative_noise(cube), fallback)
end

function estimate_channel_noise(cube::AbstractArray{<:Any, 3}; fallback::Real,
                                min_count::Int=100)
    # Estimate one noise value per velocity channel.
    _, _, nv = size(cube)
    out = fill(Float64(fallback), nv)
    for k in 1:nv
        est = robust_negative_noise(@view cube[:, :, k]; min_count=min_count)
        out[k] = _finite_positive_or(est, fallback)
    end
    return out
end

function estimate_channel_noise(path::AbstractString; fallback::Real,
                                min_count::Int=100, channels=nothing)
    # channels: the input channels to use (velocity_range_kms); nothing = all.
    f = FITS(path)
    try
        ks = channels === nothing ? (1:size(f[1], 3)) : channels
        out = fill(Float64(fallback), length(ks))
        for (i, k) in enumerate(ks)
            sl = read(f[1], :, :, k:k)[:, :, 1]
            out[i] = _finite_positive_or(robust_negative_noise(sl; min_count=min_count), fallback)
        end
        return out
    finally
        close(f)
    end
end

function estimate_spatial_noise(cube::AbstractArray{<:Any, 3}; fallback::Real,
                                tile_size::Int=128, min_count::Int=100)
    # Estimate a noise map on coarse spatial tiles, then assign each tile value
    # to all pixels in that tile.
    nx, ny, _ = size(cube)
    tile_size >= 1 || throw(ArgumentError("noise_tile_size must be >= 1; got $tile_size"))
    out = fill(Float64(fallback), nx, ny)
    for y1 in 1:tile_size:ny, x1 in 1:tile_size:nx
        x2 = min(nx, x1 + tile_size - 1)
        y2 = min(ny, y1 + tile_size - 1)
        est = robust_negative_noise(@view cube[x1:x2, y1:y2, :]; min_count=min_count)
        est = _finite_positive_or(est, fallback)
        out[x1:x2, y1:y2] .= est
    end
    return out
end

function estimate_spatial_noise(path::AbstractString; fallback::Real,
                                tile_size::Int=128, min_count::Int=100, channels=nothing)
    tile_size >= 1 || throw(ArgumentError("noise_tile_size must be >= 1; got $tile_size"))
    f = FITS(path)
    try
        nx, ny, _ = size(f[1])
        out = fill(Float64(fallback), nx, ny)
        for y1 in 1:tile_size:ny, x1 in 1:tile_size:nx
            x2 = min(nx, x1 + tile_size - 1)
            y2 = min(ny, y1 + tile_size - 1)
            tile = read(f[1], x1:x2, y1:y2, channels === nothing ? (:) : channels)
            est = _finite_positive_or(robust_negative_noise(tile; min_count=min_count), fallback)
            out[x1:x2, y1:y2] .= est
        end
        return out
    finally
        close(f)
    end
end

function separable_noise_cube(spatial::AbstractMatrix, channel::AbstractVector, scalar::Real)
    # Combine spatial and channel noise estimates while preserving the scalar
    # scale. This approximates noise(x,y,v) as spatial(x,y)*channel(v)/scalar.
    scalar > 0 || throw(ArgumentError("scalar noise must be > 0"))
    nx, ny = size(spatial)
    nv = length(channel)
    out = Array{Float64}(undef, nx, ny, nv)
    inv_scalar = 1.0 / Float64(scalar)
    @inbounds for k in 1:nv, j in 1:ny, i in 1:nx
        out[i, j, k] = spatial[i, j] * channel[k] * inv_scalar
    end
    return out
end

function _select_auto_noise_mode(channel::AbstractVector, spatial::AbstractMatrix,
                               spread_threshold_pct::Real)
    # Auto mode chooses the simplest noise model that shows meaningful variation.
    # This avoids adding noisy complexity when scalar noise is adequate.
    channel_spread = _safe_spread_pct(channel)
    spatial_spread = _safe_spread_pct(vec(spatial))
    channel_varies = isfinite(channel_spread) && channel_spread > spread_threshold_pct
    spatial_varies = isfinite(spatial_spread) && spatial_spread > spread_threshold_pct
    if channel_varies && spatial_varies
        return :separable
    elseif channel_varies
        return :channel
    elseif spatial_varies
        return :spatial
    end
    return :scalar
end

function _channel_scalar(channel::AbstractVector, fallback::Real)
    clean = collect(filter(v -> isfinite(v) && v > 0.0, channel))
    isempty(clean) && return Float64(fallback)
    return median(clean)
end

function build_noise_model(cube::AbstractArray{<:Any, 3}, scalar_noise::Real,
                           mode::Symbol; tile_size::Int=128,
                           spread_threshold_pct::Real=35.0)
    # Build a NoiseModel from an already-loaded cube, mostly used by tests and
    # small validation scripts.
    scalar = Float64(scalar_noise)
    scalar > 0 || throw(ArgumentError("scalar noise must be > 0"))
    mode === :scalar && return NoiseModel(:scalar, scalar, nothing, nothing)

    if mode === :channel
        channel = estimate_channel_noise(cube; fallback=scalar)
        return NoiseModel(:channel, _channel_scalar(channel, scalar), channel, nothing)
    elseif mode === :spatial
        spatial = estimate_spatial_noise(cube; fallback=scalar, tile_size=tile_size)
        return NoiseModel(:spatial, _estimated_scalar_noise(cube, scalar), nothing, spatial)
    elseif mode === :separable
        channel = estimate_channel_noise(cube; fallback=scalar)
        spatial = estimate_spatial_noise(cube; fallback=scalar, tile_size=tile_size)
        return NoiseModel(:separable, _channel_scalar(channel, scalar), channel, spatial)
    elseif mode === :auto
        channel = estimate_channel_noise(cube; fallback=scalar)
        spatial = estimate_spatial_noise(cube; fallback=scalar, tile_size=tile_size)
        selected = _select_auto_noise_mode(channel, spatial, spread_threshold_pct)
        estimated_scalar = _channel_scalar(channel, scalar)
        if selected === :channel
            return NoiseModel(:channel, estimated_scalar, channel, nothing)
        elseif selected === :spatial
            return NoiseModel(:spatial, estimated_scalar, nothing, spatial)
        elseif selected === :separable
            return NoiseModel(:separable, estimated_scalar, channel, spatial)
        end
        return NoiseModel(:scalar, estimated_scalar, nothing, nothing)
    end

    throw(ArgumentError("Unsupported noise_mode=$mode"))
end

function build_noise_model(path::AbstractString, scalar_noise::Real,
                           mode::Symbol; tile_size::Int=128,
                           spread_threshold_pct::Real=35.0, channels=nothing)
    # Build a NoiseModel directly from a FITS file without requiring the full cube
    # to stay in memory. channels restricts the estimates to the selected input
    # channels (velocity_range_kms), so they match a cube cut to those channels.
    scalar = Float64(scalar_noise)
    scalar > 0 || throw(ArgumentError("scalar noise must be > 0"))
    mode === :scalar && return NoiseModel(:scalar, scalar, nothing, nothing)

    if mode === :channel
        channel = estimate_channel_noise(path; fallback=scalar, channels=channels)
        return NoiseModel(:channel, _channel_scalar(channel, scalar), channel, nothing)
    elseif mode === :spatial
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size, channels=channels)
        return NoiseModel(:spatial, scalar, nothing, spatial)
    elseif mode === :separable
        channel = estimate_channel_noise(path; fallback=scalar, channels=channels)
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size, channels=channels)
        return NoiseModel(:separable, _channel_scalar(channel, scalar), channel, spatial)
    elseif mode === :auto
        channel = estimate_channel_noise(path; fallback=scalar, channels=channels)
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size, channels=channels)
        selected = _select_auto_noise_mode(channel, spatial, spread_threshold_pct)
        estimated_scalar = _channel_scalar(channel, scalar)
        if selected === :channel
            return NoiseModel(:channel, estimated_scalar, channel, nothing)
        elseif selected === :spatial
            return NoiseModel(:spatial, estimated_scalar, nothing, spatial)
        elseif selected === :separable
            return NoiseModel(:separable, estimated_scalar, channel, spatial)
        end
        return NoiseModel(:scalar, estimated_scalar, nothing, nothing)
    end

    throw(ArgumentError("Unsupported noise_mode=$mode"))
end

function noise_for_chunk(model::NoiseModel, nx::Int, y_range::UnitRange{Int}, nv::Int;
                         multiplier::Real=1.0)
    # Expand the selected noise model to match the chunk being processed. The VGT
    # kernel accepts either a scalar noise or a full per-pixel noise cube.
    factor = Float64(multiplier)
    factor > 0 || throw(ArgumentError("noise_sigma_multiplier must be > 0"))
    if model.mode === :scalar
        return model.scalar * factor
    elseif model.mode === :channel
        return (reshape(model.channel, 1, 1, :) .+ zeros(Float64, nx, length(y_range), nv)) .* factor
    elseif model.mode === :spatial
        return repeat(@view(model.spatial[:, y_range]), 1, 1, nv) .* factor
    elseif model.mode === :separable
        spatial = @view model.spatial[:, y_range]
        out = Array{Float64}(undef, nx, length(y_range), nv)
        inv_scalar = factor / model.scalar
        @inbounds for k in 1:nv, j in axes(out, 2), i in 1:nx
            out[i, j, k] = spatial[i, j] * model.channel[k] * inv_scalar
        end
        return out
    end
    throw(ArgumentError("Unsupported noise model mode=$(model.mode)"))
end

function build_noise_cube(cube::AbstractArray{<:Any, 3}, scalar_noise::Real,
                        mode::Symbol; tile_size::Int=128,
                        spread_threshold_pct::Real=35.0,
                        multiplier::Real=1.0)
    model = build_noise_model(
        cube,
        scalar_noise,
        mode;
        tile_size=tile_size,
        spread_threshold_pct=spread_threshold_pct,
    )
    return noise_for_chunk(model, size(cube, 1), 1:size(cube, 2), size(cube, 3);
                           multiplier=multiplier)
end

"""
    measure_gradient_noise(cube; gradient_mode=:central_checked)

Covariance of the gradient vector (gx, gy) over all interior pixels of `cube`
(x, y, channel), with the same stencil and normalisation as the VGT kernel
(`:central_checked`: (a - b) / 2; `:sobel_checked`: Sobel / 8;
`:isotropic_checked`: central difference times [1, 4, 1] / 6). Pass only
emission-free channels. Returns a named tuple: `cov = [sxx, syy, sxy]` (the
`gradient_noise_cov` config value), `sigma_x`, `sigma_y`, the gradient
anisotropy `anisotropy` = |sum g^2 e^(2i theta)| / sum g^2 and the preferred
gradient `direction_deg` (axis of the larger variance, counter-clockwise from
+x), and `n`. Non-finite pixels are skipped.
"""
function measure_gradient_noise(cube::AbstractArray{<:Real, 3};
                                gradient_mode::Symbol=:central_checked)
    gradient_mode in (:central_checked, :sobel_checked, :isotropic_checked) ||
        throw(ArgumentError("measure_gradient_noise supports :central_checked, :sobel_checked and :isotropic_checked; got $gradient_mode"))
    nx, ny, nv = size(cube)
    (nx >= 3 && ny >= 3) || throw(ArgumentError("cube needs at least 3 pixels per sky axis"))
    n = 0
    sx = sy = sxx = syy = sxy = 0.0
    @inbounds for k in 1:nv, j in 2:(ny - 1), i in 2:(nx - 1)
        gx, gy = _stencil_gradient(view(cube, :, :, k), i, j, gradient_mode)
        (isfinite(gx) && isfinite(gy)) || continue
        n += 1
        sx += gx; sy += gy; sxx += gx * gx; syy += gy * gy; sxy += gx * gy
    end
    n >= 100 || throw(ArgumentError("only $n finite gradient samples; need at least 100"))
    mx, my = sx / n, sy / n
    cov = [sxx / n - mx^2, syy / n - my^2, sxy / n - mx * my]
    A, dir = gradient_anisotropy(cov)
    return (; cov, sigma_x=sqrt(cov[1]), sigma_y=sqrt(cov[2]), anisotropy=A, direction_deg=dir, n)
end

"""
    gradient_anisotropy(cov) -> (A, direction_deg)

Anisotropy `A = |sum g^2 e^(2i theta)| / sum g^2 = hypot(sxx - syy, 2 sxy) / (sxx + syy)`
and preferred direction `0.5 atan(2 sxy, sxx - syy)` (degrees) of a gradient
covariance `[sxx, syy, sxy]`.
"""
function gradient_anisotropy(cov::AbstractVector{<:Real})
    sxx, syy, sxy = cov
    return hypot(sxx - syy, 2sxy) / (sxx + syy), rad2deg(0.5 * atan(2sxy, sxx - syy))
end

# Whitening matrices for covariance [sxx, syy, sxy]. With C_n = C / sqrt(det C)
# (unit determinant): S = C_n^(1/2) = (C_n + I) / sqrt(tr C_n + 2) and
# W = C_n^(-1/2) = S^(-1) = [[s22, -s12], [-s12, s11]] since det S = 1.
_gradient_whitening(cov::AbstractVector{<:Real}) = _gradient_whitening(cov[1], cov[2], cov[3])

function _gradient_whitening(sxx::Real, syy::Real, sxy::Real)
    sxx, syy, sxy = Float64(sxx), Float64(syy), Float64(sxy)
    d = sqrt(sxx * syy - sxy^2)
    axx, ayy, axy = sxx / d, syy / d, sxy / d
    norm = sqrt(axx + ayy + 2)
    s11, s22, s12 = (axx + 1) / norm, (ayy + 1) / norm, axy / norm
    return (; w11=s22, w22=s11, w12=-s12, s11, s22, s12)
end

@inline _whiten_gradient(m, gx, gy) = (m.w11 * gx + m.w12 * gy, m.w12 * gx + m.w22 * gy)

# Map a whitened-frame axial angle back to the original frame with S = C_n^(1/2).
@inline function _unwhiten_angle(m, θ)
    c, s = cos(θ), sin(θ)
    return atan(m.s12 * c + m.s22 * s, m.s11 * c + m.s12 * s)
end

# Per-pixel whitening matrices from a (nx, ny, 3) covariance map (sxx, syy, sxy).
function _whitening_maps(covmap::AbstractArray{<:Real, 3})
    size(covmap, 3) == 3 || throw(ArgumentError("covariance map needs 3 planes (sxx, syy, sxy); got $(size(covmap))"))
    nx, ny, _ = size(covmap)
    names = (:w11, :w12, :w22, :s11, :s12, :s22)
    maps = NamedTuple{names}(ntuple(_ -> Matrix{Float64}(undef, nx, ny), 6))
    @inbounds for j in 1:ny, i in 1:nx
        sxx, syy, sxy = covmap[i, j, 1], covmap[i, j, 2], covmap[i, j, 3]
        # No covariance (NaN, e.g. outside the survey coverage, where the map
        # builder had no emission-free data): the whitened gradients there are
        # NaN, so those pixels drop out instead of being used unwhitened.
        if !(isfinite(sxx) && isfinite(syy) && isfinite(sxy))
            for n in names
                maps[n][i, j] = NaN
            end
            continue
        end
        (sxx > 0 && syy > 0 && sxx * syy - sxy^2 > 0) ||
            throw(ArgumentError("covariance map must be positive definite wherever it is finite; failed at ($i, $j): [$sxx, $syy, $sxy]"))
        m = _gradient_whitening(sxx, syy, sxy)
        for n in names
            maps[n][i, j] = getfield(m, n)
        end
    end
    return maps
end

@inline function _whiten_gradient_at(maps, i, j, gx, gy)
    w11, w12, w22 = maps.w11[i, j], maps.w12[i, j], maps.w22[i, j]
    return (w11 * gx + w12 * gy, w12 * gx + w22 * gy)
end

@inline function _unwhiten_angle_at(maps, i, j, θ)
    c, s = cos(θ), sin(θ)
    s11, s12, s22 = maps.s11[i, j], maps.s12[i, j], maps.s22[i, j]
    return atan(s12 * c + s22 * s, s11 * c + s12 * s)
end

# NaN-aware (normalised) separable Gaussian smoothing with std `σ` pixels.
function _nan_gaussian_smooth(a::AbstractMatrix{<:Real}, σ::Real)
    σ > 0 || throw(ArgumentError("smoothing sigma must be > 0"))
    r = max(1, ceil(Int, 4σ))
    kern = [exp(-0.5 * (x / σ)^2) for x in -r:r]
    nx, ny = size(a)
    num = [isfinite(v) ? Float64(v) : 0.0 for v in a]
    den = [isfinite(v) ? 1.0 : 0.0 for v in a]
    function pass(x, dim)
        out = zeros(nx, ny)
        @inbounds if dim == 1
            for j in 1:ny, i in 1:nx
                acc = 0.0
                for t in max(1, i - r):min(nx, i + r)
                    acc += kern[t - i + r + 1] * x[t, j]
                end
                out[i, j] = acc
            end
        else
            for j in 1:ny, i in 1:nx
                acc = 0.0
                for t in max(1, j - r):min(ny, j + r)
                    acc += kern[t - j + r + 1] * x[i, t]
                end
                out[i, j] = acc
            end
        end
        return out
    end
    num = pass(pass(num, 1), 2)
    den = pass(pass(den, 1), 2)
    return [d > 1e-6 ? n / d : NaN for (n, d) in zip(num, den)]
end

"""
    gradient_noise_cov_map(cube; gradient_mode=:central_checked, fwhm_pixels=60) -> Array{Float64,3}

Spatially varying gradient noise covariance for `gradient_noise_cov_file`: planes
(sxx, syy, sxy) of shape (nx, ny, 3). Per pixel, the products gx^2, gy^2 and gx*gy
of the gradients of the emission-free channels in `cube` (stencil of
`gradient_mode`) are averaged over the channels, then smoothed with a NaN-aware
Gaussian of the given FWHM in pixels (use about 3 x dn). Non-finite gradients are
skipped. A larger FWHM gives a smoother, less noisy map; the noise anisotropy
must vary slowly compared with it.
"""
function gradient_noise_cov_map(cube::AbstractArray{<:Real, 3};
                                gradient_mode::Symbol=:central_checked,
                                fwhm_pixels::Real=60)
    gradient_mode in (:central_checked, :sobel_checked, :isotropic_checked) ||
        throw(ArgumentError("gradient_noise_cov_map supports :central_checked, :sobel_checked and :isotropic_checked; got $gradient_mode"))
    fwhm_pixels > 0 || throw(ArgumentError("fwhm_pixels must be > 0"))
    nx, ny, nv = size(cube)
    (nx >= 3 && ny >= 3) || throw(ArgumentError("cube needs at least 3 pixels per sky axis"))
    sums = zeros(nx, ny, 3)
    counts = zeros(nx, ny)
    @inbounds for k in 1:nv, j in 2:(ny - 1), i in 2:(nx - 1)
        gx, gy = _stencil_gradient(view(cube, :, :, k), i, j, gradient_mode)
        (isfinite(gx) && isfinite(gy)) || continue
        sums[i, j, 1] += gx * gx; sums[i, j, 2] += gy * gy; sums[i, j, 3] += gx * gy
        counts[i, j] += 1
    end
    σ = fwhm_pixels / (2 * sqrt(2 * log(2)))
    out = Array{Float64}(undef, nx, ny, 3)
    for p in 1:3
        plane = [counts[i, j] > 0 ? sums[i, j, p] / counts[i, j] : NaN for i in 1:nx, j in 1:ny]
        out[:, :, p] = _nan_gaussian_smooth(plane, σ)
    end
    return out
end

read_gradient_noise_cov_map(path::AbstractString) = Float64.(read_fits_image(path))

"Write a `gradient_noise_cov_map` result as a 3-plane FITS cube (planes: sxx, syy, sxy)."
function write_gradient_noise_cov_map(path::AbstractString, covmap::AbstractArray{<:Real, 3};
                                      gradient_mode::Symbol=:central_checked, fwhm_pixels::Real=60)
    header = FITSHeader(["PLANE1", "PLANE2", "PLANE3", "GRADMODE", "SMOOTHFW"],
                        ["sxx", "syy", "sxy", String(gradient_mode), Float64(fwhm_pixels)],
                        ["gx*gx", "gy*gy", "gx*gy", "gradient stencil", "smoothing FWHM in pixels"])
    return write_fits(path, covmap; header=header)
end
