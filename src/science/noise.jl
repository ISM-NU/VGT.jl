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
                                min_count::Int=100)
    f = FITS(path)
    try
        _, _, nv = size(f[1])
        out = fill(Float64(fallback), nv)
        for k in 1:nv
            sl = read(f[1], :, :, k:k)[:, :, 1]
            out[k] = _finite_positive_or(robust_negative_noise(sl; min_count=min_count), fallback)
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
                                tile_size::Int=128, min_count::Int=100)
    tile_size >= 1 || throw(ArgumentError("noise_tile_size must be >= 1; got $tile_size"))
    f = FITS(path)
    try
        nx, ny, _ = size(f[1])
        out = fill(Float64(fallback), nx, ny)
        for y1 in 1:tile_size:ny, x1 in 1:tile_size:nx
            x2 = min(nx, x1 + tile_size - 1)
            y2 = min(ny, y1 + tile_size - 1)
            tile = read(f[1], x1:x2, y1:y2, :)
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
                           spread_threshold_pct::Real=35.0)
    # Build a NoiseModel directly from a FITS file without requiring the full cube
    # to stay in memory.
    scalar = Float64(scalar_noise)
    scalar > 0 || throw(ArgumentError("scalar noise must be > 0"))
    mode === :scalar && return NoiseModel(:scalar, scalar, nothing, nothing)

    if mode === :channel
        channel = estimate_channel_noise(path; fallback=scalar)
        return NoiseModel(:channel, _channel_scalar(channel, scalar), channel, nothing)
    elseif mode === :spatial
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size)
        return NoiseModel(:spatial, scalar, nothing, spatial)
    elseif mode === :separable
        channel = estimate_channel_noise(path; fallback=scalar)
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size)
        return NoiseModel(:separable, _channel_scalar(channel, scalar), channel, spatial)
    elseif mode === :auto
        channel = estimate_channel_noise(path; fallback=scalar)
        spatial = estimate_spatial_noise(path; fallback=scalar, tile_size=tile_size)
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
