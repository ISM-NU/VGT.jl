# Pure helpers for memory-bounded FITS cube regridding.

"""Return the FITS pixel coordinate corresponding to a linear-axis value."""
@inline function source_axis_position(world_value::Real, crval::Real, crpix::Real, cdelt::Real)
    cdelt == 0 && throw(ArgumentError("Linear WCS axis has zero CDELT"))
    return (Float64(world_value) - Float64(crval)) / Float64(cdelt) + Float64(crpix)
end

"""Return lower/upper 1-based pixel indices and interpolation weight."""
function linear_bracket(position::Real, n::Int)
    n >= 1 || throw(ArgumentError("Axis length must be positive; got $n"))
    p = Float64(position)
    (isfinite(p) && 1.0 <= p <= n) || return nothing
    lo = clamp(floor(Int, p), 1, n)
    hi = min(lo + 1, n)
    return (lo=lo, hi=hi, w=hi == lo ? 0.0 : p - lo)
end

"""Sample a 2-D image at a 1-based fractional pixel coordinate."""
function bilinear_sample(image::AbstractMatrix, x::Real, y::Real)
    bx = linear_bracket(x, size(image, 1))
    by = linear_bracket(y, size(image, 2))
    (bx === nothing || by === nothing) && return Float32(NaN)

    v11 = image[bx.lo, by.lo]
    v21 = image[bx.hi, by.lo]
    v12 = image[bx.lo, by.hi]
    v22 = image[bx.hi, by.hi]
    (isfinite(v11) && isfinite(v21) && isfinite(v12) && isfinite(v22)) ||
        return Float32(NaN)

    tx = bx.w
    ty = by.w
    value = (1 - tx) * (1 - ty) * Float64(v11) +
            tx * (1 - ty) * Float64(v21) +
            (1 - tx) * ty * Float64(v12) +
            tx * ty * Float64(v22)
    return Float32(value)
end

"""Return target channels whose centers lie within the source linear axis."""
function velocity_overlap_indices(target_header, target_nv::Int, source_header, source_nv::Int)
    target_nv >= 1 || throw(ArgumentError("target_nv must be positive"))
    source_nv >= 1 || throw(ArgumentError("source_nv must be positive"))
    target_values = [
        Float64(target_header["CRVAL3"]) +
        (k - Float64(target_header["CRPIX3"])) * Float64(target_header["CDELT3"])
        for k in 1:target_nv
    ]
    source_values = [
        Float64(source_header["CRVAL3"]) +
        (k - Float64(source_header["CRPIX3"])) * Float64(source_header["CDELT3"])
        for k in 1:source_nv
    ]
    lo, hi = extrema(source_values)
    selected = findall(v -> lo <= v <= hi, target_values)
    isempty(selected) && throw(ArgumentError("Source and target velocity axes do not overlap"))
    return first(selected):last(selected)
end
