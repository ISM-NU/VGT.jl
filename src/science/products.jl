# Final science products from accumulated VGT Q/U/W maps (I masks sky).
#
# The VGT kernel accumulates gradient orientations. This file smooths those
# Stokes-like maps, rotates gradients into magnetic-field orientation proxies,
# and computes reliability/coherence diagnostics.

function imfilter_gaussian(d::AbstractMatrix, p::Number)
    # Gaussian smoothing pools nearby orientation evidence before converting
    # Q/U into the final angle map.
    return imfilter(d, Kernel.gaussian(p), NA())
end

function compute_stokes_intensity(cube)
    # Collapse the velocity axis. finalize_vgt_products uses this sum only to
    # mask sky; the reliability denominator is the fitted Stokes weight.
    return sum(cube; dims=3)[:, :, 1]
end

function compute_reliability(q_smooth::AbstractMatrix, u_smooth::AbstractMatrix, weight_smooth::AbstractMatrix)
    # Reliability is the fractional coherent amplitude of the smoothed
    # orientation signal: sqrt(Q^2+U^2)/W. It lies in [0, 1] (up to rounding)
    # only if W is the sum of the same non-negative weights as Q/U and all
    # three maps were smoothed over the same pixels; finalize_vgt_products
    # ensures both.
    reliability = Array{Float64}(undef, size(q_smooth))
    @inbounds for i in eachindex(reliability)
        weight_value = weight_smooth[i]
        reliability[i] = (isfinite(weight_value) && weight_value > 0.0) ?
            sqrt(q_smooth[i]^2 + u_smooth[i]^2) / weight_value :
            NaN
    end
    return reliability
end

@inline function _invalid_finalize_intensity(value)
    # Sky mask: non-finite or non-positive channel-summed intensity gets no
    # psi. This sets the psi support; it no longer normalises reliability.
    return !isfinite(value) || value <= 0
end

"""
    finalize_vgt_products(q_accum, u_accum, intensity_accum, ker,
                          zero_stokes_policy=:nan_before_smoothing,
                          output_support_mode=:smoothed_stokes_support;
                          fitted_weight=nothing)

Smooth the accumulated Q/U maps and form `psi` and the reliability map.

`intensity_accum` (the channel-summed intensity) only masks sky: pixels where it
is non-finite or non-positive get no psi. `fitted_weight` is the sum of the
Stokes weights that went into Q/U (`profile[:fitted_weight]` from
`VGT_QU_mv_2d_profiled`, or the stitched chunk map). It is smoothed over the
same pixels as Q/U and is the reliability denominator, so
`reliability = sqrt(Q^2+U^2)/W` is in [0, 1] up to rounding. Without
`fitted_weight` the reliability map is all NaN; psi does not depend on it.

Returns `(; q_smooth, u_smooth, psi, reliability, q_reliability, u_reliability,
weight_smooth)`. `q_reliability`/`u_reliability` are the Q/U behind the
reliability map: the same arrays as `q_smooth`/`u_smooth` unless Q or U is
exactly zero at a pixel with fitted weight (see the comment in the body).
"""
function finalize_vgt_products(q_accum::AbstractMatrix, u_accum::AbstractMatrix, intensity_accum::AbstractMatrix,
                               ker::Int,
                               zero_stokes_policy::Symbol=:nan_before_smoothing,
                               output_support_mode::Symbol=:smoothed_stokes_support;
                               edge_mask::Bool=false,
                               fitted_weight::Union{Nothing, AbstractMatrix}=nothing)
    # Convert raw accumulated Q/U maps into the final psi and reliability
    # products. The policy arguments should not be changed casually, because they alter output support and masks.
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))
    fitted_weight === nothing || size(fitted_weight) == size(q_accum) ||
        throw(ArgumentError("fitted_weight size $(size(fitted_weight)) must match Q size $(size(q_accum))"))
    q_before_smoothing = copy(q_accum)
    u_before_smoothing = copy(u_accum)
    if zero_stokes_policy === :nan_before_smoothing
        @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing, q_accum, u_accum)
            # A pixel with no fit has Q == U == 0; mask both together. Masking
            # each map on its own zeros dropped every pixel of a field whose
            # fitted angles are exactly 0 or 90 deg (U == 0 or Q == 0 there),
            # leaving psi all NaN for a perfectly axis-aligned field.
            if q_accum[idx] == 0 && u_accum[idx] == 0
                q_before_smoothing[idx] = NaN
                u_before_smoothing[idx] = NaN
            end
        end
    end
    @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_before_smoothing[idx] = NaN
            u_before_smoothing[idx] = NaN
        end
    end
    q_smooth = imfilter_gaussian(q_before_smoothing, ker)
    u_smooth = imfilter_gaussian(u_before_smoothing, ker)
    @inbounds for idx in eachindex(q_smooth, u_smooth, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_smooth[idx] = NaN
            u_smooth[idx] = NaN
        end
    end
    if edge_mask
        # The kernel marked pixels whose window footprint reaches a blank or the
        # array border with NaN Q/U. Smoothing would fill them from neighbours,
        # so they stay NaN: no direction there.
        @inbounds for idx in eachindex(q_smooth, u_smooth, q_accum, u_accum)
            if !isfinite(q_accum[idx]) || !isfinite(u_accum[idx])
                q_smooth[idx] = NaN
                u_smooth[idx] = NaN
            end
        end
    end
    if output_support_mode === :direct_fit_support
        @inbounds for idx in eachindex(q_smooth, u_smooth, q_accum, u_accum)
            if q_accum[idx] == 0.0 && u_accum[idx] == 0.0
                q_smooth[idx] = NaN
                u_smooth[idx] = NaN
            end
        end
    end
    # atan(U,Q)/2 recovers the axial gradient orientation. VGT uses gradients as
    # a magnetic-field proxy after a 90-degree rotation, hence +pi/2.
    psi = 0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)

    if fitted_weight === nothing
        reliability = fill(NaN, size(psi))
        return (; q_smooth, u_smooth, psi, reliability, q_reliability=q_smooth, u_reliability=u_smooth,
                weight_smooth=fill(NaN, size(psi)))
    end
    # The normalised (NaN-aware) filter divides each map by the kernel mass of
    # its own finite pixels, so Q, U and W must be NaN at exactly the same
    # pixels for P <= 1 to hold after smoothing. W follows the zero policy:
    # under :nan_before_smoothing pixels without fitted weight are NaN.
    w_before_smoothing = copy(fitted_weight)
    @inbounds for idx in eachindex(w_before_smoothing, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx]) ||
           (zero_stokes_policy === :nan_before_smoothing && fitted_weight[idx] == 0)
            w_before_smoothing[idx] = NaN
        end
    end
    weight_smooth = imfilter_gaussian(w_before_smoothing, ker)
    # With weights >= 0, W == 0 implies Q == U == 0, so in practice the masks
    # already agree and the psi Q/U are reused. They differ where Q or U is
    # exactly zero at a fitted pixel (e.g. a fitted angle of exactly 0 gives
    # U == 0); then Q/U are smoothed again on W's mask for the reliability only,
    # and psi keeps its historical masks.
    same_mask = all(eachindex(w_before_smoothing)) do idx
        w_nan = isnan(w_before_smoothing[idx])
        isnan(q_before_smoothing[idx]) == w_nan && isnan(u_before_smoothing[idx]) == w_nan
    end
    q_reliability, u_reliability = if same_mask
        q_smooth, u_smooth
    else
        q_r = copy(q_accum)
        u_r = copy(u_accum)
        @inbounds for idx in eachindex(q_r, u_r, w_before_smoothing)
            if isnan(w_before_smoothing[idx])
                q_r[idx] = NaN
                u_r[idx] = NaN
            end
        end
        q_r = imfilter_gaussian(q_r, ker)
        u_r = imfilter_gaussian(u_r, ker)
        # Reliability is reported only where psi is.
        @inbounds for idx in eachindex(q_r, u_r, psi)
            if isnan(psi[idx])
                q_r[idx] = NaN
                u_r[idx] = NaN
            end
        end
        q_r, u_r
    end
    reliability = compute_reliability(q_reliability, u_reliability, weight_smooth)
    return (; q_smooth, u_smooth, psi, reliability, q_reliability, u_reliability, weight_smooth)
end

function headless_angle_difference(a, b)
    # Difference between two axial angles, wrapped so 0 and pi represent the
    # same orientation.
    return mod.(a .- b .+ pi / 2, pi) .- pi / 2
end

"""
    block_mean(a, b; valid=isfinite.(a)) -> Matrix{Float64}

Mean of `a` over non-overlapping `b` x `b` blocks, using only entries where
`valid` is true; a block with no valid entry is NaN. Partial blocks at the
right and top edges are dropped, so output pixel (i, j) covers input pixels
(i-1)b+1:ib by (j-1)b+1:jb and the binned WCS stays exact.
"""
function block_mean(a::AbstractMatrix, b::Int; valid::AbstractMatrix{Bool}=isfinite.(a))
    b >= 1 || throw(ArgumentError("block size must be positive; got $b"))
    nx, ny = size(a) .÷ b
    out = Matrix{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        total = 0.0
        count = 0
        for jj in (j - 1) * b + 1:j * b, ii in (i - 1) * b + 1:i * b
            valid[ii, jj] || continue
            total += a[ii, jj]
            count += 1
        end
        out[i, j] = count == 0 ? NaN : total / count
    end
    return out
end

"""
    bin_header(header, b) -> header

Copy of a cube or map header for a grid of `b` x `b` blocks of axes 1 and 2:
`NAXIS1/2`, `CRPIX1/2`, `CDELT1/2` and `CDi_j` (j = 1, 2) are rescaled so each
binned pixel centre has the world coordinates of its block centre. `PC` and
the reference value are unchanged. Returns `nothing` for `nothing`.
"""
bin_header(::Nothing, ::Int) = nothing
function bin_header(header, b::Int)
    b == 1 && return header
    out = deepcopy(header)
    for axis in 1:2
        n = "NAXIS$axis"
        haskey(out, n) && (out[n] = Int(out[n]) ÷ b)
        # Input pixel x lies at binned coordinate (x - 0.5) / b + 0.5.
        c = "CRPIX$axis"
        haskey(out, c) && (out[c] = (Float64(out[c]) - 0.5) / b + 0.5)
        d = "CDELT$axis"
        haskey(out, d) && (out[d] = Float64(out[d]) * b)
        for i in 1:2
            cd = "CD$(i)_$(axis)"
            haskey(out, cd) && (out[cd] = Float64(out[cd]) * b)
        end
    end
    return out
end

"""
    bin_products(finalized, fit_rbar2, fit_support, header, b)

Products on the `b` x `b` binned grid: Q and U are averaged over the pixels
where both are valid and psi is recomputed from the averages (never by
averaging angles). The reliability Q/U and fitted weight W are averaged over
the pixels where all three are valid, so the binned reliability stays in
[0, 1]. The fit diagnostics are plain block means. Returns
`(; psi, reliability, fit_rbar2, fit_support, header)`.
"""
function bin_products(finalized, fit_rbar2::AbstractMatrix, fit_support::AbstractMatrix, header, b::Int)
    valid = isfinite.(finalized.q_smooth) .& isfinite.(finalized.u_smooth)
    q = block_mean(finalized.q_smooth, b; valid)
    u = block_mean(finalized.u_smooth, b; valid)
    psi = 0.5 .* atan.(u, q) .+ (pi / 2)
    valid_r = isfinite.(finalized.q_reliability) .& isfinite.(finalized.u_reliability) .&
              isfinite.(finalized.weight_smooth)
    reliability = compute_reliability(block_mean(finalized.q_reliability, b; valid=valid_r),
                                      block_mean(finalized.u_reliability, b; valid=valid_r),
                                      block_mean(finalized.weight_smooth, b; valid=valid_r))
    return (; psi, reliability,
            fit_rbar2=block_mean(fit_rbar2, b), fit_support=block_mean(fit_support, b),
            header=bin_header(header, b))
end
