# Final science products from accumulated VGT Q/U/I maps.
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
    # Collapse the velocity axis to make the intensity normalization map.
    return sum(cube; dims=3)[:, :, 1]
end

function compute_reliability(q_smooth::AbstractMatrix, u_smooth::AbstractMatrix, intensity_smooth::AbstractMatrix)
    # Reliability is the fractional polarized/coherent amplitude of the
    # smoothed orientation signal: sqrt(Q^2+U^2)/I.
    reliability = Array{Float64}(undef, size(q_smooth))
    @inbounds for i in eachindex(reliability)
        intensity_value = intensity_smooth[i]
        reliability[i] = (isfinite(intensity_value) && intensity_value > 0.0) ?
            sqrt(q_smooth[i]^2 + u_smooth[i]^2) / intensity_value :
            NaN
    end
    return reliability
end

@inline function _invalid_finalize_intensity(value)
    # The mask rejects non-finite and non-positive intensity, because
    # reliability is normalized by intensity.
    return !isfinite(value) || value <= 0
end

function finalize_vgt_products(q_accum::AbstractMatrix, u_accum::AbstractMatrix, intensity_accum::AbstractMatrix,
                               ker::Int,
                               zero_stokes_policy::Symbol=:nan_before_smoothing,
                               output_support_mode::Symbol=:smoothed_stokes_support)
    # Convert raw accumulated Q/U/I maps into the final psi and reliability
    # products. The policy arguments should not be changed casually, because they alter output support and masks.
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))
    q_before_smoothing = copy(q_accum)
    u_before_smoothing = copy(u_accum)
    intensity_before_smoothing = copy(intensity_accum)
    if zero_stokes_policy === :nan_before_smoothing
        @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing, q_accum, u_accum)
            q_accum[idx] == 0 && (q_before_smoothing[idx] = NaN)
            u_accum[idx] == 0 && (u_before_smoothing[idx] = NaN)
        end
    end
    @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing,
                                   intensity_before_smoothing, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_before_smoothing[idx] = NaN
            u_before_smoothing[idx] = NaN
            intensity_before_smoothing[idx] = NaN
        end
    end
    q_smooth = imfilter_gaussian(q_before_smoothing, ker)
    u_smooth = imfilter_gaussian(u_before_smoothing, ker)
    intensity_smooth = imfilter_gaussian(intensity_before_smoothing, ker)
    @inbounds for idx in eachindex(q_smooth, u_smooth, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_smooth[idx] = NaN
            u_smooth[idx] = NaN
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
    reliability = compute_reliability(q_smooth, u_smooth, intensity_smooth)
    return (; q_smooth, u_smooth, intensity_smooth, psi, reliability)
end

function headless_angle_difference(a, b)
    # Difference between two axial angles, wrapped so 0 and pi represent the
    # same orientation.
    return mod.(a .- b .+ pi / 2, pi) .- pi / 2
end
