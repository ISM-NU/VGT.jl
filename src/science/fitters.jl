# Local orientation fitters for VGT windows.
#
# A window contains many gradient vectors (Gx, Gy). This file converts those
# vectors into a single headless orientation plus diagnostics. The production
# path is the unbinned axial doubled-angle estimator.

Base.@kwdef struct WindowFitResult
    # angle is the fitted gradient orientation in radians. The remaining fields
    # describe how concentrated or uncertain the local angle distribution was.
    angle::Float64 = NaN
    rbar2::Float64 = NaN
    effective_n::Int = 0
    low_confidence::Bool = false
    fit_ok::Bool = false
    failure_reason::Symbol = :not_run
end

@inline function _axial_gradient_angle(ax, ay)
    # atan(y, x) keeps quadrant information. The result is wrapped into
    # [-pi/2, pi/2), because orientations are headless rather than arrows.
    isfinite(ax) && isfinite(ay) || return NaN
    angle = _wrap_axial_angle(atan(Float64(ay), Float64(ax)))
    return isfinite(angle) ? angle : NaN
end

function _collect_gradient_angles(grad_x::AbstractArray, grad_y::AbstractArray)
    # Convert finite gradient pairs into window angles. Gx and Gy from the same
    # pixel are kept together; a pair with any non-finite component is dropped.
    gradient_angles = Float64[]
    sizehint!(gradient_angles, length(grad_x))
    @inbounds for idx in eachindex(grad_x, grad_y)
        angle = _axial_gradient_angle(grad_x[idx], grad_y[idx])
        isnan(angle) || push!(gradient_angles, angle)
    end
    return gradient_angles
end

function _wrap_axial_angle(angle::Float64)
    # Wrap a direction onto the axial interval. Angles separated by pi are
    # physically equivalent for a plane-of-sky field orientation.
    return mod(angle + pi / 2, pi) - pi / 2
end

function _axial_circular_mean(gradient_angles::AbstractVector{<:Real})
    # Doubled-angle averaging handles headless orientations correctly. Averaging
    # raw angles directly would fail near the -pi/2 and +pi/2 boundary.
    s = 0.0
    c = 0.0
    @inbounds for angle in gradient_angles
        alpha = 2.0 * Float64(angle)
        s += sin(alpha)
        c += cos(alpha)
    end
    return 0.5 * atan(s, c)
end

function _axial_vonmises_angle(gradient_angles::AbstractVector{<:Real})
    return _axial_resultant_stats(gradient_angles).angle
end

function _axial_resultant_stats(gradient_angles::AbstractVector{<:Real})
    # rbar2 close to 1 means angles in the window mostly agree. Values near 0
    # mean the window has weak or conflicting orientation evidence.
    n = length(gradient_angles)
    n == 0 && return (; angle=NaN, rbar2=NaN, effective_n=0, resultant=NaN)

    s = 0.0
    c = 0.0
    @inbounds for angle in gradient_angles
        alpha = 2.0 * Float64(angle)
        s += sin(alpha)
        c += cos(alpha)
    end

    resultant = hypot(s, c)
    rbar2 = resultant / n
    angle = resultant <= sqrt(eps(Float64)) ? NaN : _wrap_axial_angle(0.5 * atan(s, c))
    angle = isfinite(angle) ? angle : NaN
    return (; angle, rbar2, effective_n=n, resultant)
end

function _fit_result(angle::Real, stats; failure_reason::Symbol=:fit_failed)
    fit_ok = isfinite(angle)
    low_confidence = isfinite(stats.rbar2) && stats.rbar2 < FIT_LOW_CONFIDENCE_RBAR2_THRESHOLD
    return WindowFitResult(
        angle=Float64(angle),
        rbar2=Float64(stats.rbar2),
        effective_n=stats.effective_n,
        low_confidence=low_confidence,
        fit_ok=fit_ok,
        failure_reason=fit_ok ? :none : failure_reason,
    )
end

function _failed_fit_result(reason::Symbol, stats=(; rbar2=NaN, effective_n=0))
    return WindowFitResult(
        rbar2=Float64(stats.rbar2),
        effective_n=stats.effective_n,
        low_confidence=isfinite(stats.rbar2) && stats.rbar2 < FIT_LOW_CONFIDENCE_RBAR2_THRESHOLD,
        failure_reason=reason,
    )
end

function fit_axial_vonmises(grad_x::AbstractArray, grad_y::AbstractArray)
    # Production-style unbinned axial estimator. "vonmises" here reflects the
    # circular-statistics model and diagnostics.
    gradient_angles = _collect_gradient_angles(grad_x, grad_y)
    isempty(gradient_angles) && return NaN, NaN

    angle = _axial_vonmises_angle(gradient_angles)
    return isfinite(angle) ? (angle, NaN) : (NaN, NaN)
end

# Same sums as _collect_gradient_angles + _axial_resultant_stats for the
# production path (paired finite samples), in one pass and without the
# temporary angle vector. Bitwise identical; it avoids millions of small
# allocations per run, which stall every thread during garbage collection.
function _axial_resultant_stats_paired(grad_x::AbstractArray, grad_y::AbstractArray)
    s = 0.0
    c = 0.0
    n = 0
    @inbounds for idx in eachindex(grad_x, grad_y)
        angle = _axial_gradient_angle(grad_x[idx], grad_y[idx])
        isnan(angle) && continue
        alpha = 2.0 * Float64(angle)
        s += sin(alpha)
        c += cos(alpha)
        n += 1
    end
    n == 0 && return (; angle=NaN, rbar2=NaN, effective_n=0, resultant=NaN)
    resultant = hypot(s, c)
    rbar2 = resultant / n
    angle = resultant <= sqrt(eps(Float64)) ? NaN : _wrap_axial_angle(0.5 * atan(s, c))
    angle = isfinite(angle) ? angle : NaN
    return (; angle, rbar2, effective_n=n, resultant)
end

# Gradient-amplitude (inverse-variance) weighting, opt-in via gradient_weight =
# :amplitude2. Each pair's doubled angle is weighted by |g|^2, so the sums are those
# of the structure tensor: sum (gx^2 - gy^2) and sum 2 gx gy. With isotropic gradient
# noise of variance s^2 per component, a gradient's doubled-angle variance is about
# 2 s^2 / |G|^2 where it is signal dominated, so |g|^2 is its inverse-variance weight,
# and the noise adds 2 s^2 to E|g|^2 but nothing to E(gx^2 - gy^2) or E(gx gy): the
# angle is unbiased. Anisotropic noise does bias it (by the full covariance anisotropy,
# twice the unit-weight pull), so combine with gradient_whiten there. rbar2 is the
# weighted resultant length |sum w e^(2i theta)| / sum w; effective_n still counts pairs.
function _axial_resultant_stats_paired_amp2(grad_x::AbstractArray, grad_y::AbstractArray)
    s = 0.0
    c = 0.0
    wsum = 0.0
    n = 0
    @inbounds for idx in eachindex(grad_x, grad_y)
        gx = Float64(grad_x[idx])
        gy = Float64(grad_y[idx])
        (isfinite(gx) && isfinite(gy)) || continue
        c += gx * gx - gy * gy
        s += 2.0 * gx * gy
        wsum += gx * gx + gy * gy
        n += 1
    end
    n == 0 && return (; angle=NaN, rbar2=NaN, effective_n=0, resultant=NaN)
    resultant = hypot(s, c)
    rbar2 = wsum > 0 ? resultant / wsum : NaN
    # Scale-free: a fixed floor would drop every window when |g| is small (units of 1e-6).
    angle = resultant <= sqrt(eps(Float64)) * wsum ? NaN : _wrap_axial_angle(0.5 * atan(s, c))
    angle = isfinite(angle) ? angle : NaN
    # resultant is reported on the unit-weight scale (rbar2 * n) so the degenerate test below matches.
    return (; angle, rbar2, effective_n=n, resultant=isfinite(rbar2) ? rbar2 * n : 0.0)
end

function fit_window_result(mode::Symbol, grad_x::AbstractArray, grad_y::AbstractArray;
                           gradient_weight::Symbol=:unit)
    # Central dispatcher used by the VGT kernel. It returns one consistent
    # WindowFitResult for the (single) supported fitter mode.
    mode === :axial_vonmises || throw(ArgumentError("Unsupported fitter_mode=$mode"))
    gradient_weight in GRADIENT_WEIGHT_MODES ||
        throw(ArgumentError("Unsupported gradient_weight=$gradient_weight; use one of $(GRADIENT_WEIGHT_MODES)"))
    stats = gradient_weight === :unit ? _axial_resultant_stats_paired(grad_x, grad_y) :
            _axial_resultant_stats_paired_amp2(grad_x, grad_y)
    stats.effective_n == 0 && return _failed_fit_result(:empty_sample)
    reason = stats.resultant <= sqrt(eps(Float64)) ? :degenerate_resultant : :fit_failed
    return _fit_result(stats.angle, stats; failure_reason=reason)
end

function fit_window_angle(mode::Symbol, grad_x::AbstractArray, grad_y::AbstractArray)
    result = fit_window_result(mode, grad_x, grad_y)
    return result.angle, NaN
end
