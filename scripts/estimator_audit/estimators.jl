# Candidate per-window orientation estimators for the estimator audit.
#
# Every estimator takes the gradient window views (gx, gy) that the production
# kernel passes to fit_window_result and returns a VGT.WindowFitResult. The
# rbar2 field carries each estimator's own per-window confidence, so the
# psi_fit_rbar2 diagnostic map can be compared between estimators:
#   :vm      mean resultant length R of unit doubled-angle vectors (production)
#   :st1/2   |sum w e^{2i theta}| / sum w with w = |g| or |g|^2
#   :ms*     R of the final kernel-weighted iteration
#   :vmu     fitted peak fraction of the von Mises + uniform mixture
#   :gauss   R (the Gaussian baseline reports R like the old repository did)
# Hyperparameters were fixed before any run: mean-shift kappa 2 and 4 on the
# doubled angle (half weight at 24.6 and 17.0 deg from the current mode), three
# iterations started from the plain mean; mixture: 8 EM iterations, kappa <= 50.

using VGT: WindowFitResult

const TINY = sqrt(eps(Float64))

@inline _wrap(a) = mod(a + pi / 2, pi) - pi / 2

@inline function _result(s, c, wsum, n, conf)
    r = hypot(s, c)
    (n == 0 || wsum <= 0) && return WindowFitResult(failure_reason=:empty_sample)
    r <= TINY * wsum && return WindowFitResult(rbar2=0.0, effective_n=n, failure_reason=:degenerate_resultant)
    a = _wrap(0.5 * atan(s, c))
    return WindowFitResult(angle=a, rbar2=conf, effective_n=n, low_confidence=conf < 0.5,
                           fit_ok=isfinite(a), failure_reason=isfinite(a) ? :none : :fit_failed)
end

# Unit doubled-angle vector of one gradient without trig:
# cos 2t = (gx^2 - gy^2)/|g|^2, sin 2t = 2 gx gy/|g|^2. Exact zero gradients are
# skipped here (production keeps them as angle 0, see the audit report).
@inline function _unit2(x, y)
    m2 = x * x + y * y
    return (x * x - y * y) / m2, 2 * x * y / m2, m2
end

@inline _ok(x, y) = isfinite(x) && isfinite(y) && (x != 0 || y != 0)

# Weighted doubled-angle mean. p = 0: unit vectors (production up to the zero
# gradient rule), p = 1: |g| weights, p = 2: |g|^2 weights (structure tensor).
function est_weighted(gx, gy, p::Int)
    s = 0.0; c = 0.0; w = 0.0; n = 0
    @inbounds for k in eachindex(gx, gy)
        x = Float64(gx[k]); y = Float64(gy[k])
        _ok(x, y) || continue
        if p == 2
            c += x * x - y * y; s += 2 * x * y; w += x * x + y * y
        else
            cu, su, m2 = _unit2(x, y)
            wt = p == 1 ? sqrt(m2) : 1.0
            c += wt * cu; s += wt * su; w += wt
        end
        n += 1
    end
    return _result(s, c, w, n, hypot(s, c) / max(w, TINY))
end

# Mean shift to the mode of a von Mises kernel density of doubled angles:
# weights exp(kappa (cos(2(t - mu)) - 1)), optionally times |g|^2.
function est_meanshift(gx, gy, kappa::Float64; iters::Int=3, p::Int=0, init::Symbol=:mean)
    s = 0.0; c = 0.0; n = 0
    @inbounds for k in eachindex(gx, gy)
        x = Float64(gx[k]); y = Float64(gy[k])
        _ok(x, y) || continue
        cu, su, m2 = _unit2(x, y)
        wt = p == 2 ? m2 : 1.0
        c += wt * cu; s += wt * su; n += 1
    end
    n == 0 && return WindowFitResult(failure_reason=:empty_sample)
    r0 = hypot(s, c)
    r0 <= TINY && return WindowFitResult(rbar2=0.0, effective_n=n, failure_reason=:degenerate_resultant)
    mc = c / r0; ms = s / r0
    if init === :hist
        # Start from the peak of a coarse 18-bin (10 deg) histogram smoothed with
        # [1, 2, 1], so the iteration starts on a mode rather than between two.
        h = zeros(18)
        @inbounds for k in eachindex(gx, gy)
            x = Float64(gx[k]); y = Float64(gy[k])
            _ok(x, y) || continue
            t = mod(atan(y, x), pi)
            h[min(18, floor(Int, t / pi * 18) + 1)] += 1
        end
        hs = [h[mod1(i - 1, 18)] + 2h[i] + h[mod1(i + 1, 18)] for i in 1:18]
        t0 = (argmax(hs) - 0.5) * pi / 18
        mc = cos(2t0); ms = sin(2t0)
    end
    conf = 0.0
    for _ in 1:iters
        s2 = 0.0; c2 = 0.0; w = 0.0; wu = 0.0
        @inbounds for k in eachindex(gx, gy)
            x = Float64(gx[k]); y = Float64(gy[k])
            _ok(x, y) || continue
            cu, su, m2 = _unit2(x, y)
            base = p == 2 ? m2 : 1.0
            wt = base * exp(kappa * (cu * mc + su * ms - 1))
            c2 += wt * cu; s2 += wt * su; w += wt
        end
        r = hypot(s2, c2)
        r <= TINY * max(w, TINY) && break
        mc = c2 / r; ms = s2 / r
        conf = r / w
    end
    a = _wrap(0.5 * atan(ms, mc))
    return WindowFitResult(angle=a, rbar2=conf, effective_n=n, low_confidence=conf < 0.5,
                           fit_ok=isfinite(a), failure_reason=:none)
end

# Inverse of A1(kappa) = I1/I0 (Best & Fisher 1981 approximation).
@inline function a1inv(r)
    r < 0.53 && return 2r + r^3 + 5r^5 / 6
    r < 0.85 && return -0.4 + 1.39r + 0.43 / (1 - r)
    return 1 / (r^3 - 4r^2 + 3r)
end

# Maximum likelihood von Mises (on doubled angles) plus uniform background via EM.
# The analogue of the prototype's Gaussian + constant histogram model, unbinned.
function est_vmu(gx, gy; iters::Int=8, kmax::Float64=50.0)
    s = 0.0; c = 0.0; n = 0
    @inbounds for k in eachindex(gx, gy)
        x = Float64(gx[k]); y = Float64(gy[k])
        _ok(x, y) || continue
        cu, su, _ = _unit2(x, y)
        c += cu; s += su; n += 1
    end
    n == 0 && return WindowFitResult(failure_reason=:empty_sample)
    r0 = hypot(s, c)
    r0 <= TINY && return WindowFitResult(rbar2=0.0, effective_n=n, failure_reason=:degenerate_resultant)
    mc = c / r0; ms = s / r0
    rb = r0 / n
    pk = 0.5
    kap = clamp(a1inv(min(rb / pk, 0.98)), 0.5, kmax)
    for _ in 1:iters
        # log normaliser of the von Mises on the doubled-angle circle: I0(kap);
        # use the scaled form exp(kap (cos - 1)) / I0e(kap).
        i0e = _besseli0e(kap)
        s2 = 0.0; c2 = 0.0; wsum = 0.0
        @inbounds for k in eachindex(gx, gy)
            x = Float64(gx[k]); y = Float64(gy[k])
            _ok(x, y) || continue
            cu, su, _ = _unit2(x, y)
            fv = pk * exp(kap * (cu * mc + su * ms - 1)) / i0e
            rsp = fv / (fv + (1 - pk))
            c2 += rsp * cu; s2 += rsp * su; wsum += rsp
        end
        r = hypot(s2, c2)
        (wsum <= 1e-9 || r <= TINY) && break
        mc = c2 / r; ms = s2 / r
        pk = clamp(wsum / n, 1e-3, 1.0)
        kap = clamp(a1inv(min(r / wsum, 0.995)), 0.05, kmax)
    end
    a = _wrap(0.5 * atan(ms, mc))
    return WindowFitResult(angle=a, rbar2=pk, effective_n=n, low_confidence=pk < 0.5,
                           fit_ok=isfinite(a), failure_reason=:none)
end

# exp(-x) I0(x), polynomial approximation (Abramowitz & Stegun 9.8.1-2).
function _besseli0e(x::Float64)
    ax = abs(x)
    if ax < 3.75
        t = (x / 3.75)^2
        return exp(-ax) * (1 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492 + t * (0.2659732 +
                           t * (0.0360768 + t * 0.0045813))))))
    end
    t = 3.75 / ax
    return (0.39894228 + t * (0.01328592 + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281 +
            t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377)))))))) / sqrt(ax)
end

# Histogram + Gaussian(+constant) peak fit written for this audit (not the
# prototype code): bins on the axial circle, histogram recentred on its peak bin
# so the fit never straddles the wrap, LsqFit from the peak. Used only to ablate
# which ingredient of the Gaussian baseline matters (offset term, bin count).
using LsqFit
function est_hist_gauss(gx, gy; nbins::Int=100, offset::Bool=true)
    counts = zeros(nbins)
    n = 0
    @inbounds for k in eachindex(gx, gy)
        x = Float64(gx[k]); y = Float64(gy[k])
        _ok(x, y) || continue
        t = atan(y, x)
        t = mod(t + pi / 2, pi)                      # [0, pi)
        b = min(nbins, floor(Int, t / pi * nbins) + 1)
        counts[b] += 1; n += 1
    end
    n == 0 && return WindowFitResult(failure_reason=:empty_sample)
    width = pi / nbins
    pkb = argmax(counts)
    shift = div(nbins, 2) - pkb                      # move the peak bin to the middle
    y = circshift(counts, shift) ./ n
    xc = [(-pi / 2) + (i - 0.5) * width for i in 1:nbins]   # centred grid, peak near 0
    m(x, p) = offset ? p[1] .* exp.(-(x .- p[2]) .^ 2 .* p[3]) .+ p[4] :
                       p[1] .* exp.(-(x .- p[2]) .^ 2 .* p[3])
    p0 = offset ? [maximum(y), 0.0, 10.0, 0.0] : [maximum(y), 0.0, 10.0]
    local mu
    try
        f = curve_fit(m, xc, y, p0)
        mu = f.param[2]
    catch
        return WindowFitResult(failure_reason=:fit_failed, effective_n=n)
    end
    # Undo the shift: bin centre of pkb in original coordinates.
    centre0 = (-pi / 2) + (div(nbins, 2) - 0.5) * width   # where the peak bin was moved
    t0 = (pkb - 0.5) * width - pi / 2                       # original peak-bin centre on [-pi/2, pi/2)
    a = _wrap(t0 + (mu - centre0))
    return WindowFitResult(angle=a, rbar2=NaN, effective_n=n, fit_ok=isfinite(a),
                           failure_reason=isfinite(a) ? :none : :fit_failed)
end
