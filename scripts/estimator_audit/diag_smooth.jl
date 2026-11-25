# Where does the old Gaussian's advantage in the smoothed product (ker=4) sit on
# the crossing scene: clean or ambiguous pixels, near or far from crossings?
# Also reports how often psi_reliability exceeds 1.
#   AUDIT_OLD_FITTERS=<private fitters.jl> julia --project=scripts/estimator_audit -t 16 scripts/estimator_audit/diag_smooth.jl
ENV["AUDIT_SEEDS"] = "0:-1"; ENV["AUDIT_OUT"] = tempname()
include(joinpath(@__DIR__, "run_audit.jl"))
load_old_gauss(); patch_kernel!()
for (label, b) in (("crossing_filaments", MW.build_crossing_filaments), ("curved_fibers", MW.build_curved_fibers),
                   ("sparse_co_like", MW.build_sparse_co_like), ("velocity_components", MW.build_velocity_components))
    cube, truth, noise = b(; nx=96, ny=96)
    ta, tr = MW.truth_angle_and_rbar(truth.q, truth.u, truth.w)
    cut = quantile([v for v in vec(truth.w) if v > 0], 0.55)
    base = (truth.w .> cut) .& isfinite.(ta)
    pts = Tuple{Float64,Float64}[]
    if label == "crossing_filaments"
        L = [(48 + o * cos(a + pi / 2), 48 + o * sin(a + pi / 2), cos(a), sin(a)) for a in (deg2rad(25), deg2rad(105)) for o in (-26, 0, 26)]
        for l1 in L[1:3], l2 in L[4:6]
            t = [l1[3] -l2[3]; l1[4] -l2[4]] \ [l2[1] - l1[1], l2[2] - l1[2]]
            push!(pts, (l1[1] + t[1] * l1[3], l1[2] + t[1] * l1[4]))
        end
    end
    dist = isempty(pts) ? fill(Inf, 96, 96) : [minimum(hypot(i - p[1], j - p[2]) for p in pts) for i in 1:96, j in 1:96]
    I = dropdims(sum(cube; dims=3); dims=3)
    for est in (:vm, :gauss, :ms4h, :hg100, :hg100n)
        q, u, _ = Base.invokelatest(VGT_QU_mv_2d_profiled, cube, 12, 1, size(cube, 3), noise; min_valid=50, fitter_mode=est, gradient_mode=:central_checked)
        fin = finalize_vgt_products(copy(q), copy(u), copy(I), 4)
        e = reshape([base[i] && isfinite(fin.psi[i]) ? abs(rad2deg(headless_angle_difference(fin.psi[i], ta[i]))) : NaN for i in eachindex(ta)], size(ta))
        rel = filter(isfinite, vec(fin.reliability))
        @printf("%-20s %-7s smoothed clean %.2f ambiguous %.2f", label, est,
                median(e[isfinite.(e) .& (tr .>= 0.75)]), median(e[isfinite.(e) .& (tr .< 0.5)]))
        if !isempty(pts)
            for (lo, hi) in ((0, 9), (9, 18), (18, Inf))
                m = isfinite.(e) .& (tr .>= 0.75) .& (lo .<= dist .< hi)
                @printf("  clean d[%g,%g) %.2f", lo, hi, median(e[m]))
            end
        end
        @printf("  reliability max %.3f frac>1 %.4f\n", maximum(rel), count(>(1), rel) / length(rel))
    end
end
