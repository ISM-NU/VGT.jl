# Q1 diagnostic: why does a histogram peak fit beat the doubled-angle mean on
# filament scenes? Separates noise background from in-window geometry by
# rerunning with the planted noise switched off, and measures the signed pull
# of each estimator toward the other filament family (crossing scene).
#
#   julia --project=scripts/estimator_audit -t 16 scripts/estimator_audit/diag_q1.jl

using Printf, Statistics, StatsBase, LsqFit
using VGT
using VGT: VGT_QU_mv_2d_profiled, headless_angle_difference, WindowFitResult

const HERE = @__DIR__
const REPO = normpath(joinpath(HERE, "..", ".."))
module MW end
for e in Meta.parseall(read(joinpath(REPO, "scripts", "validate_mw_like_synthetic.jl"), String)).args
    e isa LineNumberNode && continue
    s = string(e)
    (occursin("quickactivate", s) || startswith(s, "using DrWatson") || startswith(s, "if abspath")) && continue
    Core.eval(MW, e)
end
include(joinpath(HERE, "estimators.jl"))

function estimate(mode::Symbol, gx, gy)
    mode === :vm && return est_weighted(gx, gy, 0)
    mode === :ms4 && return est_meanshift(gx, gy, 4.0)
    mode === :vmu && return est_vmu(gx, gy)
    mode === :hg100 && return est_hist_gauss(gx, gy; nbins=100, offset=true)
    error("unknown $mode")
end
@eval VGT fit_window_result(mode::Symbol, gx::AbstractArray, gy::AbstractArray) = Main.estimate(mode, gx, gy)

# Count exactly-zero gradients that production would keep as angle 0.
function zero_gradient_fraction(cube, noise)
    nz = 0; nvalid = 0
    nx, ny, nv = size(cube)
    for k in 1:nv, j in 2:ny-1, i in 2:nx-1
        all(cube[ii, jj, k] >= noise for (ii, jj) in ((i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1))) || continue
        nvalid += 1
        gx = (cube[i-1, j, k] - cube[i+1, j, k]) / 2
        gy = (cube[i, j-1, k] - cube[i, j+1, k]) / 2
        (gx == 0 && gy == 0) && (nz += 1)
    end
    return nz, nvalid
end

function run(label, builder; noise_scale=1.0, size_=96)
    MW.NOISE_SCALE[] = noise_scale
    cube, truth, noise = builder(; nx=size_, ny=size_)
    ta, tr = MW.truth_angle_and_rbar(truth.q, truth.u, truth.w)
    cut = quantile([v for v in vec(truth.w) if v > 0], 0.55)
    clean = (truth.w .> cut) .& isfinite.(ta) .& (tr .>= 0.75)
    nz, nvalid = zero_gradient_fraction(cube, noise)
    res = Dict{Symbol, Vector{Float64}}()
    for est in (:vm, :ms4, :vmu, :hg100)
        q, u, _ = Base.invokelatest(VGT_QU_mv_2d_profiled, cube, 12, 1, size(cube, 3), noise;
                                    min_valid=50, fitter_mode=est, gradient_mode=:central_checked)
        b = [(q[i] == 0 && u[i] == 0) ? NaN : 0.5 * atan(u[i], q[i]) + pi / 2 for i in eachindex(q)]
        # signed error in degrees, est - truth
        res[est] = [clean[i] && isfinite(b[i]) ? rad2deg(headless_angle_difference(b[i], ta[i])) : NaN for i in eachindex(b)]
    end
    @printf("%-20s noise x%.0f  zero-gradient %d/%d  ", label, noise_scale, nz, nvalid)
    for est in (:vm, :ms4, :vmu, :hg100)
        e = filter(isfinite, res[est])
        @printf("%s med|e| %.2f  ", est, median(abs.(e)))
    end
    println()
    if label == "crossing_filaments"
        # Pull toward the other family: for truth near 25 deg the other family
        # (105 deg) lies at +80 deg; for truth near 105 deg at -80 deg.
        sgn = [isfinite(ta[i]) ? (abs(rad2deg(headless_angle_difference(ta[i], deg2rad(25)))) < 20 ? 1.0 :
               abs(rad2deg(headless_angle_difference(ta[i], deg2rad(105)))) < 20 ? -1.0 : NaN) : NaN for i in eachindex(ta)]
        for est in (:vm, :ms4, :vmu, :hg100)
            p = filter(isfinite, res[est] .* sgn)
            @printf("    %-6s mean signed pull toward other family %+.2f deg (median %+.2f, n=%d)\n",
                    est, mean(p), median(p), length(p))
        end
    end
end

if get(ENV, "DIAG_CROSS", "0") == "1"
    # |error| by distance to the nearest crossing point of the two families.
    for ns in (1.0, 0.0)
        MW.NOISE_SCALE[] = ns
        cube, truth, noise = MW.build_crossing_filaments(; nx=96, ny=96)
        ta, tr = MW.truth_angle_and_rbar(truth.q, truth.u, truth.w)
        cut = quantile([v for v in vec(truth.w) if v > 0], 0.55)
        clean = (truth.w .> cut) .& isfinite.(ta) .& (tr .>= 0.75)
        pts = Tuple{Float64,Float64}[]
        # line k: p = c_k + t d_k with c_k = (48,48) + o n_k
        L = [(48 + o * cos(a + pi / 2), 48 + o * sin(a + pi / 2), cos(a), sin(a)) for a in (deg2rad(25), deg2rad(105)) for o in (-26, 0, 26)]
        for l1 in L[1:3], l2 in L[4:6]
            A = [l1[3] -l2[3]; l1[4] -l2[4]]; t = A \ [l2[1] - l1[1], l2[2] - l1[2]]
            push!(pts, (l1[1] + t[1] * l1[3], l1[2] + t[1] * l1[4]))
        end
        dist = [minimum(hypot(i - p[1], j - p[2]) for p in pts) for i in 1:96, j in 1:96]
        for est in (:vm, :ms4, :hg100)
            q, u, _ = Base.invokelatest(VGT_QU_mv_2d_profiled, cube, 12, 1, size(cube, 3), noise;
                                        min_valid=50, fitter_mode=est, gradient_mode=:central_checked)
            b = [(q[i] == 0 && u[i] == 0) ? NaN : 0.5 * atan(u[i], q[i]) + pi / 2 for i in eachindex(q)]
            e = [clean[i] && isfinite(b[i]) ? abs(rad2deg(headless_angle_difference(b[i], ta[i]))) : NaN for i in eachindex(b)]
            @printf("cross noise x%.0f %-6s", ns, est)
            for (lo, hi) in ((0, 6), (6, 9), (9, 12), (12, 18), (18, 1e9))
                m = [isfinite(e[i]) && lo <= dist[i] < hi for i in eachindex(e)]
                @printf("  d[%g,%g) med %.2f (n=%d)", lo, hi, median(e[m]), count(m))
            end
            println()
        end
    end
    exit()
end

for (label, b) in (("many_filaments", MW.build_many_filaments), ("crossing_filaments", MW.build_crossing_filaments),
                   ("curved_fibers", MW.build_curved_fibers), ("shells_arcs", MW.build_shells_arcs),
                   ("sparse_co_like", MW.build_sparse_co_like))
    for ns in (1.0, 0.0)
        run(label, b; noise_scale=ns)
    end
end
