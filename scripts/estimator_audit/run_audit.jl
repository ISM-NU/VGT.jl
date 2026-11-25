# Estimator audit harness (scratch; not part of the package).
#
# Runs the production kernel VGT_QU_mv_2d_profiled with candidate window
# estimators swapped in through fit_window_result, on the MW-like morphology
# suite and the uniform/sparse/varying suite, over several seeds, and writes
# per-(scenario, seed, estimator, region) error rows plus confidence-map
# calibration rows.
#
#   julia --project=scripts/estimator_audit -t 16 scripts/estimator_audit/run_audit.jl
#
# Environment:
#   AUDIT_SUITE      mw | uni                         (default mw)
#   AUDIT_SCEN       comma list of scenario labels     (default all of the suite)
#   AUDIT_SEEDS      e.g. 0:9 ; seed 0 = the scripts' documented default seed
#   AUDIT_EST        comma list (vm,gauss,st1,st2,ms2,ms4,st2ms2,vmu,hg100,hg100n,hg36)
#   AUDIT_SIZE       map size for the MW suite (default 96)
#   AUDIT_OUT        output prefix (default validation_output/estimator_audit/run)
#   AUDIT_OLD_FITTERS path to the archived private repository's src/science/fitters.jl;
#                    needed only for the `gauss` baseline. That file is loaded at
#                    run time and is deliberately not copied into this repository
#                    (docs/provenance.md).
#   AUDIT_PARITY=1   check that `vm` through the patched path is bitwise equal to
#                    the unpatched package on the first cube.

using Printf, Random, Statistics, StatsBase, LsqFit
using VGT
using VGT: VGT_QU_mv_2d_profiled, headless_angle_difference, finalize_vgt_products, WindowFitResult

const HERE = @__DIR__
const REPO = normpath(joinpath(HERE, "..", ".."))

# Load builder code from the validation scripts without their DrWatson activation.
function load_script_defs(mod::Module, path::AbstractString)
    src = read(path, String)
    ex = Meta.parseall(src)
    for e in ex.args
        e isa LineNumberNode && continue
        s = string(e)
        (occursin("quickactivate", s) || startswith(s, "using DrWatson") ||
         startswith(s, "if abspath(PROGRAM_FILE)") || s == "main()") && continue
        Core.eval(mod, e)
    end
end

module MW end
module UNI end
load_script_defs(MW, joinpath(REPO, "scripts", "validate_mw_like_synthetic.jl"))
load_script_defs(UNI, joinpath(REPO, "scripts", "validate_synthetic.jl"))

include(joinpath(HERE, "estimators.jl"))

# Gaussian baseline from the archived private repository, loaded at run time.
module OldGauss
using LsqFit, StatsBase, Statistics
const HIST100_EDGES = collect(range(-pi / 2, stop=pi / 2, length=101))
const HIST100_AX = 0.5 .* (HIST100_EDGES[1:end-1] .+ HIST100_EDGES[2:end])
const HIST100_AX_SHIFTED = HIST100_AX .- pi / 2
const LEGACY_HIST100_EDGES = [(-pi / 2) + (i - 1) * (pi / 101) for i in 1:101]
const LEGACY_HIST100_AX = 0.5 .* (LEGACY_HIST100_EDGES[1:end-1] .+ LEGACY_HIST100_EDGES[2:end])
const LEGACY_HIST100_AX_SHIFTED = LEGACY_HIST100_AX .- pi / 2
const FIT_LOW_CONFIDENCE_RBAR2_THRESHOLD = 0.5
is_interrupt_exception(e) = e isa InterruptException
end
const HAVE_GAUSS = Ref(false)
function load_old_gauss()
    p = get(ENV, "AUDIT_OLD_FITTERS", "")
    isfile(p) || error("gauss requested but AUDIT_OLD_FITTERS=$(repr(p)) is not a file")
    Base.include(OldGauss, p)
    HAVE_GAUSS[] = true
end
function gauss_old(gx, gy)
    r = Base.invokelatest(OldGauss.fit_window_result, :lsqfit_reference, gx, gy, 100)
    return WindowFitResult(angle=r.angle, rbar2=r.rbar2, effective_n=r.effective_n,
                           low_confidence=r.low_confidence, fit_ok=r.fit_ok, failure_reason=r.failure_reason)
end

const PRODUCTION_FIT = VGT.fit_window_result  # same generic function; capture stats helper below
vm_production(gx, gy) = begin
    stats = VGT._axial_resultant_stats_paired(gx, gy)
    stats.effective_n == 0 && return VGT._failed_fit_result(:empty_sample)
    reason = stats.resultant <= sqrt(eps(Float64)) ? :degenerate_resultant : :fit_failed
    VGT._fit_result(stats.angle, stats; failure_reason=reason)
end

function estimate(mode::Symbol, gx, gy)
    mode === :vm && return vm_production(gx, gy)
    mode === :gauss && return gauss_old(gx, gy)
    mode === :st1 && return est_weighted(gx, gy, 1)
    mode === :st2 && return est_weighted(gx, gy, 2)
    mode === :ms2 && return est_meanshift(gx, gy, 2.0)
    mode === :ms4 && return est_meanshift(gx, gy, 4.0)
    mode === :st2ms2 && return est_meanshift(gx, gy, 2.0; p=2)
    mode === :vmu && return est_vmu(gx, gy)
    mode === :hg100 && return est_hist_gauss(gx, gy; nbins=100, offset=true)
    mode === :hg100n && return est_hist_gauss(gx, gy; nbins=100, offset=false)
    mode === :hg36 && return est_hist_gauss(gx, gy; nbins=36, offset=true)
    throw(ArgumentError("unknown estimator $mode"))
end

function patch_kernel!()
    @eval VGT function fit_window_result(mode::Symbol, grad_x::AbstractArray, grad_y::AbstractArray)
        return Main.estimate(mode, grad_x, grad_y)
    end
end

# ---- scenarios -------------------------------------------------------------

const MW_SCEN = [
    ("many_filaments", :build_many_filaments, 101),
    ("crossing_filaments", :build_crossing_filaments, 102),
    ("curved_fibers", :build_curved_fibers, 103),
    ("velocity_components", :build_velocity_components, 104),
    ("shells_arcs", :build_shells_arcs, 105),
    ("planck_misaligned", :build_planck_misaligned, 106),
    ("sparse_co_like", :build_sparse_co_like, 107),
]

# seed 0 -> the script default, otherwise default*1000 + s
mwseed(default, s) = s == 0 ? default : default * 1000 + s

const UNI_SCEN = [
    ("uniform_n0.05", rng -> UNI.make_uniform_cube(; noise_level=0.05, rng=rng), :const),
    ("uniform_n0.12", rng -> UNI.make_uniform_cube(; noise_level=0.12, rng=rng), :const),
    ("uniform_n0.20", rng -> UNI.make_uniform_cube(; noise_level=0.20, rng=rng), :const),
    ("weak_amp0.05", rng -> UNI.make_uniform_cube(; gradient_amplitude=0.05, noise_level=0.05, rng=rng), :const),
    ("sparse_fill0.35", rng -> UNI.make_sparse_cube(; fill_fraction=0.35, noise_level=0.05, rng=rng), :const),
    ("varying_c8e-4", rng -> UNI.make_varying_angle_cube(; noise_level=0.05, curvature=8e-4, rng=rng), :vary8),
    ("varying_c2e-3", rng -> UNI.make_varying_angle_cube(; noise_level=0.05, curvature=2.0e-3, rng=rng), :vary2),
]

# ---- scoring ---------------------------------------------------------------

function window_truth_rbar(tq, tu, tw, dn)
    # Truth coherence over the same dn x dn window the kernel uses at each centre.
    nx, ny = size(tq); h = dn ÷ 2
    out = fill(NaN, nx, ny)
    for jj in 1:ny, ii in 1:nx
        (ii < h || ii > nx - h || jj < h || jj > ny - h) && continue
        q = 0.0; u = 0.0; w = 0.0
        for j in jj-h+1:jj+h, i in ii-h+1:ii+h
            q += tq[i, j]; u += tu[i, j]; w += tw[i, j]
        end
        w > 0 && (out[ii, jj] = hypot(q, u) / w)
    end
    return out
end

function errs_in(est, truth, mask)
    v = Float64[]; elig = 0
    for idx in eachindex(est)
        mask[idx] || continue
        elig += 1
        isfinite(est[idx]) || continue
        push!(v, abs(rad2deg(headless_angle_difference(est[idx], truth[idx]))))
    end
    return v, elig
end

q_(v, p) = isempty(v) ? NaN : quantile(v, p)

function decile_medians(err, conf)
    err = vec(err); conf = vec(conf)
    ok = isfinite.(err) .& isfinite.(conf)
    e = err[ok]; c = conf[ok]
    length(e) < 50 && return fill(NaN, 10), NaN
    o = sortperm(c)
    edges = round.(Int, range(0, length(o), length=11))
    meds = [median(e[o[edges[k]+1:edges[k+1]]]) for k in 1:10]
    return meds, corspearman(e, c)
end

function main()
    suite = Symbol(get(ENV, "AUDIT_SUITE", "mw"))
    seeds = eval(Meta.parse(get(ENV, "AUDIT_SEEDS", "0:0")))
    ests = Symbol.(strip.(split(get(ENV, "AUDIT_EST", "vm,st2,ms2,vmu"), ",")))
    size_ = parse(Int, get(ENV, "AUDIT_SIZE", "96"))
    out = get(ENV, "AUDIT_OUT", joinpath(REPO, "validation_output", "estimator_audit", "run"))
    mkpath(dirname(out))
    :gauss in ests && load_old_gauss()
    scen_all = suite === :mw ? [s[1] for s in MW_SCEN] : [s[1] for s in UNI_SCEN]
    scen_sel = haskey(ENV, "AUDIT_SCEN") ? strip.(split(ENV["AUDIT_SCEN"], ",")) : scen_all
    dn = suite === :mw ? 12 : 10
    min_valid = 50
    ker = 4

    parity_done = get(ENV, "AUDIT_PARITY", "0") != "1"
    ref_q = nothing
    patched = false

    ferr = open(out * "_errors.csv", "w")
    fcal = open(out * "_calib.csv", "w")
    println(ferr, "suite,scenario,seed,estimator,reference,region,eligible,n,support,med,p90,mean,secs")
    println(fcal, "suite,scenario,seed,estimator,map,spearman," * join(["d$k" for k in 1:10], ","))

    for s in seeds, label in scen_sel
        if suite === :mw
            _, bsym, dseed = MW_SCEN[findfirst(x -> x[1] == label, MW_SCEN)]
            cube, truth, noise = getfield(MW, bsym)(; nx=size_, ny=size_, seed=mwseed(dseed, s))
            refs = [("total", MW.truth_angle_and_rbar(truth.q, truth.u, truth.w)..., truth.q, truth.u, truth.w)]
            label == "planck_misaligned" &&
                push!(refs, ("foreground", MW.truth_angle_and_rbar(truth.fg_q, truth.fg_u, truth.fg_w)...,
                             truth.fg_q, truth.fg_u, truth.fg_w))
        else
            _, fn, kind = UNI_SCEN[findfirst(x -> x[1] == label, UNI_SCEN)]
            cube = fn(MersenneTwister(1000 + s))
            noise = 1.0
            nx, ny, _ = size(cube)
            ta = kind === :const ? fill(UNI.TRUE_ANGLE, nx, ny) :
                 [UNI.varying_true_angle(i, j, nx, ny; curvature=(kind === :vary8 ? 8e-4 : 2.0e-3)) for i in 1:nx, j in 1:ny]
            refs = [("total", ta, ones(nx, ny), nothing, nothing, nothing)]
        end
        nv = size(cube, 3)
        intensity = dropdims(sum(cube; dims=3); dims=3)
        if !parity_done
            q0, u0, _ = VGT_QU_mv_2d_profiled(cube, dn, 1, nv, noise; min_valid=min_valid,
                                              fitter_mode=:axial_vonmises, gradient_mode=:central_checked)
            ref_q = (q0, u0)
        end
        if !patched
            patch_kernel!(); patched = true
        end
        for est in ests
            t = @elapsed q, u, prof = Base.invokelatest(VGT_QU_mv_2d_profiled, cube, dn, 1, nv, noise;
                                                        min_valid=min_valid, fitter_mode=est,
                                                        gradient_mode=:central_checked)
            if !parity_done && est === :vm
                ok = q == ref_q[1] && u == ref_q[2]
                println("PARITY vm patched vs package bitwise equal: ", ok)
                ok || error("parity failed")
                parity_done = true
            end
            gang = [(q[i] == 0 && u[i] == 0) ? NaN : 0.5 * atan(u[i], q[i]) for i in eachindex(q)]
            # MW truth is the filament tangent: compare the B proxy (gradient + 90 deg).
            raw = suite === :mw ? gang .+ pi / 2 : gang
            fin = finalize_vgt_products(copy(q), copy(u), copy(intensity), ker)
            smooth = suite === :mw ? fin.psi : fin.psi .- pi / 2
            fr = prof[:fit_rbar2]
            for (refname, tang, trbar, tq, tu, tw) in refs
                if suite === :mw
                    sig = tw
                    cut = q_([v for v in vec(sig) if v > 0], 0.55)
                    base = (sig .> cut) .& isfinite.(tang)
                    wr = window_truth_rbar(tq, tu, tw, dn)
                    regions = [("clean", base .& (trbar .>= 0.75)), ("all_signal", base),
                               ("ambiguous", base .& (trbar .< 0.5)),
                               ("window_clean", base .& (wr .>= 0.75))]
                else
                    h = dn ÷ 2; nx, ny = size(tang)
                    regions = [("all", [h + 1 <= i <= nx - h - 1 && h + 1 <= j <= ny - h - 1 for i in 1:nx, j in 1:ny])]
                end
                for (rname, m) in regions
                    v, elig = errs_in(raw, tang, m)
                    @printf(ferr, "%s,%s,%d,%s,%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.2f\n", suite, label, s, est,
                            refname, rname, elig, length(v), elig == 0 ? NaN : length(v) / elig,
                            q_(v, 0.5), q_(v, 0.9), isempty(v) ? NaN : mean(v), t)
                end
                # Calibration on the broad region (all_signal or all).
                m = regions[suite === :mw ? 2 : 1][2]
                eraw = [m[i] && isfinite(raw[i]) ? abs(rad2deg(headless_angle_difference(raw[i], tang[i]))) : NaN for i in eachindex(raw)]
                esm = [m[i] && isfinite(smooth[i]) ? abs(rad2deg(headless_angle_difference(smooth[i], tang[i]))) : NaN for i in eachindex(raw)]
                for (mapname, e, c) in (("fit_conf_vs_raw_err", eraw, fr),
                                        ("reliability_vs_smoothed_err", esm, fin.reliability),
                                        ("support_vs_raw_err", eraw, prof[:fit_support]))
                    meds, rho = decile_medians(e, c)
                    @printf(fcal, "%s,%s,%d,%s,%s,%.4f,%s\n", suite, label, s, est, mapname, rho,
                            join([@sprintf("%.3f", x) for x in meds], ","))
                end
                # Also the smoothed-map errors in the main region, for the report.
                v, elig = errs_in(smooth, tang, m)
                @printf(ferr, "%s,%s,%d,%s,%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.2f\n", suite, label, s, est,
                        refname, "smoothed_" * (suite === :mw ? "all_signal" : "all"), elig, length(v),
                        elig == 0 ? NaN : length(v) / elig, q_(v, 0.5), q_(v, 0.9), isempty(v) ? NaN : mean(v), t)
            end
            flush(ferr); flush(fcal)
            @printf("%-20s seed %2d %-7s %6.1fs\n", label, s, est, t)
        end
    end
    close(ferr); close(fcal)
end

main()
