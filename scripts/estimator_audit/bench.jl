# Per-window cost of each estimator on a 12x12 window (single thread).
#   AUDIT_OLD_FITTERS=<private fitters.jl> julia --project=scripts/estimator_audit scripts/estimator_audit/bench.jl
using Printf, Random, Statistics, LsqFit, StatsBase
using VGT
include(joinpath(@__DIR__, "estimators.jl"))
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
haskey(ENV, "AUDIT_OLD_FITTERS") && Base.include(OldGauss, ENV["AUDIT_OLD_FITTERS"])

rng = MersenneTwister(1)
nw = 400
wins = map(1:nw) do _
    t = rand(rng) * pi
    gx = cos(t) .+ 0.6 .* randn(rng, 12, 12); gy = sin(t) .+ 0.6 .* randn(rng, 12, 12)
    gx[rand(rng, 1:144, 30)] .= NaN
    (gx, gy)
end
fns = [
    ("vm (production)", (x, y) -> VGT.fit_window_result(:axial_vonmises, x, y)),
    ("st2", (x, y) -> est_weighted(x, y, 2)),
    ("ms2", (x, y) -> est_meanshift(x, y, 2.0)),
    ("ms4", (x, y) -> est_meanshift(x, y, 4.0)),
    ("vmu", (x, y) -> est_vmu(x, y)),
    ("hg100 (own)", (x, y) -> est_hist_gauss(x, y)),
]
isdefined(OldGauss, :fit_window_result) &&
    push!(fns, ("gauss (old repo)", (x, y) -> Base.invokelatest(OldGauss.fit_window_result, :lsqfit_reference, x, y, 100)))
for (name, f) in fns
    f(wins[1]...)
    reps = occursin("g", split(name)[1]) && occursin("h", name) || occursin("gauss", name) ? 1 : 50
    t = @elapsed for _ in 1:reps, w in wins
        f(w...)
    end
    @printf("%-18s %9.2f us/window\n", name, 1e6 * t / (reps * nw))
end
