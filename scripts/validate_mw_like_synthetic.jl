using DrWatson
quickactivate(@__DIR__, "VGT")

using Printf
using Random
using Statistics
using VGT: VGT_QU_mv_2d_profiled, headless_angle_difference

const DEG = pi / 180

function selected_symbols(env_name::AbstractString, default)
    raw = strip(get(ENV, env_name, ""))
    isempty(raw) && return default
    return Tuple(Symbol(strip(part)) for part in split(raw, ",") if !isempty(strip(part)))
end

const DEFAULT_MIN_VALID_CANDIDATES = (
    ("frac50", 72),
    ("frac35", 50),
    ("frac25", 36),
)
const MIN_VALID_LABELS = Dict(value => label for (label, value) in DEFAULT_MIN_VALID_CANDIDATES)

function selected_min_valid_candidates()
    raw = strip(get(ENV, "VGT_MW_MIN_VALID", ""))
    isempty(raw) && return DEFAULT_MIN_VALID_CANDIDATES
    values = [parse(Int, strip(part)) for part in split(raw, ",") if !isempty(strip(part))]
    return Tuple((get(MIN_VALID_LABELS, value, "min$(value)"), value) for value in values)
end

const MODES = selected_symbols("VGT_MW_MODES", (:nan_edge_backward, :central_checked, :sobel_checked))
const FITTERS = selected_symbols("VGT_MW_FITTERS", (:axial_vonmises,))
const MIN_VALID_CANDIDATES = selected_min_valid_candidates()
const MAP_SIZE = parse(Int, get(ENV, "VGT_MW_SIZE", "128"))

Base.@kwdef struct TruthMaps
    q::Matrix{Float64}
    u::Matrix{Float64}
    w::Matrix{Float64}
    fg_q::Matrix{Float64}
    fg_u::Matrix{Float64}
    fg_w::Matrix{Float64}
end

function TruthMaps(nx::Int, ny::Int)
    return TruthMaps(
        q=zeros(nx, ny),
        u=zeros(nx, ny),
        w=zeros(nx, ny),
        fg_q=zeros(nx, ny),
        fg_u=zeros(nx, ny),
        fg_w=zeros(nx, ny),
    )
end

function add_truth!(truth::TruthMaps, signal::AbstractMatrix, angle::AbstractMatrix;
                    weight::Real=1.0, foreground::Bool=false)
    @inbounds for idx in eachindex(signal, angle)
        s = Float64(weight) * signal[idx]
        s > 0 || continue
        θ = angle[idx]
        truth.q[idx] = muladd(s, cos(2θ), truth.q[idx])
        truth.u[idx] = muladd(s, sin(2θ), truth.u[idx])
        truth.w[idx] += s
        if foreground
            truth.fg_q[idx] = muladd(s, cos(2θ), truth.fg_q[idx])
            truth.fg_u[idx] = muladd(s, sin(2θ), truth.fg_u[idx])
            truth.fg_w[idx] += s
        end
    end
    return truth
end

function truth_angle_and_rbar(q::AbstractMatrix, u::AbstractMatrix, w::AbstractMatrix)
    angle = fill(NaN, size(q))
    rbar = fill(NaN, size(q))
    @inbounds for idx in eachindex(q, u, w)
        if w[idx] > 0
            angle[idx] = 0.5 * atan(u[idx], q[idx])
            rbar[idx] = sqrt(q[idx]^2 + u[idx]^2) / w[idx]
        end
    end
    return angle, rbar
end

function add_component!(cube, truth::TruthMaps, signal, angle, channels;
                        foreground::Bool=false)
    add_truth!(truth, signal, angle; weight=sum(channels), foreground=foreground)
    @inbounds for k in axes(cube, 3)
        w = channels[k]
        w == 0.0 && continue
        for j in axes(cube, 2), i in axes(cube, 1)
            cube[i, j, k] += w * signal[i, j]
        end
    end
    return nothing
end

function add_noise!(cube, sigma::Float64, rng)
    @inbounds for idx in eachindex(cube)
        cube[idx] += sigma * randn(rng)
    end
    return cube
end

function gaussian_segment(nx, ny; cx, cy, angle, length, width, amp, gap_period=0.0)
    signal = zeros(nx, ny)
    bangle = fill(angle, nx, ny)
    ca = cos(angle)
    sa = sin(angle)
    @inbounds for j in 1:ny, i in 1:nx
        x = i - cx
        y = j - cy
        s = x * ca + y * sa
        d = -x * sa + y * ca
        longitudinal = exp(-(s / length)^8)
        transverse = exp(-(d^2) / (2 * width^2))
        gap = gap_period > 0 ? (0.45 + 0.55 * max(0.0, sin(2pi * s / gap_period))) : 1.0
        signal[i, j] = amp * longitudinal * transverse * gap
    end
    return signal, bangle
end

function ring_arc(nx, ny; cx, cy, radius, width, amp, phi_min, phi_max)
    signal = zeros(nx, ny)
    angle = fill(NaN, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        x = i - cx
        y = j - cy
        r = hypot(x, y)
        phi = atan(y, x)
        in_arc = phi_min <= phi <= phi_max
        in_arc || continue
        s = amp * exp(-((r - radius)^2) / (2 * width^2))
        signal[i, j] = s
        angle[i, j] = phi + pi / 2
    end
    return signal, angle
end

function finite_angle(signal, angle, fallback)
    out = fill(fallback, size(signal))
    @inbounds for idx in eachindex(signal, angle)
        isfinite(angle[idx]) && (out[idx] = angle[idx])
    end
    return out
end

function chan_all(nv)
    return ones(nv)
end

function chan_range(nv, lo, hi)
    w = zeros(nv)
    w[lo:hi] .= 1.0
    return w
end

function build_many_filaments(; nx=128, ny=128, nv=24, seed=101)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    base = 33 * DEG
    for n in 1:9
        angle = base + randn(rng) * 10 * DEG
        cx = rand(rng, 18:nx-18)
        cy = rand(rng, 18:ny-18)
        signal, bangle = gaussian_segment(
            nx, ny;
            cx=cx, cy=cy, angle=angle,
            length=42 + 12 * rand(rng),
            width=1.8 + 1.6 * rand(rng),
            amp=2.0 + 1.5 * rand(rng),
        )
        add_component!(cube, truth, signal, bangle, chan_all(nv); foreground=true)
    end
    return add_noise!(cube, 0.20, rng), truth, 0.20
end

function build_crossing_filaments(; nx=128, ny=128, nv=24, seed=102)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    for angle in (25 * DEG, 105 * DEG)
        for offset in (-26, 0, 26)
            cx = nx / 2 + offset * cos(angle + pi / 2)
            cy = ny / 2 + offset * sin(angle + pi / 2)
            signal, bangle = gaussian_segment(
                nx, ny;
                cx=cx, cy=cy, angle=angle,
                length=58, width=2.2, amp=2.8,
            )
            add_component!(cube, truth, signal, bangle, chan_all(nv); foreground=true)
        end
    end
    return add_noise!(cube, 0.18, rng), truth, 0.18
end

function build_curved_fibers(; nx=128, ny=128, nv=24, seed=103)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    arcs = [
        (cx=34.0, cy=32.0, radius=58.0, phi_min=-0.15, phi_max=1.55, amp=2.7),
        (cx=96.0, cy=38.0, radius=52.0, phi_min=1.55, phi_max=3.05, amp=2.4),
        (cx=64.0, cy=110.0, radius=62.0, phi_min=-2.65, phi_max=-0.60, amp=2.1),
    ]
    for arc in arcs
        signal, angle = ring_arc(
            nx, ny;
            cx=arc.cx, cy=arc.cy, radius=arc.radius,
            width=2.5, amp=arc.amp,
            phi_min=arc.phi_min, phi_max=arc.phi_max,
        )
        add_component!(cube, truth, signal, finite_angle(signal, angle, 0.0), chan_all(nv); foreground=true)
    end
    return add_noise!(cube, 0.18, rng), truth, 0.18
end

function build_velocity_components(; nx=128, ny=128, nv=24, seed=104)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    for (angle, channels, amp) in ((28 * DEG, chan_range(nv, 1, 10), 3.0),
                                   (108 * DEG, chan_range(nv, 15, 24), 2.6))
        for offset in (-18, 18)
            signal, bangle = gaussian_segment(
                nx, ny;
                cx=64 + offset * cos(angle + pi / 2),
                cy=64 + offset * sin(angle + pi / 2),
                angle=angle, length=58, width=2.6, amp=amp,
            )
            add_component!(cube, truth, signal, bangle, channels; foreground=true)
        end
    end
    return add_noise!(cube, 0.20, rng), truth, 0.20
end

function build_shells_arcs(; nx=128, ny=128, nv=24, seed=105)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    for arc in ((cx=56.0, cy=60.0, radius=28.0, phi_min=-2.8, phi_max=2.2, amp=3.0),
                (cx=92.0, cy=88.0, radius=20.0, phi_min=-3.0, phi_max=0.9, amp=2.2))
        signal, angle = ring_arc(
            nx, ny;
            cx=arc.cx, cy=arc.cy, radius=arc.radius,
            width=2.3, amp=arc.amp,
            phi_min=arc.phi_min, phi_max=arc.phi_max,
        )
        add_component!(cube, truth, signal, finite_angle(signal, angle, 0.0), chan_all(nv); foreground=true)
    end
    return add_noise!(cube, 0.18, rng), truth, 0.18
end

function build_planck_misaligned(; nx=128, ny=128, nv=24, seed=106)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    fg_angle = 35 * DEG
    contaminant_angle = 115 * DEG
    for offset in (-22, 0, 22)
        signal, bangle = gaussian_segment(
            nx, ny;
            cx=64 + offset * cos(fg_angle + pi / 2),
            cy=64 + offset * sin(fg_angle + pi / 2),
            angle=fg_angle, length=62, width=3.0, amp=2.0,
        )
        add_component!(cube, truth, signal, bangle, chan_range(nv, 1, 12); foreground=true)
    end
    for offset in (-14, 14)
        signal, bangle = gaussian_segment(
            nx, ny;
            cx=66 + offset * cos(contaminant_angle + pi / 2),
            cy=66 + offset * sin(contaminant_angle + pi / 2),
            angle=contaminant_angle, length=55, width=3.4, amp=2.7,
        )
        add_component!(cube, truth, signal, bangle, chan_range(nv, 13, 24); foreground=false)
    end
    return add_noise!(cube, 0.20, rng), truth, 0.20
end

function build_sparse_co_like(; nx=128, ny=128, nv=24, seed=107)
    rng = MersenneTwister(seed)
    cube = zeros(nx, ny, nv)
    truth = TruthMaps(nx, ny)
    for n in 1:10
        angle = (20 + 130 * rand(rng)) * DEG
        signal, bangle = gaussian_segment(
            nx, ny;
            cx=rand(rng, 15:nx-15),
            cy=rand(rng, 15:ny-15),
            angle=angle,
            length=24 + 16 * rand(rng),
            width=0.9 + 0.8 * rand(rng),
            amp=3.0 + 2.0 * rand(rng),
            gap_period=10 + 8 * rand(rng),
        )
        channels = chan_range(nv, rand(rng, 1:18), rand(rng, 19:24))
        add_component!(cube, truth, signal, bangle, channels; foreground=true)
    end
    return add_noise!(cube, 0.24, rng), truth, 0.24
end

function estimated_b_angle(q, u)
    out = fill(NaN, size(q))
    @inbounds for idx in eachindex(q, u)
        (q[idx] == 0.0 && u[idx] == 0.0) && continue
        out[idx] = 0.5 * atan(u[idx], q[idx]) + pi / 2
    end
    return out
end

function qtile(vals, q)
    isempty(vals) && return NaN
    return quantile(sort(vals), q)
end

function evaluate_region(est, truth_angle, truth_rbar, signal; region::Symbol)
    vals = Float64[]
    eligible = 0
    signal_cut = qtile([v for v in vec(signal) if v > 0], 0.55)
    @inbounds for idx in eachindex(est, truth_angle, truth_rbar, signal)
        signal[idx] > signal_cut || continue
        isfinite(truth_angle[idx]) || continue
        if region === :clean
            truth_rbar[idx] >= 0.75 || continue
        elseif region === :ambiguous
            truth_rbar[idx] < 0.50 || continue
        elseif region === :all_signal
            # keep all selected signal pixels
        else
            throw(ArgumentError("unknown region $region"))
        end
        eligible += 1
        isfinite(est[idx]) || continue
        push!(vals, abs(rad2deg(headless_angle_difference(est[idx], truth_angle[idx]))))
    end
    return (
        eligible=eligible,
        n=length(vals),
        support=eligible == 0 ? NaN : length(vals) / eligible,
        med=qtile(vals, 0.50),
        p90=qtile(vals, 0.90),
    )
end

function run_mode(cube, noise, mode, fitter, min_valid::Int)
    q, u, profile = VGT_QU_mv_2d_profiled(
        cube, 12, 1, size(cube, 3), noise;
        min_valid=min_valid,
        fitter_mode=fitter,
        gradient_mode=mode,
    )
    est = estimated_b_angle(q, u)
    rbar_vals = [Float64(v) for v in vec(profile[:fit_rbar2]) if isfinite(v)]
    support_vals = [Float64(v) for v in vec(profile[:fit_support]) if v > 0]
    return est, profile, qtile(rbar_vals, 0.50), qtile(support_vals, 0.50)
end

function print_result(scenario, reference, region, mode, fitter, min_valid_label, min_valid, metrics,
                      median_rbar2, median_fit_support, attempts)
    @printf("%-26s %-16s %-12s %-18s %-18s %-8s %5d %7d %8d %8.3f %10.3f %10.3f %10.3f %10.3f %10d\n",
            scenario, reference, region, mode, fitter, min_valid_label, min_valid,
            metrics.n, metrics.eligible, metrics.support, metrics.med, metrics.p90,
            median_rbar2, median_fit_support, attempts)
end

function run_scenario(label, builder)
    cube, truth, noise = builder(; nx=MAP_SIZE, ny=MAP_SIZE)
    total_angle, total_rbar = truth_angle_and_rbar(truth.q, truth.u, truth.w)
    fg_angle, fg_rbar = truth_angle_and_rbar(truth.fg_q, truth.fg_u, truth.fg_w)
    for mode in MODES
        for fitter in FITTERS
            for (min_valid_label, min_valid) in MIN_VALID_CANDIDATES
                est, profile, median_rbar2, median_fit_support = run_mode(cube, noise, mode, fitter, min_valid)
                for region in (:clean, :all_signal, :ambiguous)
                    print_result(label, "total", region, mode, fitter, min_valid_label, min_valid,
                                 evaluate_region(est, total_angle, total_rbar, truth.w; region=region),
                                 median_rbar2, median_fit_support, profile[:fit_attempts])
                end
                if label == "planck_misaligned"
                    for region in (:clean, :all_signal)
                        print_result(label, "foreground", region, mode, fitter, min_valid_label, min_valid,
                                     evaluate_region(est, fg_angle, fg_rbar, truth.fg_w; region=region),
                                     median_rbar2, median_fit_support, profile[:fit_attempts])
                    end
                end
            end
        end
    end
end

function main()
    println("Milky-Way-like synthetic VGT morphology validation")
    min_valid_text = join(["$(label)=$(value)" for (label, value) in MIN_VALID_CANDIDATES], ", ")
    println("Fitters: $(join(string.(FITTERS), ", ")); map=$(MAP_SIZE)x$(MAP_SIZE); dn=12; min_valid candidates: $(min_valid_text)")
    println("Gradient modes: $(join(string.(MODES), ", "))")
    println("Angles are final VGT magnetic-field proxies compared to planted filament tangents.")
    println("Low truth coherence regions are reported separately because crossings/mixes have no unique answer.")
    println()
    @printf("%-26s %-16s %-12s %-18s %-18s %-8s %5s %7s %8s %8s %10s %10s %10s %10s %10s\n",
            "scenario", "reference", "region", "mode", "fitter", "minlbl", "minv", "n", "eligible",
            "support", "med_deg", "p90_deg", "med_rbar2", "med_fit_support", "attempts")
    println("-" ^ 192)
    scenarios = [
        ("many_filaments", build_many_filaments),
        ("crossing_filaments", build_crossing_filaments),
        ("curved_fibers", build_curved_fibers),
        ("velocity_components", build_velocity_components),
        ("shells_arcs", build_shells_arcs),
        ("planck_misaligned", build_planck_misaligned),
        ("sparse_co_like", build_sparse_co_like),
    ]
    for (label, builder) in scenarios
        run_scenario(label, builder)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
