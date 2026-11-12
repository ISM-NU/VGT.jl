using DrWatson
quickactivate(@__DIR__, "VGT")

using Printf
using Random
using Statistics
using VGT: build_noise_cube, fit_window_result, headless_angle_difference

const MAD_NORMAL = 0.6744897501960817
const TRUE_ANGLE_DEG = 37.0
const TRUE_ANGLE = deg2rad(TRUE_ANGLE_DEG)

robust_negative_noise(values) = begin
    vals = [abs(Float64(v)) for v in values if isfinite(v) && v < 0.0]
    length(vals) < 100 && return NaN
    return median(vals) / MAD_NORMAL
end

function make_patchy_mask(nx::Int, ny::Int, fill_fraction::Float64, rng)
    target = clamp(round(Int, fill_fraction * nx * ny), 1, nx * ny)
    mask = falses(nx, ny)
    n_clumps = max(1, round(Int, 8 * fill_fraction))
    base_radius = max(4.0, sqrt(target / (pi * n_clumps)))
    attempts = 0
    while count(mask) < target && attempts < 120 * n_clumps
        attempts += 1
        cx = rand(rng, 1:nx)
        cy = rand(rng, 1:ny)
        radius = base_radius * (0.75 + 0.6 * rand(rng))
        r2 = radius^2
        @inbounds for j in max(1, floor(Int, cy - radius)):min(ny, ceil(Int, cy + radius))
            for i in max(1, floor(Int, cx - radius)):min(nx, ceil(Int, cx + radius))
                (i - cx)^2 + (j - cy)^2 <= r2 && (mask[i, j] = true)
            end
        end
    end
    return mask
end

function sigma_value(kind::Symbol, i::Int, j::Int, k::Int, nx::Int, ny::Int, nv::Int)
    sx = nx == 1 ? 0.0 : (i - 1) / (nx - 1)
    sy = ny == 1 ? 0.0 : (j - 1) / (ny - 1)
    sk = nv == 1 ? 0.0 : (k - 1) / (nv - 1)
    spatial = 0.11 * (1.0 + 1.05 * sx + 0.35 * sy)
    channel = 0.11 * (1.0 + 1.25 * sk)
    if kind === :stationary
        return 0.17
    elseif kind === :channel
        return channel
    elseif kind === :spatial
        return spatial
    elseif kind === :mixed
        return spatial * channel / 0.17
    end
    throw(ArgumentError("unknown synthetic noise kind: $kind"))
end

function make_foreground_cube(kind::Symbol; nx = 72, ny = 72, nv = 18,
                              fill_fraction = 0.45, rng = MersenneTwister(42))
    cube = Array{Float64}(undef, nx, ny, nv)
    sigma = Array{Float64}(undef, nx, ny, nv)
    mask = make_patchy_mask(nx, ny, fill_fraction, rng)
    ic = (nx + 1) / 2.0
    jc = (ny + 1) / 2.0
    gx = cos(TRUE_ANGLE)
    gy = sin(TRUE_ANGLE)

    @inbounds for k in 1:nv
        env = 0.9 + 0.25 * sin(2pi * k / nv)
        for j in 1:ny, i in 1:nx
            sig = sigma_value(kind, i, j, k, nx, ny, nv)
            sigma[i, j, k] = sig
            line = 0.0
            if mask[i, j]
                x = i - ic
                y = j - jc
                line = 0.62 * env + 0.014 * (gx * x + gy * y) * env
            end
            cube[i, j, k] = line + sig * randn(rng)
        end
    end
    return cube, sigma, mask
end

function estimate_scalar(cube)
    return robust_negative_noise(cube)
end

function estimate_channel(cube; fallback = estimate_scalar(cube))
    _, _, nv = size(cube)
    est = fill(fallback, nv)
    for k in 1:nv
        v = robust_negative_noise(@view cube[:, :, k])
        isfinite(v) && (est[k] = v)
    end
    return est
end

function estimate_spatial(cube; tile = 18, fallback = estimate_scalar(cube))
    nx, ny, _ = size(cube)
    est = fill(fallback, nx, ny)
    for y1 in 1:tile:ny, x1 in 1:tile:nx
        x2 = min(nx, x1 + tile - 1)
        y2 = min(ny, y1 + tile - 1)
        v = robust_negative_noise(@view cube[x1:x2, y1:y2, :])
        isfinite(v) || (v = fallback)
        est[x1:x2, y1:y2] .= v
    end
    return est
end

function build_thresholds(cube, true_sigma)
    scalar = estimate_scalar(cube)
    return Dict(
        :scalar => fill(scalar, size(cube)),
        :channel => build_noise_cube(cube, scalar, :channel; tile_size = 18),
        :spatial => build_noise_cube(cube, scalar, :spatial; tile_size = 18),
        :separable => build_noise_cube(cube, scalar, :separable; tile_size = 18),
        :oracle => true_sigma,
    )
end

function threshold_error(est_sigma, true_sigma)
    rel = Float64[]
    @inbounds for idx in eachindex(est_sigma, true_sigma)
        t = true_sigma[idx]
        t > 0 || continue
        push!(rel, abs(est_sigma[idx] - t) / t)
    end
    return median(rel), quantile(rel, 0.90)
end

function sampled_pixel_psi(cube, sigma_est, ii::Int, jj::Int;
                           dn = 12, min_valid = 72, fitter = :axial_vonmises)
    nx, ny, nv = size(cube)
    half_dn = dn ÷ 2
    ii < half_dn && return NaN
    jj < half_dn && return NaN
    ii > nx - half_dn && return NaN
    jj > ny - half_dn && return NaN

    is = ii - half_dn + 1
    ie = ii + half_dn
    js = jj - half_dn + 1
    je = jj + half_dn
    q = 0.0
    u = 0.0
    successes = 0
    ax = Matrix{Float64}(undef, dn, dn)
    ay = Matrix{Float64}(undef, dn, dn)

    @inbounds for k in 1:nv
        valid = 0
        for (lj, j) in enumerate(js:je), (li, i) in enumerate(is:ie)
            val = cube[i, j, k]
            thresh = sigma_est[i, j, k]
            if i == 1 || j == 1 || !isfinite(val) || val < thresh
                ax[li, lj] = NaN
                ay[li, lj] = NaN
            else
                ax[li, lj] = cube[i - 1, j, k] - val
                ay[li, lj] = cube[i, j - 1, k] - val
                valid += 1
            end
        end
        valid >= min_valid || continue

        result = fit_window_result(fitter, ax, ay)
        result.fit_ok || continue
        weight = max(cube[ii, jj, k], 0.0)
        q = muladd(weight, cos(2.0 * result.angle), q)
        u = muladd(weight, sin(2.0 * result.angle), u)
        successes += 1
    end

    return successes == 0 || (q == 0.0 && u == 0.0) ? NaN : 0.5 * atan(u, q)
end

function evaluate_strategy(cube, sigma_est; dn = 12, stride = 9)
    nx, ny, _ = size(cube)
    half_dn = dn ÷ 2
    errs = Float64[]
    attempted = 0
    for j in (half_dn + 2):stride:(ny - half_dn - 2), i in (half_dn + 2):stride:(nx - half_dn - 2)
        attempted += 1
        psi = sampled_pixel_psi(cube, sigma_est, i, j; dn = dn)
        isfinite(psi) || continue
        push!(errs, abs(rad2deg(headless_angle_difference(psi, TRUE_ANGLE))))
    end
    return (
        median_error = isempty(errs) ? NaN : median(errs),
        p90_error = isempty(errs) ? NaN : quantile(errs, 0.90),
        valid_fraction = attempted == 0 ? NaN : length(errs) / attempted,
        n = length(errs),
        attempted = attempted,
    )
end

function run_scenario_once(kind::Symbol, seed::Int)
    cube, true_sigma, _ = make_foreground_cube(kind; rng = MersenneTwister(seed))
    estimates = build_thresholds(cube, true_sigma)
    order = [:scalar, :channel, :spatial, :separable, :oracle]
    rows = []
    for name in order
        med_rel, p90_rel = threshold_error(estimates[name], true_sigma)
        metrics = evaluate_strategy(cube, estimates[name])
        push!(rows, (name, med_rel, p90_rel, metrics))
    end
    return rows
end

function aggregate_metric(values)
    clean = collect(filter(isfinite, values))
    isempty(clean) && return NaN
    return median(clean)
end

function run_scenario(kind::Symbol; seeds = 7301:7307)
    order = [:scalar, :channel, :spatial, :separable, :oracle]
    raw = Dict(name => [] for name in order)
    for seed in seeds
        for row in run_scenario_once(kind, seed)
            push!(raw[row[1]], row)
        end
    end

    rows = []
    for name in order
        rs = raw[name]
        med_rel = aggregate_metric([r[2] for r in rs])
        p90_rel = aggregate_metric([r[3] for r in rs])
        med_err = aggregate_metric([r[4].median_error for r in rs])
        p90_err = aggregate_metric([r[4].p90_error for r in rs])
        valid_fraction = aggregate_metric([r[4].valid_fraction for r in rs])
        n = sum(r[4].n for r in rs)
        attempted = sum(r[4].attempted for r in rs)
        metrics = (
            median_error = med_err,
            p90_error = p90_err,
            valid_fraction = valid_fraction,
            n = n,
            attempted = attempted,
        )
        push!(rows, (name, med_rel, p90_rel, metrics))
    end
    return rows
end

function print_scenario(kind::Symbol, rows)
    println()
    println("=" ^ 86)
    println("Noise model: $kind")
    println("=" ^ 86)
    @printf("%-11s  %-13s  %-13s  %-12s  %-12s  %-9s\n",
            "strategy", "sigma_med%", "sigma_p90%", "psi_med_deg", "psi_p90_deg", "valid")
    println("-" ^ 86)
    for (name, med_rel, p90_rel, metrics) in rows
        @printf("%-11s  %-13.2f  %-13.2f  %-12.3f  %-12.3f  %3d/%-3d\n",
                String(name), 100 * med_rel, 100 * p90_rel,
                metrics.median_error, metrics.p90_error, metrics.n, metrics.attempted)
    end
end

function winner(rows)
    candidates = filter(r -> r[1] !== :oracle && isfinite(r[4].median_error), rows)
    isempty(candidates) && return :none
    best_med = minimum(r[4].median_error for r in candidates)
    best_p90 = minimum(r[4].p90_error for r in candidates)
    tied = filter(r -> r[4].median_error <= best_med + 0.5 &&
                       r[4].p90_error <= best_p90 + 1.0, candidates)
    complexity = Dict(
        :scalar => 1,
        :channel => 2,
        :spatial => 2,
        :separable => 4,
    )
    sort!(tied, by = r -> (r[2], complexity[r[1]], -r[4].valid_fraction))
    return tied[1][1]
end

function main()
    println("Synthetic noise-threshold evaluation for VGT foreground-style cubes")
    println("True field angle: $(TRUE_ANGLE_DEG) deg; fitter: axial_vonmises; threshold = noise")
    println("Each row is the median over seven deterministic synthetic realizations.")
    println("Oracle is included as an upper-bound diagnostic, not an available estimator.")

    kinds = [:stationary, :channel, :spatial, :mixed]
    wins = Dict{Symbol, Int}()
    for kind in kinds
        rows = run_scenario(kind)
        print_scenario(kind, rows)
        w = winner(rows)
        wins[w] = get(wins, w, 0) + 1
        println("Best available strategy for $kind: $w")
    end

    println()
    println("=" ^ 86)
    println("Summary")
    println("=" ^ 86)
    for name in [:scalar, :channel, :spatial, :separable, :none]
        get(wins, name, 0) > 0 && println("  $name wins: $(wins[name]) / $(length(kinds))")
    end
    println()
    println("Interpretation:")
    println("  Use scalar noise when channel and spatial spread are small.")
    println("  Use per-channel noise when channel spread dominates.")
    println("  Use coarse spatial noise when spatial spread dominates.")
    println("  Use full separable spatial-channel noise only when both spreads are large and validated.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
