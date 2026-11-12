using DrWatson
quickactivate(@__DIR__, "VGT")

using Printf
using Random
using Statistics
using VGT: VGT_QU_mv_2d_profiled, fit_window_result, headless_angle_difference

# ---------------------------------------------------------------------------
# Synthetic cube design:
#   cube[i,j,k] = BASE + A*(cos(θ)*(i-ic) + sin(θ)*(j-jc))*env(k) + noise
#
# Backward-difference gradients (gradient_mode = :nan_edge_backward):
#   grad_x[i,j] = cube[i-1,j,k] - cube[i,j,k] ≈ -A*cos(θ)*env
#   grad_y[i,j] = cube[i,j-1,k] - cube[i,j,k] ≈ -A*sin(θ)*env
#
# Gradient angle = atan(-sin(θ), -cos(θ)) = θ+π → doubled-angle axial trick → θ ✓
#
# SNR_grad = A / (noise * √2):
#   noise=0.01 → SNR≈10   noise=0.15 → SNR≈0.7   noise=0.50 → SNR≈0.2
# BASE = 20 → all pixels far above pipeline thresh = 3*1.0 = 3.
# ---------------------------------------------------------------------------

const TRUE_ANGLE_DEG = 37.0
const TRUE_ANGLE     = deg2rad(TRUE_ANGLE_DEG)

# ─── Cube constructors ──────────────────────────────────────────────────────

function make_uniform_cube(;
        nx = 64, ny = 64, nv = 32,
        noise_level = 0.05, gradient_amplitude = 0.15,
        angle = TRUE_ANGLE, rng = MersenneTwister(42))
    BASE = 20.0; ic = (nx + 1) / 2.0; jc = (ny + 1) / 2.0
    gx = cos(angle); gy = sin(angle)
    cube = Array{Float64}(undef, nx, ny, nv)
    @inbounds for k in 1:nv
        env = 1.0 + 0.4 * sin(2π * k / nv)
        for j in 1:ny, i in 1:nx
            cube[i, j, k] = BASE + gradient_amplitude * (gx * (i - ic) + gy * (j - jc)) * env +
                             noise_level * randn(rng)
        end
    end
    return cube
end

function make_patchy_mask(nx::Int, ny::Int, fill_fraction::Float64, rng)
    target = clamp(round(Int, fill_fraction * nx * ny), 1, nx * ny)
    mask = falses(nx, ny)
    n_clumps = max(1, round(Int, 10 * fill_fraction))
    base_radius = max(3.0, sqrt(target / (π * n_clumps)))
    attempts = 0
    while count(mask) < target && attempts < 100 * n_clumps
        attempts += 1
        cx = rand(rng, 1:nx)
        cy = rand(rng, 1:ny)
        radius = base_radius * (0.75 + 0.5 * rand(rng))
        r2 = radius^2
        @inbounds for j in max(1, floor(Int, cy - radius)):min(ny, ceil(Int, cy + radius))
            for i in max(1, floor(Int, cx - radius)):min(nx, ceil(Int, cx + radius))
                (i - cx)^2 + (j - cy)^2 <= r2 && (mask[i, j] = true)
            end
        end
    end
    return mask
end

# Approximate fraction fill_fraction of pixels retain signal in coherent patches;
# the rest are zeroed below threshold. Patchy masks avoid the single-pixel edge
# gradients produced by salt-and-pepper random dropout.
function make_sparse_cube(;
        nx = 128, ny = 128, nv = 32, noise_level = 0.05,
        fill_fraction = 0.5, rng = MersenneTwister(42))
    cube = make_uniform_cube(; nx = nx, ny = ny, nv = nv, noise_level = noise_level, rng = rng)
    mask = make_patchy_mask(nx, ny, fill_fraction, rng)
    @inbounds for k in 1:nv, j in 1:ny, i in 1:nx
        mask[i, j] || (cube[i, j, k] = 0.0)
    end
    return cube
end

# Curved scalar field with a known finite-difference gradient orientation.
# This avoids scoring against an arbitrary angle field that is not actually the
# gradient of the planted intensity cube.
function make_varying_angle_cube(;
        nx = 128, ny = 128, nv = 32, noise_level = 0.05,
        gradient_amplitude = 0.15, curvature = 8e-4,
        rng = MersenneTwister(42))
    BASE = 20.0; ic = (nx + 1) / 2.0; jc = (ny + 1) / 2.0
    gx = cos(TRUE_ANGLE); gy = sin(TRUE_ANGLE)
    cube = Array{Float64}(undef, nx, ny, nv)
    @inbounds for k in 1:nv
        env = 1.0 + 0.4 * sin(2π * k / nv)
        for j in 1:ny, i in 1:nx
            x = Float64(i) - ic
            y = Float64(j) - jc
            signal = gradient_amplitude * (gx * x + gy * y) + curvature * x * y
            cube[i, j, k] = BASE + signal * env + noise_level * randn(rng)
        end
    end
    return cube
end

function varying_true_angle(i, j, nx, ny; gradient_amplitude = 0.15, curvature = 8e-4)
    ic = (nx + 1) / 2.0
    jc = (ny + 1) / 2.0
    gx = gradient_amplitude * cos(TRUE_ANGLE) + curvature * (Float64(j) - jc)
    gy = gradient_amplitude * sin(TRUE_ANGLE) + curvature * (Float64(i) - ic)
    return atan(gy, gx)
end

# Left half (i ≤ nx/2) has TRUE_ANGLE; right half has TRUE_ANGLE + delta_deg.
# Windows straddling the seam see a bimodal gradient distribution.
function make_boundary_cube(;
        nx = 128, ny = 128, nv = 32, noise_level = 0.05,
        delta_deg = 90.0, rng = MersenneTwister(42))
    BASE = 20.0; ic = (nx + 1) / 2.0; jc = (ny + 1) / 2.0
    θ1 = TRUE_ANGLE; θ2 = TRUE_ANGLE + deg2rad(delta_deg)
    cube = Array{Float64}(undef, nx, ny, nv)
    @inbounds for k in 1:nv
        env = 1.0 + 0.4 * sin(2π * k / nv)
        for j in 1:ny, i in 1:nx
            θ = i <= nx ÷ 2 ? θ1 : θ2
            cube[i, j, k] = BASE + 0.15 * (cos(θ) * (i - ic) + sin(θ) * (j - jc)) * env +
                             noise_level * randn(rng)
        end
    end
    return cube
end

# Axial circular mean of two equal-weight headless orientations.
function circular_mean_2(θ1, θ2)
    s = sin(2θ1) + sin(2θ2)
    c = cos(2θ1) + cos(2θ2)
    return 0.5 * atan(s, c)
end

# ─── Measurement helpers ────────────────────────────────────────────────────

# Suppresses per-channel progress lines from the VGT implementation.
function run_fitter(cube, fitter; dn = 10, min_valid = 50, pipeline_noise = 1.0)
    nv = size(cube, 3)
    local q, u, profile
    redirect_stdout(devnull) do
        q, u, profile = VGT_QU_mv_2d_profiled(
            cube, dn, 1, nv, pipeline_noise;
            min_valid = min_valid, fitter_mode = fitter)
    end
    return q, u, profile
end

# Returns per-pixel |Δψ| in degrees; skips pixels with Q=U=0 (not computed).
function collect_errors(q, u, angle_fn)
    errs = Float64[]
    nx, ny = size(q)
    @inbounds for j in 1:ny, i in 1:nx
        (q[i, j] == 0.0 && u[i, j] == 0.0) && continue
        psi = 0.5 * atan(u[i, j], q[i, j])
        e = abs(rad2deg(headless_angle_difference(psi, angle_fn(i, j))))
        isfinite(e) && push!(errs, e)
    end
    return errs
end

med(v) = isempty(v) ? NaN : median(v)
noisee(v) = isempty(v) ? NaN : sqrt(mean(v .^ 2))

function meanfinite(v)
    finite = filter(isfinite, v)
    return isempty(finite) ? NaN : mean(finite)
end

function quantilefinite(v, p)
    finite = sort(filter(isfinite, v))
    return isempty(finite) ? NaN : quantile(finite, p)
end

# Single-pixel evaluation (backward-difference gradients) for targeted stress tests.
# It uses the same backward-difference gradients and fit_window_result call as
# VGT_QU_mv_2d_profiled, but avoids running a full map when only a few windows
# are scientifically relevant.
function sampled_pixel_psi(cube, ii::Int, jj::Int, fitter::Symbol;
                           dn = 10, pipeline_noise = 1.0, min_valid = 50)
    return sampled_pixel_psi_variant(
        cube, ii, jj, fitter;
        dn = dn,
        pipeline_noise = pipeline_noise,
        min_valid = min_valid,
        rbar2_weight = false,
        min_rbar2 = 0.0,
    )
end

function sampled_pixel_psi_variant(cube, ii::Int, jj::Int, fitter::Symbol;
                                   dn = 10, pipeline_noise = 1.0,
                                   min_valid = 50, rbar2_weight = false,
                                   min_rbar2 = 0.0)
    nx, ny, nv = size(cube)
    half_dn = dn ÷ 2
    if ii < half_dn || ii > nx - half_dn || jj < half_dn || jj > ny - half_dn
        return NaN
    end

    is = ii - half_dn + 1
    ie = ii + half_dn
    js = jj - half_dn + 1
    je = jj + half_dn
    thresh = pipeline_noise
    q = 0.0
    u = 0.0
    successes = 0
    grad_x_window = Matrix{Float64}(undef, dn, dn)
    grad_y_window = Matrix{Float64}(undef, dn, dn)

    @inbounds for k in 1:nv
        valid = 0
        for (lj, j) in enumerate(js:je), (li, i) in enumerate(is:ie)
            val = cube[i, j, k]
            if i == 1 || j == 1 || val < thresh
                grad_x_window[li, lj] = NaN
                grad_y_window[li, lj] = NaN
            else
                grad_x_window[li, lj] = cube[i - 1, j, k] - val
                grad_y_window[li, lj] = cube[i, j - 1, k] - val
                valid += 1
            end
        end
        valid >= min_valid || continue

        result = fit_window_result(fitter, grad_x_window, grad_y_window)
        result.fit_ok || continue
        if min_rbar2 > 0.0 && (!isfinite(result.rbar2) || result.rbar2 < min_rbar2)
            continue
        end
        weight = cube[ii, jj, k]
        if rbar2_weight
            isfinite(result.rbar2) || continue
            weight *= result.rbar2
        end
        q = muladd(weight, cos(2.0 * result.angle), q)
        u = muladd(weight, sin(2.0 * result.angle), u)
        successes += 1
    end

    return (successes == 0 || (q == 0.0 && u == 0.0)) ? NaN : 0.5 * atan(u, q)
end

function sampled_centers(nx::Int, ny::Int; dn = 10, stride = 8)
    half_dn = dn ÷ 2
    is = collect((half_dn + 2):stride:(nx - half_dn - 2))
    js = collect((half_dn + 2):stride:(ny - half_dn - 2))
    isempty(is) && (is = [nx ÷ 2])
    isempty(js) && (js = [ny ÷ 2])
    return is, js
end

function collect_sampled_errors(cube, fitter::Symbol, angle_fn;
                                dn = 10, min_valid = 50, stride = 8)
    nx, ny, _ = size(cube)
    is, js = sampled_centers(nx, ny; dn = dn, stride = stride)
    errs = Float64[]
    @inbounds for j in js, i in is
        psi = sampled_pixel_psi(cube, i, j, fitter; dn = dn, min_valid = min_valid)
        isfinite(psi) || continue
        push!(errs, abs(rad2deg(headless_angle_difference(psi, angle_fn(i, j)))))
    end
    return errs
end

function collect_sampled_errors_variant(cube, fitter::Symbol, angle_fn;
                                        dn = 10, min_valid = 50, stride = 8,
                                        rbar2_weight = false, min_rbar2 = 0.0)
    nx, ny, _ = size(cube)
    is, js = sampled_centers(nx, ny; dn = dn, stride = stride)
    errs = Float64[]
    @inbounds for j in js, i in is
        psi = sampled_pixel_psi_variant(
            cube, i, j, fitter;
            dn = dn,
            min_valid = min_valid,
            rbar2_weight = rbar2_weight,
            min_rbar2 = min_rbar2,
        )
        isfinite(psi) || continue
        push!(errs, abs(rad2deg(headless_angle_difference(psi, angle_fn(i, j)))))
    end
    return errs
end

function hanning3_velocity(cube)
    nx, ny, nv = size(cube)
    out = similar(cube, Float64)
    @inbounds for k in 1:nv
        if k == 1
            out[:, :, k] .= 0.75 .* cube[:, :, k] .+ 0.25 .* cube[:, :, k + 1]
        elseif k == nv
            out[:, :, k] .= 0.25 .* cube[:, :, k - 1] .+ 0.75 .* cube[:, :, k]
        else
            out[:, :, k] .= 0.25 .* cube[:, :, k - 1] .+ 0.50 .* cube[:, :, k] .+
                             0.25 .* cube[:, :, k + 1]
        end
    end
    return out
end

@inline above_threshold(cube, i::Int, j::Int, k::Int, thresh::Real) =
    isfinite(cube[i, j, k]) && cube[i, j, k] >= thresh

function experimental_gradient_pair(cube, i::Int, j::Int, k::Int,
                                    thresh::Real, mode::Symbol)
    nx, ny, _ = size(cube)
    if mode === :backward_current
        (i == 1 || j == 1 || !above_threshold(cube, i, j, k, thresh)) &&
            return NaN, NaN
        return cube[i - 1, j, k] - cube[i, j, k],
               cube[i, j - 1, k] - cube[i, j, k]
    elseif mode === :backward_checked
        (i == 1 || j == 1) && return NaN, NaN
        (above_threshold(cube, i, j, k, thresh) &&
         above_threshold(cube, i - 1, j, k, thresh) &&
         above_threshold(cube, i, j - 1, k, thresh)) || return NaN, NaN
        return cube[i - 1, j, k] - cube[i, j, k],
               cube[i, j - 1, k] - cube[i, j, k]
    elseif mode === :central_checked
        (i == 1 || i == nx || j == 1 || j == ny) && return NaN, NaN
        (above_threshold(cube, i, j, k, thresh) &&
         above_threshold(cube, i - 1, j, k, thresh) &&
         above_threshold(cube, i + 1, j, k, thresh) &&
         above_threshold(cube, i, j - 1, k, thresh) &&
         above_threshold(cube, i, j + 1, k, thresh)) || return NaN, NaN
        return (cube[i - 1, j, k] - cube[i + 1, j, k]) / 2,
               (cube[i, j - 1, k] - cube[i, j + 1, k]) / 2
    elseif mode === :sobel_checked
        (i == 1 || i == nx || j == 1 || j == ny) && return NaN, NaN
        for jj in (j - 1):(j + 1), ii in (i - 1):(i + 1)
            above_threshold(cube, ii, jj, k, thresh) || return NaN, NaN
        end
        gx = (
            cube[i - 1, j - 1, k] + 2 * cube[i - 1, j, k] + cube[i - 1, j + 1, k] -
            cube[i + 1, j - 1, k] - 2 * cube[i + 1, j, k] - cube[i + 1, j + 1, k]
        ) / 8
        gy = (
            cube[i - 1, j - 1, k] + 2 * cube[i, j - 1, k] + cube[i + 1, j - 1, k] -
            cube[i - 1, j + 1, k] - 2 * cube[i, j + 1, k] - cube[i + 1, j + 1, k]
        ) / 8
        return gx, gy
    end
    throw(ArgumentError("unsupported experimental gradient mode $mode"))
end

function sampled_pixel_psi_experimental(cube, ii::Int, jj::Int, fitter::Symbol;
                                        dn = 10, pipeline_noise = 1.0,
                                        min_valid = 50,
                                        gradient_mode = :backward_current,
                                        rbar2_weight = false,
                                        min_rbar2 = 0.0)
    nx, ny, nv = size(cube)
    half_dn = dn ÷ 2
    if ii < half_dn || ii > nx - half_dn || jj < half_dn || jj > ny - half_dn
        return NaN
    end

    is = ii - half_dn + 1
    ie = ii + half_dn
    js = jj - half_dn + 1
    je = jj + half_dn
    q = 0.0
    u = 0.0
    successes = 0
    grad_x_window = Matrix{Float64}(undef, dn, dn)
    grad_y_window = Matrix{Float64}(undef, dn, dn)

    @inbounds for k in 1:nv
        valid = 0
        for (lj, j) in enumerate(js:je), (li, i) in enumerate(is:ie)
            gx, gy = experimental_gradient_pair(cube, i, j, k, pipeline_noise, gradient_mode)
            grad_x_window[li, lj] = gx
            grad_y_window[li, lj] = gy
            isfinite(gx) && isfinite(gy) && (valid += 1)
        end
        valid >= min_valid || continue

        result = fit_window_result(fitter, grad_x_window, grad_y_window)
        result.fit_ok || continue
        if min_rbar2 > 0.0 && (!isfinite(result.rbar2) || result.rbar2 < min_rbar2)
            continue
        end
        weight = cube[ii, jj, k]
        if rbar2_weight
            isfinite(result.rbar2) || continue
            weight *= result.rbar2
        end
        q = muladd(weight, cos(2.0 * result.angle), q)
        u = muladd(weight, sin(2.0 * result.angle), u)
        successes += 1
    end

    return (successes == 0 || (q == 0.0 && u == 0.0)) ? NaN : 0.5 * atan(u, q)
end

function collect_sampled_errors_experimental(cube, angle_fn;
                                             dn = 10, min_valid = 50, stride = 8,
                                             gradient_mode = :backward_current,
                                             rbar2_weight = false,
                                             min_rbar2 = 0.0,
                                             spectral_hanning = false)
    work = spectral_hanning ? hanning3_velocity(cube) : cube
    nx, ny, _ = size(work)
    is, js = sampled_centers(nx, ny; dn = dn, stride = stride)
    errs = Float64[]
    @inbounds for j in js, i in is
        psi = sampled_pixel_psi_experimental(
            work, i, j, :axial_vonmises;
            dn = dn,
            min_valid = min_valid,
            gradient_mode = gradient_mode,
            rbar2_weight = rbar2_weight,
            min_rbar2 = min_rbar2,
        )
        isfinite(psi) || continue
        push!(errs, abs(rad2deg(headless_angle_difference(psi, angle_fn(i, j)))))
    end
    return errs
end

function error_summary(errs)
    return (
        n = length(errs),
        med = med(errs),
        p90 = quantilefinite(errs, 0.90),
        rms = noisee(errs),
    )
end

# ─── Test 1: Noise sweep ────────────────────────────────────────────────────

function test_noise_sweep()
    println("=" ^ 72)
    println("Test 1: Noise sweep")
    println("=" ^ 72)
    println("Uniform angle=$(TRUE_ANGLE_DEG)°, A=0.15, nx=64, ny=64, nv=32, dn=10, min_valid=50")
    println("SNR_grad = 0.15 / (noise·√2)")
    println()
    @printf("%-7s  %-10s  %-10s  %-8s\n", "noise", "med°", "snr_grad", "n")
    println("-" ^ 40)

    noise_levels = [0.01, 0.02, 0.05, 0.08, 0.12, 0.15, 0.20, 0.30, 0.40, 0.50]
    meds = Float64[]
    for nl in noise_levels
        cube = make_uniform_cube(; noise_level = nl)
        ea = collect_sampled_errors(cube, :axial_vonmises, (i, j) -> TRUE_ANGLE)
        push!(meds, med(ea))
        @printf("%-7.3f  %-10.3f  %-10.2f  %-8d\n", nl, last(meds), 0.15 / (nl * sqrt(2)), length(ea))
    end
    return (name = "noise sweep", med = meds)
end

# ─── Test 2: Sparse fill sweep ──────────────────────────────────────────────

function test_sparse_fill()
    println()
    println("=" ^ 72)
    println("Test 2: Sparse fill sweep  [patchy emission stress test]")
    println("=" ^ 72)
    println("Uniform angle=$(TRUE_ANGLE_DEG)°, noise=0.05, nx=128, ny=128, nv=32, dn=10, min_valid=50")
    println("fill_fraction: approximate spatial fill in coherent patches (rest zeroed below threshold)")
    println("Low-fill cases may have no valid windows; this is part of the stress test.")
    println()
    @printf("%-8s  %-9s  %-10s  %-8s\n", "fill", "min_v", "med°", "n")
    println("-" ^ 40)

    fill_fractions = [0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.00]
    dn = 10
    meds = Float64[]
    for ff in fill_fractions
        rng  = MersenneTwister(42)
        cube = ff >= 1.0 ?
            make_uniform_cube(; nx = 128, ny = 128, noise_level = 0.05, rng = rng) :
            make_sparse_cube(; fill_fraction = ff, noise_level = 0.05, rng = rng)
        mv = 50
        ea = collect_sampled_errors(cube, :axial_vonmises, (i, j) -> TRUE_ANGLE;
                                    dn = dn, min_valid = mv, stride = 12)
        push!(meds, med(ea))
        @printf("%-8.2f  %-9d  %-10.3f  %-8d\n", ff, mv, last(meds), length(ea))
    end
    return (name = "sparse fill", med = meds)
end

# ─── Test 3: Gradient amplitude sweep ───────────────────────────────────────

function test_amplitude_sweep()
    println()
    println("=" ^ 72)
    println("Test 3: Gradient amplitude sweep")
    println("=" ^ 72)
    println("Uniform angle=$(TRUE_ANGLE_DEG)°, noise=0.05, nx=64, ny=64, nv=32, dn=10, min_valid=50")
    println("SNR_grad = amplitude / (0.05·√2)")
    println()
    @printf("%-8s  %-10s  %-10s  %-8s\n", "amp", "med°", "snr_grad", "n")
    println("-" ^ 40)

    amplitudes = [0.03, 0.05, 0.08, 0.12, 0.15, 0.25, 0.40, 0.60]
    meds = Float64[]
    for amp in amplitudes
        cube = make_uniform_cube(; gradient_amplitude = amp, noise_level = 0.05)
        ea = collect_sampled_errors(cube, :axial_vonmises, (i, j) -> TRUE_ANGLE)
        push!(meds, med(ea))
        @printf("%-8.3f  %-10.3f  %-10.2f  %-8d\n", amp, last(meds), amp / (0.05 * sqrt(2)), length(ea))
    end
    return (name = "amplitude sweep", med = meds)
end

# ─── Test 4: Spatially varying angle ────────────────────────────────────────

function test_varying_angle()
    println()
    println("=" ^ 72)
    println("Test 4: Spatially varying angle  [spatial structure test]")
    println("=" ^ 72)
    println("True angle is the exact finite-difference gradient orientation of a curved scalar field")
    println("nx=128, ny=128, nv=32, dn=10, min_valid=50")
    println("Each dn=10 window spans a few degrees of angle variation in the curved field")
    println()
    @printf("%-7s  %-12s  %-12s  %-8s\n", "noise", "rms°", "med°", "n")
    println("-" ^ 44)

    noise_levels = [0.02, 0.05, 0.10, 0.20]
    nx = 128
    ny = 128
    meds = Float64[]
    for nl in noise_levels
        cube = make_varying_angle_cube(; noise_level = nl)
        ea = collect_sampled_errors(
            cube, :axial_vonmises, (i, j) -> varying_true_angle(i, j, nx, ny);
            min_valid = 50, stride = 12)
        push!(meds, med(ea))
        @printf("%-7.3f  %-12.3f  %-12.3f  %-8d\n", nl, noisee(ea), last(meds), length(ea))
    end
    return (name = "varying angle", med = meds)
end

# ─── Test 5: Boundary / bimodal window test ──────────────────────────────────

function test_boundary()
    println()
    println("=" ^ 72)
    println("Test 5: Two-region boundary  [bimodal window test]")
    println("=" ^ 72)
    println("Left half (i ≤ 64): TRUE_ANGLE=37°.  Right half: 37°+delta.")
    println("nx=128, ny=128, nv=32, dn=10, noise=0.05, min_valid=50")
    println()
    println("Interior (|i-64| > dn+4 = 14): accuracy vs local true angle.")
    println("  Boundary window (|i-64| ≤ dn): accuracy vs circular mean of both populations.")
    println("  Circular mean of θ and θ+δ flips at δ=90° (axial geometry); error spike expected there.")
    println()

    delta_degs = [0, 15, 30, 45, 60, 90]
    nx = 128; ny = 128; dn = 10
    y_samples = [32, 64, 96]
    interior_i_samples = [32, 96]
    boundary_i_samples = [60, 64, 68]

    rows = Tuple[]
    for dd in delta_degs
        cube = make_boundary_cube(; nx = nx, ny = ny, noise_level = 0.05,
                                    delta_deg = Float64(dd))

        θ1 = TRUE_ANGLE; θ2 = TRUE_ANGLE + deg2rad(Float64(dd))
        θ_bnd = circular_mean_2(θ1, θ2)

        int_e = Float64[]
        bnd_e = Float64[]

        @inbounds for j in y_samples, i in interior_i_samples
            θ_true = i <= nx ÷ 2 ? θ1 : θ2
            psi = sampled_pixel_psi(cube, i, j, :axial_vonmises; min_valid = 50)
            isfinite(psi) && push!(int_e, abs(rad2deg(headless_angle_difference(psi, θ_true))))
        end
        @inbounds for j in y_samples, i in boundary_i_samples
            psi = sampled_pixel_psi(cube, i, j, :axial_vonmises; min_valid = 50)
            isfinite(psi) && push!(bnd_e, abs(rad2deg(headless_angle_difference(psi, θ_bnd))))
        end
        push!(rows, (dd, med(int_e), med(bnd_e)))
    end

    println("Interior accuracy and boundary-window mean-orientation behavior")
    println("(boundary reference = circular mean of the two populations):")
    @printf("  %-8s  %-14s  %-14s\n", "delta°", "interior_med°", "boundary_med°")
    println("  " * "-" ^ 40)
    for (dd, mi, mb) in rows
        @printf("  %-8d  %-14.3f  %-14.3f\n", dd, mi, mb)
    end
    return (name = "boundary", med = Float64[r[2] for r in rows])
end

# ─── Test 6: Accuracy-work experiments ──────────────────────────────────────

function candidate_errors(cube, angle_fn; fitter = :axial_vonmises, stride = 12,
                          min_valid = 50, rbar2_weight = false, min_rbar2 = 0.0)
    return collect_sampled_errors_variant(
        cube, fitter, angle_fn;
        stride = stride,
        min_valid = min_valid,
        rbar2_weight = rbar2_weight,
        min_rbar2 = min_rbar2,
    )
end

function test_accuracy_candidates()
    println()
    println("=" ^ 72)
    println("Test 6: Accuracy-work candidates  [do not assume adoption]")
    println("=" ^ 72)
    println("Candidates are evaluated against the current axial_vonmises baseline.")
    println("A candidate must reduce median error without severe coverage loss.")
    println()
    @printf("%-24s  %-18s  %-10s  %-10s  %-8s  %-8s  %s\n",
            "scenario", "candidate", "base_med°", "cand_med°", "base_n", "cand_n", "verdict")
    println("-" ^ 96)

    seeds = 201:206
    scenarios = [
        ("uniform noisy", seed -> make_uniform_cube(; noise_level = 0.12, rng = MersenneTwister(seed)),
         (i, j) -> TRUE_ANGLE, 8),
        ("sparse fill", seed -> make_sparse_cube(; fill_fraction = 0.35, noise_level = 0.05, rng = MersenneTwister(seed)),
         (i, j) -> TRUE_ANGLE, 12),
        ("varying angle", seed -> make_varying_angle_cube(; noise_level = 0.05, rng = MersenneTwister(seed)),
         (i, j) -> varying_true_angle(i, j, 128, 128), 12),
    ]
    candidates = [
        ("rbar2_weight", (; min_valid = 50, rbar2_weight = true, min_rbar2 = 0.0)),
        ("rbar2_gate_0.50", (; min_valid = 50, rbar2_weight = false, min_rbar2 = 0.50)),
    ]

    adopted = Tuple[]
    rejected = Tuple[]
    for (scenario_label, cube_fn, angle_fn, stride) in scenarios
        base_all = Float64[]
        cand_all = Dict(label => Float64[] for (label, _) in candidates)
        cand_counts = Dict(label => 0 for (label, _) in candidates)
        base_counts = 0
        for seed in seeds
            cube = cube_fn(seed)
            base = candidate_errors(cube, angle_fn; stride = stride)
            append!(base_all, base)
            base_counts += length(base)
            for (label, kwargs) in candidates
                cand = candidate_errors(
                    cube, angle_fn;
                    stride = stride,
                    min_valid = kwargs.min_valid,
                    rbar2_weight = kwargs.rbar2_weight,
                    min_rbar2 = kwargs.min_rbar2,
                )
                append!(cand_all[label], cand)
                cand_counts[label] += length(cand)
            end
        end

        base_med = med(base_all)
        for (label, _) in candidates
            cand_med = med(cand_all[label])
            coverage_ratio = base_counts == 0 ? 0.0 : cand_counts[label] / base_counts
            improved = isfinite(base_med) && isfinite(cand_med) && cand_med < base_med
            acceptable_coverage = coverage_ratio >= 0.85
            verdict = improved && acceptable_coverage ? "keep investigating" :
                (improved ? "reject: coverage loss" : "reject: no accuracy gain")
            @printf("%-24s  %-18s  %-10.3f  %-10.3f  %-8d  %-8d  %s\n",
                    scenario_label, label, base_med, cand_med,
                    base_counts, cand_counts[label], verdict)
            if verdict == "keep investigating"
                push!(adopted, (scenario_label, label, base_med, cand_med, base_counts, cand_counts[label]))
            else
                push!(rejected, (scenario_label, label, base_med, cand_med, base_counts, cand_counts[label], verdict))
            end
        end
    end

    println()
    println("Accuracy-candidate policy: none of these changes production behavior.")
    println("Only candidates marked \"keep investigating\" deserve follow-up validation.")
    return (adopted = adopted, rejected = rejected)
end

# ─── Test 7: Method improvement search ──────────────────────────────────────

function test_method_improvement_search()
    println()
    println("=" ^ 72)
    println("Test 7: Method improvement search  [accuracy vs physics risk]")
    println("=" ^ 72)
    println("Literature-motivated candidates: checked finite-difference stencils,")
    println("velocity-channel Hanning filtering, dn/min_valid changes, and rbar2 use.")
    println("Baseline is previous production semantics: backward_current, dn=10, min_valid=50.")
    println()

    seeds = 301:308
    scenarios = [
        ("uniform noisy",
         seed -> make_uniform_cube(; noise_level = 0.12, rng = MersenneTwister(seed)),
         (i, j) -> TRUE_ANGLE, 8),
        ("weak gradient",
         seed -> make_uniform_cube(; gradient_amplitude = 0.05, noise_level = 0.05, rng = MersenneTwister(seed)),
         (i, j) -> TRUE_ANGLE, 8),
        ("sparse fill",
         seed -> make_sparse_cube(; fill_fraction = 0.35, noise_level = 0.05, rng = MersenneTwister(seed)),
         (i, j) -> TRUE_ANGLE, 12),
        ("varying mild",
         seed -> make_varying_angle_cube(; noise_level = 0.05, curvature = 8e-4, rng = MersenneTwister(seed)),
         (i, j) -> varying_true_angle(i, j, 128, 128; curvature = 8e-4), 12),
        ("varying sharp",
         seed -> make_varying_angle_cube(; noise_level = 0.05, curvature = 2.0e-3, rng = MersenneTwister(seed)),
         (i, j) -> varying_true_angle(i, j, 128, 128; curvature = 2.0e-3), 12),
    ]

    candidates = [
        ("central_checked",
         (; dn = 10, min_valid = 50, gradient_mode = :central_checked,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = false)),
        ("sobel_checked",
         (; dn = 10, min_valid = 50, gradient_mode = :sobel_checked,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = false)),
        ("backward_checked",
         (; dn = 10, min_valid = 50, gradient_mode = :backward_checked,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = false)),
        ("dn14_frac50",
         (; dn = 14, min_valid = 98, gradient_mode = :backward_current,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = false)),
        ("dn14_central",
         (; dn = 14, min_valid = 98, gradient_mode = :central_checked,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = false)),
        ("central_rbar_weight",
         (; dn = 10, min_valid = 50, gradient_mode = :central_checked,
          rbar2_weight = true, min_rbar2 = 0.0, spectral_hanning = false)),
        ("central_gate_0.50",
         (; dn = 10, min_valid = 50, gradient_mode = :central_checked,
          rbar2_weight = false, min_rbar2 = 0.50, spectral_hanning = false)),
        ("hanning3_baseline",
         (; dn = 10, min_valid = 50, gradient_mode = :backward_current,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = true)),
        ("hanning3_central",
         (; dn = 10, min_valid = 50, gradient_mode = :central_checked,
          rbar2_weight = false, min_rbar2 = 0.0, spectral_hanning = true)),
    ]

    candidate_wins = Dict(label => 0 for (label, _) in candidates)
    candidate_losses = Dict(label => 0 for (label, _) in candidates)
    candidate_rejects = Dict(label => String[] for (label, _) in candidates)
    rows = []

    @printf("%-16s  %-20s  %-9s  %-9s  %-7s  %-7s  %-7s  %s\n",
            "scenario", "candidate", "base_med", "cand_med", "base_n", "cand_n", "ratio", "verdict")
    println("-" ^ 100)

    for (scenario_label, cube_fn, angle_fn, stride) in scenarios
        base_errs = Float64[]
        cand_errs = Dict(label => Float64[] for (label, _) in candidates)
        for seed in seeds
            cube = cube_fn(seed)
            append!(base_errs, collect_sampled_errors_experimental(
                cube,
                angle_fn;
                stride = stride,
                dn = 10,
                min_valid = 50,
                gradient_mode = :backward_current,
            ))
            for (label, kwargs) in candidates
                append!(cand_errs[label], collect_sampled_errors_experimental(
                    cube,
                    angle_fn;
                    stride = stride,
                    dn = kwargs.dn,
                    min_valid = kwargs.min_valid,
                    gradient_mode = kwargs.gradient_mode,
                    rbar2_weight = kwargs.rbar2_weight,
                    min_rbar2 = kwargs.min_rbar2,
                    spectral_hanning = kwargs.spectral_hanning,
                ))
            end
        end

        base = error_summary(base_errs)
        for (label, _) in candidates
            cand = error_summary(cand_errs[label])
            coverage_ratio = base.n == 0 ? 0.0 : cand.n / base.n
            med_delta = cand.med - base.med
            p90_delta = cand.p90 - base.p90
            improved = isfinite(med_delta) && med_delta < -0.05
            severe_loss = coverage_ratio < 0.85
            sharp_blur = scenario_label == "varying sharp" &&
                         isfinite(med_delta) && med_delta > 0.10
            verdict = severe_loss ? "reject: coverage" :
                      (sharp_blur ? "reject: sharp-field bias" :
                       (improved ? "improves" : "no gain"))
            if verdict == "improves"
                candidate_wins[label] += 1
            elseif startswith(verdict, "reject")
                push!(candidate_rejects[label], "$scenario_label $verdict")
            else
                candidate_losses[label] += 1
            end
            @printf("%-16s  %-20s  %-9.3f  %-9.3f  %-7d  %-7d  %-7.2f  %s\n",
                    scenario_label, label, base.med, cand.med, base.n, cand.n,
                    coverage_ratio, verdict)
            push!(rows, (scenario_label, label, base, cand, coverage_ratio, verdict))
            if isfinite(p90_delta) && p90_delta > 1.0
                push!(candidate_rejects[label], "$scenario_label p90 worsened by $(round(p90_delta, digits=2)) deg")
            end
        end
    end

    println()
    println("Candidate aggregate:")
    @printf("%-20s  %-8s  %-8s  %s\n", "candidate", "improves", "no_gain", "reject notes")
    println("-" ^ 76)
    for (label, _) in candidates
        notes = candidate_rejects[label]
        note_text = isempty(notes) ? "" : join(notes, "; ")
        @printf("%-20s  %-8d  %-8d  %s\n",
                label, candidate_wins[label], candidate_losses[label], note_text)
    end

    return (rows = rows, wins = candidate_wins,
            losses = candidate_losses, rejects = candidate_rejects)
end

# ─── Summary ────────────────────────────────────────────────────────────────

function print_summary(sweeps, r6, r7)
    println()
    println("=" ^ 72)
    println("SUMMARY — axial_vonmises on synthetic cubes")
    println("=" ^ 72)
    println()

    @printf("%-20s  %-10s  %-10s\n", "Test", "worst_med°", "best_med°")
    println("-" ^ 44)
    for r in sweeps
        finite = filter(isfinite, r.med)
        isempty(finite) && continue
        @printf("%-20s  %-10.3f  %-10.3f\n", r.name, maximum(finite), minimum(finite))
    end
    println()
    println("Representativeness notes:")
    println("  Test 1 — fundamental SNR scaling; covers the current HI and CO benchmark noise range")
    println("  Test 2 — patchy sparse-emission stress test; low-fill cases may legitimately fail")
    println("  Test 3 — separates gradient magnitude from noise; tests weak-signal regime")
    println("  Test 4 — spatial structure: exact truth is the finite-difference gradient")
    println("           orientation of the curved scalar field")
    println("  Test 5 — bimodal windows have no unique truth at the boundary; boundary rows")
    println("           report mean-orientation behavior, not a universal correctness criterion")
    println("  Test 6 — accuracy candidates are experimental; see verdict rows before adopting")
    println("  Test 7 — method improvements must help noisy/sparse cases without")
    println("           blurring sharp field structure or losing too much coverage")
    println()
    if isempty(r6.adopted)
        println("Accuracy candidates: no production change recommended from this run.")
    else
        println("Accuracy candidates needing follow-up:")
        for row in r6.adopted
            println("  $(row[1]) / $(row[2]): $(round(row[3], digits=3))° → $(round(row[4], digits=3))°")
        end
    end
    println()
    println("Method-improvement aggregate:")
    for (label, wins) in sort(collect(r7.wins); by = x -> (-x[2], x[1]))
        reject_count = length(r7.rejects[label])
        println("  $label: improves $wins scenarios; reject notes $reject_count")
    end
end

function main()
    println("Synthetic validation of the axial_vonmises fitter")
    println("True angle: $(TRUE_ANGLE_DEG)°  gradient_amplitude=0.15  pipeline_noise=1.0")
    println()

    sweeps = [test_noise_sweep(), test_sparse_fill(), test_amplitude_sweep(),
              test_varying_angle(), test_boundary()]
    r6 = test_accuracy_candidates()
    r7 = test_method_improvement_search()
    print_summary(sweeps, r6, r7)
end

main()
