using DrWatson
quickactivate(@__DIR__, "VGT")

# Synthetic validation of opt-in gradient whitening (`gradient_whiten`,
# `gradient_noise_cov`). Synthetic data only.
#
# Noise: white Gaussian noise convolved with an elliptical Gaussian (FWHM 3 px,
# axis ratio 1.05, gradient anisotropy A ~ 0.05, long axis 45 deg from the
# truth). Signal on a baseline of 1000, so the intensity mask never acts and
# only the gradient S/N matters. Per-pixel Q/U is smoothed with ker = 4 and
# direct-fit support, as in the production pipeline. dn = 20, min_valid = 100.
#
# Cases (VGT_GW_CASES, comma separated; default "ramp"):
#   ramp        plane ramp, whitening off / on over the gradient S/N grid
#   ramp_mis    ramp, covariance mis-estimated by -20% / +20% in A
#   isotropic   isotropic correlated noise: whitened result must equal off
#   varying     noise anisotropy direction rotating 30 deg left to right: off /
#               one global C / per-pixel covariance maps (FWHM 3 x dn and 1.5 x dn)
#   structured  MW-like filament morphologies (validate_mw_like_synthetic.jl
#               generators), off / on / on with A mis-estimated by 20%
# Other environment: VGT_GW_SIZE (128), VGT_GW_NV (8), VGT_GW_SEEDS (4),
# VGT_GW_REF_NV (48) and VGT_GW_SPAN (30 deg) for the varying case.
#
# Run: LOCAL_OK=1 VGT_GW_CASES=ramp julia --project=. scripts/validate_gradient_whitening.jl

import VGT
using Printf
using Random
using Statistics
using VGT: VGT_QU_mv_2d_profiled, finalize_vgt_products, headless_angle_difference,
           measure_gradient_noise, gradient_anisotropy, gradient_noise_cov_map

const SIZE = parse(Int, get(ENV, "VGT_GW_SIZE", "128"))
const NV = parse(Int, get(ENV, "VGT_GW_NV", "8"))
const NSEEDS = parse(Int, get(ENV, "VGT_GW_SEEDS", "4"))
const DN = 20
const REF_NV = parse(Int, get(ENV, "VGT_GW_REF_NV", "48"))   # emission-free channels used to measure the covariance
const SPAN_DEG = parse(Float64, get(ENV, "VGT_GW_SPAN", "30")) # rotation of the noise axis across the map (varying case)
const MIN_VALID = 100
const KER = 4
const BASE = 1000.0
const FWHM = 3.0
const AXIS_RATIO = 1.05
const TRUE_DEG = 20.0                        # true gradient axis of the ramp (pixel frame)
const LONG_AXIS_DEG = TRUE_DEG + 45.0        # kernel long axis
const RAMP_LEVELS = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0)   # G / sigma_g
const STRUCT_LEVELS = (0.25, 0.5, 1.0, 2.0, 4.0)      # median signal |g| / sigma_g

function correlated_noise(nx, ny, nv, rng; fwhm=FWHM, ratio=AXIS_RATIO, long_axis=deg2rad(LONG_AXIS_DEG))
    σ_geo = fwhm / (2 * sqrt(2 * log(2)))
    σ_long = σ_geo * sqrt(ratio)
    σ_short = σ_geo / sqrt(ratio)
    r = ceil(Int, 4 * σ_long)
    c, s = cos(long_axis), sin(long_axis)
    kern = [exp(-0.5 * (((x * c + y * s) / σ_long)^2 + ((-x * s + y * c) / σ_short)^2))
            for x in -r:r, y in -r:r]
    kern ./= sqrt(sum(abs2, kern))
    out = Array{Float64}(undef, nx, ny, nv)
    for k in 1:nv
        white = randn(rng, nx + 2r, ny + 2r)
        for j in 1:ny, i in 1:nx
            acc = 0.0
            @inbounds for jj in -r:r, ii in -r:r
                acc += kern[ii + r + 1, jj + r + 1] * white[i + r + ii, j + r + jj]
            end
            out[i, j, k] = acc
        end
    end
    return out
end

function ramp_cube(nx, ny, nv, slope, θ)
    ic, jc = (nx + 1) / 2, (ny + 1) / 2
    cube = Array{Float64}(undef, nx, ny, nv)
    for k in 1:nv
        env = 1 + 0.4 * sin(2π * k / nv)
        for j in 1:ny, i in 1:nx
            cube[i, j, k] = BASE + slope * (cos(θ) * (i - ic) + sin(θ) * (j - jc)) * env
        end
    end
    return cube
end

# Covariance with the anisotropy A scaled by `factor` (trace and direction kept).
function scaled_anisotropy_cov(cov, factor)
    A, dir = gradient_anisotropy(cov)
    tr = cov[1] + cov[2]
    d = factor * A * tr * cos(2 * deg2rad(dir))
    sxy = factor * A * tr * sin(2 * deg2rad(dir)) / 2
    return [(tr + d) / 2, (tr - d) / 2, sxy]
end

# Production-like angle map psi (field frame, radians): kernel Q/U, then Q/U
# smoothing with direct-fit support.
function vgt_psi(cube; cov=nothing, covmap=nothing)
    local q, u
    redirect_stdout(devnull) do
        q, u, _ = VGT_QU_mv_2d_profiled(cube, DN, 1, size(cube, 3), 1.0;
            min_valid=MIN_VALID, gradient_whiten=cov !== nothing || covmap !== nothing,
            gradient_noise_cov=cov === nothing ? Float64[] : cov,
            gradient_noise_cov_map=covmap)
    end
    intensity = dropdims(sum(cube; dims=3); dims=3)
    return finalize_vgt_products(q, u, intensity, KER, :nan_before_smoothing, :direct_fit_support).psi
end

function interior_mask(nx, ny)
    h = div(DN, 2)
    m = falses(nx, ny)
    m[h:(nx - h), h:(ny - h)] .= true
    return m
end

# Doubled-angle mean offset from the truth (deg) and rms scatter about it.
function bias_scatter(est, truth)
    isempty(est) && return NaN, NaN
    z = mean(cis.(2 .* (est .- truth)))
    bias = rad2deg(0.5 * angle(z))
    scatter = sqrt(mean(abs2, rad2deg.(headless_angle_difference.(est, truth .+ deg2rad(bias)))))
    return bias, scatter
end

meanfinite(v) = (x = filter(isfinite, v); isempty(x) ? NaN : mean(x))
fmt(x; d=2) = isfinite(x) ? @sprintf("%+.*f", d, x) : "n/a"
fmtu(x; d=2) = isfinite(x) ? @sprintf("%.*f", d, x) : "n/a"

function reference_noise(ratio)
    ref = correlated_noise(2 * SIZE, 2 * SIZE, 4, MersenneTwister(1); ratio=ratio)
    return measure_gradient_noise(ref)
end

# ---- ramp -------------------------------------------------------------------
function ramp_case(levels, variants; ratio=AXIS_RATIO, title)
    m = reference_noise(ratio)
    θ_t = deg2rad(TRUE_DEG)
    truth_B = θ_t + π / 2
    D = rad2deg(headless_angle_difference(deg2rad(m.direction_deg), θ_t))
    σ_g = sqrt((m.sigma_x^2 + m.sigma_y^2) / 2)
    println("### $title")
    @printf("noise: sigma_gx=%.4f sigma_gy=%.4f  A=%.4f  preferred gradient direction=%.1f deg  D=%+.1f deg  cov=[%.5f, %.5f, %.5f]\n",
            m.sigma_x, m.sigma_y, m.anisotropy, m.direction_deg, D, m.cov...)
    println("| G/σ_g | whitening | bias (deg) | bias/D | scatter (deg) | coverage |")
    println("|---|---|---|---|---|---|")
    mask = interior_mask(SIZE, SIZE)
    for s in levels
        acc = Dict(name => (Float64[], Float64[], Float64[]) for (name, _) in variants)
        for seed in 1:NSEEDS
            noise = correlated_noise(SIZE, SIZE, NV, MersenneTwister(1000 + seed); ratio=ratio)
            cube = ramp_cube(SIZE, SIZE, NV, s * σ_g, θ_t) .+ noise
            for (name, cov) in variants
                psi = vgt_psi(cube; cov=cov)
                sel = mask .& isfinite.(psi)
                b, sc = bias_scatter(psi[sel], truth_B)
                push!(acc[name][1], b); push!(acc[name][2], sc); push!(acc[name][3], count(sel) / count(mask))
            end
        end
        for (name, _) in variants
            b, sc, cv = acc[name]
            ok = filter(isfinite, b)
            se = length(ok) > 1 ? std(ok) / sqrt(length(ok)) : NaN
            @printf("| %.2f | %s | %s ± %s | %s | %s | %s |\n", s, name, fmt(meanfinite(b)), fmtu(se),
                    fmt(meanfinite(b) / D), fmtu(meanfinite(sc)), fmtu(mean(cv); d=3))
        end
    end
    println()
    return m
end

# ---- spatially varying noise anisotropy --------------------------------------
# Noise whose kernel long axis rotates by `span_deg` from the left to the right
# edge: independent fields with long axes at equal steps, blended with
# sqrt(hat) weights (sum of squared weights = 1), so the covariance is the
# linear interpolation of the fields' covariances and the variance is kept.
function varying_noise(nx, ny, nv, rng; span_deg=SPAN_DEG, nfields=4)
    out = zeros(nx, ny, nv)
    spacing = (nx - 1) / (nfields - 1)
    for f in 1:nfields
        φ = deg2rad(LONG_AXIS_DEG - span_deg / 2 + span_deg * (f - 1) / (nfields - 1))
        field = correlated_noise(nx, ny, nv, rng; long_axis=φ)
        xf = 1 + spacing * (f - 1)
        for i in 1:nx
            w = sqrt(max(0.0, 1 - abs(i - xf) / spacing))
            w == 0 && continue
            @views out[i, :, :] .+= w .* field[i, :, :]
        end
    end
    return out
end

function thirds(nx)
    h = div(DN, 2)
    edges = round.(Int, range(h, nx - h; length=4))
    return [(edges[b] + (b > 1 ? 1 : 0)):edges[b + 1] for b in 1:3]
end

function varying_case()
    θ_t = deg2rad(TRUE_DEG); truth_B = θ_t + π / 2
    ref = varying_noise(SIZE, SIZE, REF_NV, MersenneTwister(1))
    glob = measure_gradient_noise(ref)
    maps = (("map FWHM 3 dn", gradient_noise_cov_map(ref; fwhm_pixels=3 * DN)),
            ("map FWHM 1.5 dn", gradient_noise_cov_map(ref; fwhm_pixels=1.5 * DN)))
    σ_g = sqrt((glob.sigma_x^2 + glob.sigma_y^2) / 2)
    bands = thirds(SIZE)
    println("### noise anisotropy direction rotating $(SPAN_DEG) deg left to right; $(REF_NV) emission-free channels for the covariance")
    println("local noise anisotropy A and preferred gradient direction in the left / middle / right third (measured on the reference cube):")
    for (b, r) in enumerate(bands)
        mb = measure_gradient_noise(ref[r, :, :])
        @printf("  third %d: A=%.4f direction=%.1f deg (D=%+.1f deg from the truth)\n", b, mb.anisotropy, mb.direction_deg,
                rad2deg(headless_angle_difference(deg2rad(mb.direction_deg), θ_t)))
    end
    @printf("global: A=%.4f direction=%.1f deg\n\n", glob.anisotropy, glob.direction_deg)

    # (1) Preferred direction left in pure-noise gradients: unit-weight resultant of 2*angle.
    println("residual preferred direction of whitened noise gradients (unit-weight resultant of the doubled angle, 0 = none) on an independent noise cube:")
    println("| correction | left third | middle third | right third |")
    println("|---|---|---|---|")
    test_noise = varying_noise(SIZE, SIZE, 16, MersenneTwister(77))
    gw = VGT._gradient_whitening(glob.cov)
    mm = Dict(name => VGT._whitening_maps(m) for (name, m) in maps)
    function band_resultant(r, how)
        z = 0.0im; n = 0
        for k in 1:16, j in 2:(SIZE - 1), i in r
            gx = (test_noise[i - 1, j, k] - test_noise[i + 1, j, k]) / 2
            gy = (test_noise[i, j - 1, k] - test_noise[i, j + 1, k]) / 2
            wx, wy = how(i, j, gx, gy)
            z += cis(2 * atan(wy, wx)); n += 1
        end
        return abs(z) / n
    end
    hows = Any[("off", (i, j, gx, gy) -> (gx, gy)),
               ("global C", (i, j, gx, gy) -> VGT._whiten_gradient(gw, gx, gy))]
    for (name, _) in maps
        push!(hows, (name, let m = mm[name]; (i, j, gx, gy) -> VGT._whiten_gradient_at(m, i, j, gx, gy) end))
    end
    for (name, how) in hows
        @printf("| %s | %s |\n", name, join([@sprintf("%.4f", band_resultant(r, how)) for r in bands], " | "))
    end
    println()

    # (2) Angle bias per third on the ramp.
    println("| G/σ_g | correction | bias left (deg) | bias middle | bias right | mean abs bias | scatter (deg) | coverage |")
    println("|---|---|---|---|---|---|---|---|")
    mask = interior_mask(SIZE, SIZE)
    variants = Any[("off", (nothing, nothing)), ("global C", (glob.cov, nothing))]
    for (name, m) in maps
        push!(variants, (name, (nothing, m)))
    end
    for s in (0.1, 0.25, 0.5, 1.0, 2.0)
        acc = Dict(name => (zeros(3, 0), Float64[], Float64[]) for (name, _) in variants)
        for seed in 1:NSEEDS
            noise = varying_noise(SIZE, SIZE, NV, MersenneTwister(1000 + seed))
            cube = ramp_cube(SIZE, SIZE, NV, s * σ_g, θ_t) .+ noise
            for (name, (cov, cmap)) in variants
                psi = vgt_psi(cube; cov=cov, covmap=cmap)
                bb = Float64[]
                for r in bands
                    sel = falses(SIZE, SIZE); sel[r, :] .= true
                    sel .&= mask .& isfinite.(psi)
                    push!(bb, bias_scatter(psi[sel], truth_B)[1])
                end
                sel = mask .& isfinite.(psi)
                a = acc[name]
                acc[name] = (hcat(a[1], bb), push!(a[2], bias_scatter(psi[sel], truth_B)[2]), push!(a[3], count(sel) / count(mask)))
            end
        end
        for (name, _) in variants
            bb, sc, cv = acc[name]
            per = [meanfinite(bb[b, :]) for b in 1:3]
            @printf("| %.2f | %s | %s | %s | %s | %s | %s | %s |\n", s, name, fmt(per[1]), fmt(per[2]), fmt(per[3]),
                    fmtu(mean(abs.(per))), fmtu(meanfinite(sc)), fmtu(mean(cv); d=3))
        end
    end
    println()
end

# ---- isotropic noise: whitened must equal off -------------------------------
function isotropic_case()
    m = reference_noise(1.0)
    println("### isotropic correlated noise (axis ratio 1): whitened vs off")
    @printf("measured A=%.4f (cov=[%.5f, %.5f, %.5f])\n", m.anisotropy, m.cov...)
    println("| G/σ_g | whitening | max abs difference in psi vs off (deg) | bias (deg) | scatter (deg) |")
    println("|---|---|---|---|---|")
    θ_t = deg2rad(TRUE_DEG); truth_B = θ_t + π / 2
    σ_g = sqrt((m.sigma_x^2 + m.sigma_y^2) / 2)
    mask = interior_mask(SIZE, SIZE)
    for s in (1.0, 4.0)
        noise = correlated_noise(SIZE, SIZE, NV, MersenneTwister(1001); ratio=1.0)
        cube = ramp_cube(SIZE, SIZE, NV, s * σ_g, θ_t) .+ noise
        off = vgt_psi(cube)
        for (name, cov) in (("off", nothing), ("on, C = identity", [1.0, 1.0, 0.0]), ("on, measured C", m.cov))
            psi = cov === nothing ? off : vgt_psi(cube; cov=cov)
            sel = mask .& isfinite.(psi) .& isfinite.(off)
            dmax = maximum(abs.(rad2deg.(headless_angle_difference.(psi[sel], off[sel]))))
            b, sc = bias_scatter(psi[sel], truth_B)
            @printf("| %.2f | %s | %.2e | %s | %s |\n", s, name, dmax, fmt(b), fmtu(sc))
        end
    end
    println()
end

# ---- structured emission ----------------------------------------------------
function structured_case()
    mw = Module(:MWLike)
    Base.include(mw, joinpath(@__DIR__, "validate_mw_like_synthetic.jl"))
    # The included functions are newer than this call's world age.
    return Base.invokelatest(structured_case, mw)
end

function structured_case(mw::Module)
    mw.NOISE_SCALE[] = 0.0                      # planted structure only
    m = reference_noise(AXIS_RATIO)
    σ_g = sqrt((m.sigma_x^2 + m.sigma_y^2) / 2)
    θ_n_B = deg2rad(m.direction_deg) + π / 2    # noise direction in the field frame
    println("### structured emission (MW-like generators), ker=$KER")
    @printf("noise: A=%.4f preferred direction (field frame)=%.1f deg\n", m.anisotropy, rad2deg(θ_n_B))
    println("Eligible pixels: signal above its 55th percentile and truth coherence >= 0.75 (as in validate_mw_like_synthetic.jl).")
    println("toward-noise bias = mean signed error with the sign set toward the noise axis; coverage = finite fraction of eligible pixels.")
    println("| scenario | S/N | whitening | median err (deg) | p90 err (deg) | toward-noise bias (deg) | coverage |")
    println("|---|---|---|---|---|---|---|")
    mask = interior_mask(SIZE, SIZE)
    variants = (("off", nothing), ("on", m.cov), ("on, A +20%", scaled_anisotropy_cov(m.cov, 1.2)),
                ("on, A -20%", scaled_anisotropy_cov(m.cov, 0.8)))
    for (label, builder) in (("many_filaments", mw.build_many_filaments), ("curved_fibers", mw.build_curved_fibers))
        sig, truth, _ = builder(; nx=SIZE, ny=SIZE, nv=NV)
        tang, rbar = mw.truth_angle_and_rbar(truth.q, truth.u, truth.w)
        cut = mw.qtile([v for v in vec(truth.w) if v > 0], 0.55)
        elig = mask .& (truth.w .> cut) .& isfinite.(tang) .& (rbar .>= 0.75)
        # Median central-difference gradient magnitude of the planted signal on eligible pixels.
        mags = Float64[]
        for k in 1:NV, j in 2:(SIZE - 1), i in 2:(SIZE - 1)
            elig[i, j] || continue
            push!(mags, hypot((sig[i - 1, j, k] - sig[i + 1, j, k]) / 2, (sig[i, j - 1, k] - sig[i, j + 1, k]) / 2))
        end
        med_g = median(mags)
        pull = sign.(headless_angle_difference.(θ_n_B, tang))
        for s in STRUCT_LEVELS
            α = s * σ_g / med_g
            acc = Dict(name => (Float64[], Float64[], Float64[], Float64[]) for (name, _) in variants)
            for seed in 1:NSEEDS
                noise = correlated_noise(SIZE, SIZE, NV, MersenneTwister(2000 + seed))
                cube = BASE .+ α .* sig .+ noise
                for (name, cov) in variants
                    psi = vgt_psi(cube; cov=cov)
                    err = rad2deg.(headless_angle_difference.(psi, tang))
                    good = elig .& isfinite.(psi)
                    a = abs.(err[good])
                    push!(acc[name][1], median(a)); push!(acc[name][2], quantile(a, 0.9))
                    push!(acc[name][3], mean(err[good] .* pull[good])); push!(acc[name][4], count(good) / count(elig))
                end
            end
            for (name, _) in variants
                a = acc[name]
                @printf("| %s | %.2f | %s | %s | %s | %s | %s |\n", label, s, name, fmtu(mean(a[1])),
                        fmtu(mean(a[2])), fmt(mean(a[3])), fmtu(mean(a[4]); d=3))
            end
        end
    end
    println()
end

function main()
    cases = split(get(ENV, "VGT_GW_CASES", "ramp"), ",")
    println("Gradient whitening validation: $(SIZE)x$(SIZE)x$NV, dn=$DN, min_valid=$MIN_VALID, ker=$KER, $NSEEDS seeds")
    println()
    if "ramp" in cases
        ramp_case(RAMP_LEVELS, (("off", nothing), ("on", reference_noise(AXIS_RATIO).cov)); title="ramp, anisotropic noise")
    end
    if "ramp_mis" in cases
        cov = reference_noise(AXIS_RATIO).cov
        ramp_case((0.25, 0.5, 1.0, 2.0),
                  (("on, A -20%", scaled_anisotropy_cov(cov, 0.8)), ("on, A +20%", scaled_anisotropy_cov(cov, 1.2)));
                  title="ramp, covariance mis-estimated in A")
    end
    "isotropic" in cases && isotropic_case()
    "varying" in cases && varying_case()
    "structured" in cases && structured_case()
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
