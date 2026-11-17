using DrWatson
quickactivate(@__DIR__, "VGT")

using CairoMakie
using FFTW
using FITSIO
using Printf
using LinearAlgebra: pinv
using Random
using Statistics
using VGT

# Grid-locked angle bias of VGT: the mean VGT angle minus the true field angle as a function of the true angle, for
# the gradient stencils :central_checked (default), :sobel_checked and :isotropic_checked.
#
# Why: foreground templates fitted to dust polarization need mean directions accurate to about a degree.
# The VGT kernel is symmetric under the square grid's mirrors (test/angle_calibration.jl), so it has no uniform
# offset; a stencil can only bias angles in a grid-locked pattern b(t) = A4 sin(4t) + A8 sin(8t) + ... (t = field
# angle counter-clockwise from +axis 1). Plane waves predict A4 = +k^2/24 rad for central differences (towards the
# diagonals), -k^2/48 for Sobel, O(k^4) for the [1, 4, 1]/6 cross-smoothed stencil; the size on real maps depends on
# the wavenumbers that carry the gradients (beam, spectrum, noise, window dispersion), hence this measurement.
#
# Declared before running (2026-11-17):
# - Material bias: |A4| > 1 deg at dn 20, ker 4, 3 px beam (the window settings of the foreground templates on
#   GASKAP 40" and GALFA 80" pixels; both have a 3 px beam, so in pixel units the two differ only in S/N).
#   If |A4| <= 1 deg for :central_checked at every S/N, the stencil is not a material bias and no correction is made.
#   If |A4| > 1 deg and :isotropic_checked brings it below 0.3 deg without losing alignment (AM within 0.01), the
#   new mode is the correction, used opt-in.
# - Checks: the cos(4t) term and the constant must be zero within errors (symmetry); a nonzero constant means a bug
#   or a realisation offset (CATS), not a calibratable bias.
#
# Part A, Gaussian random fields (exact truth, no interpolation): every channel is an independent field generated in
# Fourier space with power k^-beta, anisotropy exp(-zeta cos^2(phi_k - t)) (wave vectors along the field suppressed,
# so structures run along t), a Gaussian beam of FWHM 3 px, a uniform level 60 sigma_n, and beam-correlated noise of
# unit rms. S/N r = rms of the signal after removing its 20 px boxcar mean, in noise units (r = Inf: no noise).
# Same random phases at every angle, so the angle scan of one seed is one realisation viewed at many angles.
# Main grid: beta 8/3 (thin velocity channel of Kolmogorov turbulence, Lazarian & Pogosyan 2000), zeta 3, r in
# {Inf, 4, 2, 1}. Side runs at r = 2: beta 11/3, zeta 1 and 8, dn 10 and 40 (min_valid kept at 25 % of dn^2).
# Bias per run: half the angle of sum over interior pixels of exp(2i (psi - t)), all pixels weighted equally.
# A4, B4 (cos 4t), A8, c0 fitted per seed by least squares over t = 0, 5, ..., 85 deg (the 90 deg period follows
# from the symmetry); errors are the standard error over seeds.
#
# Part B, CATS MHD (Burkhart et al. 2020), realism check: optically thin PPV from scripts/synthetic_ppv_from_mhd.jl
# (los 3, as the calibration runs), beam 3 px, then rotated by t: the periodic, beam-smoothed cube is upsampled 4x by
# zero-padding its FFT and sampled with Keys cubic interpolation at the rotated positions (the beam leaves almost no
# power near Nyquist, so this is close to exact; checked by a 0 deg identity and a +45/-45 round trip). Noise
# (beam-correlated, a new realisation at each angle) is added after the rotation. The true angle is the rotated beam-smoothed Stokes field of the
# density-weighted plane-of-sky B. Bias as in A against the true map, over the disk of radius 96 px about the
# rotation centre (the same gas at every angle); per-angle error from a 32 px block bootstrap.
# The constant c0 here is this snapshot's own VGT-truth offset, which rotates with the sky; only A4, A8 are grid bias.
#
#   julia --project=. -t 32 scripts/validate_angle_bias.jl --out validation_output/angle_bias \
#       [--part grf|cats|both] [--seeds 8] [--cats-dir data/external/cats] [--cats-runs b1p1,b.1p1] [--quick 1]
# --configs TAG,...: only these GRF configurations (tags main, beta11/3, zeta1, zeta8, dn10, dn40).
# --seed-offset S: GRF seeds S+1..S+seeds, to split the seeds over jobs (summarize_angle_bias.py takes several CSVs).
# --quick 1: smoke test (angles 0 and 20, a 128 px GRF, one configuration, CATS with 40 channels).
#
# Writes <out>/grf_runs.csv, <out>/grf_fit.csv, <out>/cats_runs.csv, <out>/cats_fit.csv and PNG plots.

const MODES = (:central_checked, :sobel_checked, :isotropic_checked)
ANGLES = 0:5:85            # replaced by --quick
QUICK = false
const BEAM = 3.0
const KER = 4

function parse_args(args)
    o = Dict("out" => "validation_output/angle_bias", "part" => "both", "seeds" => "8",
             "cats-dir" => "data/external/cats", "cats-runs" => "b1p1,b.1p1", "cats-nv" => "218", "cats-sth" => "1.0")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("bad argument $(args[i])"))
        o[args[i][3:end]] = args[i + 1]
        i += 2
    end
    return o
end

beam_tf(kx, ky, fwhm) = exp(-0.5 * (fwhm / sqrt(8log(2)))^2 * (kx^2 + ky^2))
kgrid(n) = 2pi .* fftfreq(n)

# 20 px boxcar high-pass transfer (what the S/N r refers to).
function highpass_tf(kx, ky; L=20)
    d(k) = k == 0 ? 1.0 : sin(L * k / 2) / (L * sin(k / 2))
    return 1 - d(kx) * d(ky)
end

function beam_noise(rng, n, fwhm)
    k = kgrid(n)
    w = real.(ifft(fft(randn(rng, n, n)) .* [beam_tf(a, b, fwhm) for a in k, b in k]))
    return w ./ std(w)
end

# One channel of the anisotropic field at angle t (rad); `phase` is the FFT of white noise (fixed per seed/channel).
function grf_channel(phase, t, beta, zeta, r)
    n = size(phase, 1)
    k = kgrid(n)
    amp = [begin
               kk = hypot(a, b)
               kk == 0 ? 0.0 : kk^(-beta / 2) * exp(-0.5 * zeta * cos(atan(b, a) - t)^2) * beam_tf(a, b, BEAM)
           end for a in k, b in k]
    f = fft(real.(ifft(phase .* amp)))
    hp = real.(ifft(f .* [highpass_tf(a, b) for a in k, b in k]))
    s = real.(ifft(f))
    return s .* ((isinf(r) ? 2.0 : r) / std(hp))
end

axial(a) = mod(a + pi / 2, pi) - pi / 2

function run_vgt(cube, dn, noise, mode)
    nv = size(cube, 3)
    q, u, prof = VGT.VGT_QU_mv_2d_profiled(cube, dn, 1, nv, noise; min_valid=dn^2 ÷ 4, gradient_mode=mode)
    intensity = dropdims(sum(cube; dims=3); dims=3)
    psi = VGT.finalize_vgt_products(q, u, intensity, KER).psi
    return psi, prof[:fit_rbar2]
end

# Bias (deg) and alignment AM of psi against truth (scalar t or a map) over interior pixels: the square at least
# `edge` from the border, or with disk=true the disk of radius n/2 - edge about the rotation centre (CATS: the same
# gas is inside the disk at every rotation angle, so the angle scan sees one fixed sample).
in_region(i, j, n, edge, disk) = disk ? hypot(i - (n + 1) / 2, j - (n + 1) / 2) <= n / 2 - edge :
                                 (edge < i <= n - edge && edge < j <= n - edge)

function bias_stats(psi, truth, edge; disk=false)
    n = size(psi, 1)
    z = 0.0im
    am = 0.0
    m = 0
    for j in 1:n, i in 1:n
        in_region(i, j, n, edge, disk) || continue
        p = psi[i, j]
        tr = truth isa Number ? truth : truth[i, j]
        (isfinite(p) && isfinite(tr)) || continue
        z += cis(2 * (p - tr))
        am += cos(2 * (p - tr))
        m += 1
    end
    return rad2deg(0.5 * angle(z)), am / max(m, 1), m
end

# Least squares of b(t) on [1, sin 4t, cos 4t, sin 8t]; returns (c0, A4, B4, A8).
function harmonics(ts_deg, b; w=ones(length(b)))
    t = deg2rad.(collect(ts_deg))
    X = hcat(ones(length(t)), sin.(4t), cos.(4t), sin.(8t))
    W = sqrt.(w)
    return (X .* W) \ (b .* W)
end

function grf_part(o)
    out = o["out"]
    seeds = parse(Int, o["seeds"])
    seed0 = parse(Int, get(o, "seed-offset", "0"))   # run seeds seed0+1..seed0+seeds (split over jobs)
    n, nv = QUICK ? 128 : 512, 8
    configs = Any[]
    for r in (Inf, 4.0, 2.0, 1.0)
        push!(configs, (; tag="main", beta=8 / 3, zeta=3.0, r, dn=20))
    end
    push!(configs, (; tag="beta11/3", beta=11 / 3, zeta=3.0, r=2.0, dn=20))
    push!(configs, (; tag="zeta1", beta=8 / 3, zeta=1.0, r=2.0, dn=20))
    push!(configs, (; tag="zeta8", beta=8 / 3, zeta=8.0, r=2.0, dn=20))
    push!(configs, (; tag="dn10", beta=8 / 3, zeta=3.0, r=2.0, dn=10))
    push!(configs, (; tag="dn40", beta=8 / 3, zeta=3.0, r=2.0, dn=40))
    QUICK && (configs = configs[3:3])
    if haskey(o, "configs")        # e.g. --configs beta11/3,zeta8: only these tags
        keep = split(o["configs"], ",")
        configs = [c for c in configs if c.tag in keep]
    end
    rows = open(joinpath(out, "grf_runs.csv"), "w")
    println(rows, "tag,beta,zeta,r,dn,mode,seed,angle_deg,bias_deg,am,rbar2,npix")
    results = Dict{Any, Array{Float64}}()   # (config, mode) => [seed, angle] bias
    ams = Dict{Any, Array{Float64}}()
    rbs = Dict{Any, Array{Float64}}()
    for c in configs, mode in MODES
        results[(c, mode)] = fill(NaN, seeds, length(ANGLES))
        ams[(c, mode)] = fill(NaN, seeds, length(ANGLES))
        rbs[(c, mode)] = fill(NaN, seeds, length(ANGLES))
    end
    for seed in 1:seeds
        rng = MersenneTwister(1000 + seed0 + seed)
        phases = [fft(randn(rng, n, n)) for _ in 1:nv]
        noises = [beam_noise(rng, n, BEAM) for _ in 1:nv]
        for c in configs, (ia, a) in enumerate(ANGLES)
            t = deg2rad(a)
            cube = Array{Float64}(undef, n, n, nv)
            for k in 1:nv
                cube[:, :, k] .= 60 .+ grf_channel(phases[k], t, c.beta, c.zeta, c.r) .+ (isinf(c.r) ? 0.0 : 1.0) .* noises[k]
            end
            for mode in MODES
                psi, rb = run_vgt(cube, c.dn, 3.0, mode)
                edge = c.dn + 3KER
                b, am, m = bias_stats(psi, t, edge)
                rbar = mean(filter(isfinite, rb[(edge + 1):(n - edge), (edge + 1):(n - edge)]))
                results[(c, mode)][seed, ia] = b
                ams[(c, mode)][seed, ia] = am
                rbs[(c, mode)][seed, ia] = rbar
                @printf(rows, "%s,%.4f,%.1f,%.1f,%d,%s,%d,%d,%.5f,%.5f,%.5f,%d\n", c.tag, c.beta, c.zeta, c.r, c.dn,
                        mode, seed0 + seed, a, b, am, rbar, m)
            end
        end
        flush(rows)
        @info "GRF seed done" seed
    end
    close(rows)
    fitf = open(joinpath(out, "grf_fit.csv"), "w")
    println(fitf, "tag,beta,zeta,r,dn,mode,c0,c0_err,A4,A4_err,B4,B4_err,A8,A8_err,AM,rbar2")
    for c in configs, mode in MODES
        B = results[(c, mode)]
        H = reduce(hcat, [harmonics(ANGLES, B[s, :]) for s in 1:seeds])   # 4 x seeds
        m, e = vec(mean(H; dims=2)), vec(std(H; dims=2)) ./ sqrt(seeds)
        @printf(fitf, "%s,%.4f,%.1f,%.1f,%d,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", c.tag, c.beta,
                c.zeta, c.r, c.dn, mode, m[1], e[1], m[2], e[2], m[3], e[3], m[4], e[4], mean(ams[(c, mode)]),
                mean(rbs[(c, mode)]))
        @printf("GRF %-9s r=%-4s dn=%2d %-18s A4 %+6.3f ± %.3f  B4 %+6.3f ± %.3f  A8 %+6.3f ± %.3f  c0 %+6.3f ± %.3f  AM %.3f  rbar2 %.3f\n",
                c.tag, string(c.r), c.dn, mode, m[2], e[2], m[3], e[3], m[4], e[4], m[1], e[1],
                mean(ams[(c, mode)]), mean(rbs[(c, mode)]))
    end
    close(fitf)
    plot_bias(joinpath(out, "grf_bias.png"), "GRF main (beta 8/3, zeta 3, dn 20)",
              [("S/N $(c.r)", mode, vec(mean(results[(c, mode)]; dims=1)),
                vec(std(results[(c, mode)]; dims=1)) ./ sqrt(seeds)) for c in configs if c.tag == "main" for mode in MODES])
end

function plot_bias(path, title, curves)
    fig = Figure(; size=(1100, 420))
    rs = unique(first.(curves))
    for (ir, r) in enumerate(rs)
        ax = Axis(fig[1, ir]; title="$title, $r", xlabel="true angle (deg)", ylabel="VGT - true (deg)")
        for (r2, mode, b, e) in curves
            r2 == r || continue
            errorbars!(ax, collect(ANGLES), b, e)
            scatterlines!(ax, collect(ANGLES), b; label=String(mode))
        end
        hlines!(ax, [0.0]; color=:gray)
        ir == 1 && axislegend(ax; position=:lb)
    end
    save(path, fig)
end

# ---------------------------------------------------------------- CATS
const UP = 4

keys_w(x) = (x = abs(x); x < 1 ? 1.5x^3 - 2.5x^2 + 1 : x < 2 ? -0.5x^3 + 2.5x^2 - 4x + 2 : 0.0)

# 4x trigonometric upsampling of a periodic n x n image (the Nyquist row/column is dropped; the 3 px beam leaves
# exp(-8) of the power there).
function upsample(img)
    n = size(img, 1)
    F = fft(img)
    m = UP * n
    G = zeros(ComplexF64, m, m)
    h = n ÷ 2
    for j in 1:n, i in 1:n
        ki, kj = i <= h ? i - 1 : i - 1 - n, j <= h ? j - 1 : j - 1 - n
        (abs(ki) == h || abs(kj) == h) && continue
        G[mod(ki, m) + 1, mod(kj, m) + 1] = F[i, j]
    end
    return real.(ifft(G)) .* UP^2
end

# Sample the upsampled periodic image at map positions rotated by t about the map centre: out(x) = in(R(-t)(x - c) + c).
function rotate_sample(fine, n, t)
    m = size(fine, 1)
    c = (n + 1) / 2
    ct, st = cos(t), sin(t)
    out = Array{Float64}(undef, n, n)
    Threads.@threads for j in 1:n
        for i in 1:n
            x = ct * (i - c) + st * (j - c) + c
            y = -st * (i - c) + ct * (j - c) + c
            fx, fy = (x - 1) * UP, (y - 1) * UP          # 0-based fine coordinates
            ix, iy = floor(Int, fx), floor(Int, fy)
            acc = 0.0
            for dy in -1:2
                wy = keys_w(fy - (iy + dy))
                wy == 0 && continue
                for dx in -1:2
                    wx = keys_w(fx - (ix + dx))
                    wx == 0 && continue
                    acc += wx * wy * fine[mod(ix + dx, m) + 1, mod(iy + dy, m) + 1]
                end
            end
            out[i, j] = acc
        end
    end
    return out
end

function beam2d(a, fwhm)
    k = kgrid(size(a, 1))
    return real.(ifft(fft(a) .* [beam_tf(x, y, fwhm) for x in k, y in k]))
end

rd(dir, name) = Float64.(FITS(f -> read(f[1]), joinpath(dir, "$(name)_t800.fits")))

function block_bootstrap_bias(psi, truth, edge; B=32, nboot=300, rng=MersenneTwister(5))
    n = size(psi, 1)
    blocks = Dict{Tuple{Int, Int}, ComplexF64}()
    for j in 1:n, i in 1:n
        in_region(i, j, n, edge, true) || continue
        (isfinite(psi[i, j]) && isfinite(truth[i, j])) || continue
        key = ((i - 1) ÷ B, (j - 1) ÷ B)
        blocks[key] = get(blocks, key, 0.0im) + cis(2 * (psi[i, j] - truth[i, j]))
    end
    z = collect(values(blocks))
    draws = [rad2deg(0.5 * angle(sum(z[rand(rng, 1:length(z), length(z))]))) for _ in 1:nboot]
    return std(draws), length(z)
end

function cats_part(o)
    out = o["out"]
    root = @__DIR__
    nvs, sth = o["cats-nv"], o["cats-sth"]
    rows = open(joinpath(out, "cats_runs.csv"), "w")
    println(rows, "run,snr,mode,angle_deg,bias_deg,bias_err,am,nblocks")
    fitf = open(joinpath(out, "cats_fit.csv"), "w")
    println(fitf, "run,snr,mode,c0,c0_err,A4,A4_err,B4,B4_err,A8,A8_err,AM")
    curves = Any[]
    for run in split(o["cats-runs"], ",")
        dir = joinpath(o["cats-dir"], run)
        ppvdir = joinpath(out, "cats_$(run)_beam")
        if !isfile(joinpath(ppvdir, "ppv.fits"))
            Base.run(`$(Base.julia_cmd()) --project=$(dirname(root)) -t $(Threads.nthreads()) $(joinpath(root, "synthetic_ppv_from_mhd.jl")) --dir $dir --snap t800 --los 3 --sigma-th $sth --nv $nvs --beam-px $(BEAM) --snr 0 --out $ppvdir`)
        end
        ppv = Float64.(FITS(f -> read(f[1]), joinpath(ppvdir, "ppv.fits")))
        n, _, nv = size(ppv)
        # true Stokes of the density-weighted plane-of-sky field (as synthetic_ppv_from_mhd.jl), beam-smoothed
        rho, b1, b2 = rd(dir, "dens"), rd(dir, "magx"), rd(dir, "magy")
        q = beam2d(dropdims(sum(rho .* (b1 .^ 2 .- b2 .^ 2); dims=3); dims=3), BEAM)
        u = beam2d(dropdims(sum(rho .* 2 .* b1 .* b2; dims=3); dims=3), BEAM)
        rho = b1 = b2 = nothing
        sigma23 = median(maximum(ppv; dims=3))
        fine = [upsample(ppv[:, :, k]) for k in 1:nv]
        fq, fu = upsample(q), upsample(u)
        # interpolation checks
        id0 = maximum(abs.(rotate_sample(fine[nv ÷ 2], n, 0.0) .- ppv[:, :, nv ÷ 2])) / maximum(abs.(ppv[:, :, nv ÷ 2]))
        r45 = rotate_sample(upsample(rotate_sample(fine[nv ÷ 2], n, deg2rad(45))), n, deg2rad(-45))
        c = n ÷ 4:3n ÷ 4
        rt = maximum(abs.(r45[c, c] .- ppv[c, c, nv ÷ 2])) / maximum(abs.(ppv[:, :, nv ÷ 2]))
        @printf("CATS %s: interpolation identity at 0 deg %.2e, +45/-45 round trip (centre) %.2e of the peak\n", run, id0, rt)
        for snr in (QUICK ? (46.0,) : (Inf, 46.0, 8.0))
            sig_n = sigma23 / (isinf(snr) ? 46.0 : snr)   # the noiseless run keeps the S/N 46 mask threshold
            B = Dict(mode => fill(NaN, length(ANGLES)) for mode in MODES)
            E = Dict(mode => fill(NaN, length(ANGLES)) for mode in MODES)
            A = Dict(mode => fill(NaN, length(ANGLES)) for mode in MODES)
            for (ia, a) in enumerate(ANGLES)
                t = deg2rad(a)
                cube = Array{Float64}(undef, n, n, nv)
                rng = MersenneTwister(77 + ia)   # a new noise realisation per angle (a fixed one rotates with the grid,
                                                 # not the sky, and makes a sawtooth b(t) where the noise dominates)
                for k in 1:nv
                    cube[:, :, k] .= rotate_sample(fine[k], n, t)
                    isinf(snr) || (cube[:, :, k] .+= sig_n .* beam_noise(rng, n, BEAM))
                end
                truth = 0.5 .* atan.(rotate_sample(fu, n, t), rotate_sample(fq, n, t)) .+ t
                for mode in MODES
                    psi, _ = run_vgt(cube, 20, 3sig_n, mode)
                    edge = 20 + 3KER
                    b, am, _ = bias_stats(psi, truth, edge; disk=true)
                    e, nb = block_bootstrap_bias(psi, truth, edge)
                    B[mode][ia], E[mode][ia], A[mode][ia] = b, e, am
                    @printf(rows, "%s,%.0f,%s,%d,%.5f,%.5f,%.5f,%d\n", run, snr, mode, a, b, e, am, nb)
                end
                flush(rows)
            end
            for mode in MODES
                h = harmonics(ANGLES, B[mode]; w=1 ./ E[mode] .^ 2)
                # error of the harmonics from the per-angle bootstrap errors (angles treated as independent)
                t = deg2rad.(collect(ANGLES))
                X = hcat(ones(length(t)), sin.(4t), cos.(4t), sin.(8t))
                C = pinv(X' * (X ./ E[mode] .^ 2))
                he = sqrt.(max.(0.0, [C[i, i] for i in 1:4]))
                @printf(fitf, "%s,%.0f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", run, snr, mode, h[1], he[1],
                        h[2], he[2], h[3], he[3], h[4], he[4], mean(A[mode]))
                @printf("CATS %-6s S/N %-4s %-18s A4 %+6.3f ± %.3f  B4 %+6.3f ± %.3f  A8 %+6.3f ± %.3f  c0 %+6.3f  AM %.3f\n",
                        run, string(snr), mode, h[2], he[2], h[3], he[3], h[4], he[4], h[1], mean(A[mode]))
                push!(curves, ("$run $(snr)", mode, B[mode] .- h[1], E[mode]))
            end
            flush(fitf)
        end
    end
    close(rows)
    close(fitf)
    plot_bias(joinpath(out, "cats_bias.png"), "CATS (minus own c0)", curves)
end

function main(args)
    o = parse_args(args)
    if get(o, "quick", "0") == "1"
        global ANGLES = (0, 20)
        global QUICK = true
        o["cats-nv"] = "40"
        o["cats-runs"] = first(split(o["cats-runs"], ","))
    end
    mkpath(o["out"])
    o["part"] in ("grf", "both") && grf_part(o)
    o["part"] in ("cats", "both") && cats_part(o)
end

main(ARGS)
