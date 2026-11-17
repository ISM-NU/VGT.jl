using DrWatson
quickactivate(@__DIR__, "VGT")

using FFTW
using FITSIO
using Printf
using Random
using Statistics
using VGT

# Faint-gas VGT on MHD cubes with known field: does weighting gradients by |g|^2 (V1,
# gradient_weight = :amplitude2) or a Gaussian pre-filter of each channel before the
# gradients (V2, the Lazarian group's noise recipe; mask and Stokes weights stay on the
# unfiltered channel) bring a foreground template's direction closer to the true
# plane-of-sky field where the gradients are near the noise? Validation in
# docs/validation-faint-gas.md.
#
# Declared before any result (2026-11-17):
#   Cubes: CATS (Burkhart et al. 2020) b.1p.32, b.1p1, b1p.32, b1p1 at t800, lines of
#     sight 2 and 3 (the mean field is along axis 1, so it lies in the sky plane for
#     both): 8 cubes. Optically thin PPV, nv 48 over +-4 std(v_los), beam 3 px FWHM,
#     beam-correlated noise; peak S/N (median over the map of the spectrum peak / sigma)
#     1, 2, 4, 8, 23; noise seeds 1..NSEED, the same noise for every variant (paired).
#   VGT as in foreground runs on a 2' H I cube with 40" pixels: 4-channel windows sliding by one channel, all of
#     whose channels have |v| <= 2.5 std(v_los); dn 20, min_valid 100, central_checked,
#     scalar mask at 3 sigma, ker 4.
#   Template as the two-layer foreground template of the experiments: per window t, ZS_t = column-
#     weighted Stokes of psi_t at FWHM 15 px (10' on a 40" grid) and U_t = unit
#     Stokes of psi_t at 180 px (120'). Layers: "all" (every window) and "wing" (windows
#     whose centre is below -1 std(v_los): a fainter layer on one side of the line, as
#     a faint foreground layer). Direction at 120' = arg sum_t |ZS_t| U_t; at 10' =
#     arg sum_t ZS_t.
#   Truth: velocity-resolved density-weighted Stokes of the sky-plane field, summed over
#     the layer's channels, beam 3 px, smoothed to the same scale.
#   Score: AM = mean cos 2(psi - psi_true) over pixels finite for every variant of that
#     cube / S/N / seed; written per 32 px block for a paired block bootstrap
#     (scripts/figures/faint_gas_summary.py).
#   Variants: plain; V1 amplitude2; V2 pre-filter FWHM 4 px (chosen before running: the
#     group uses sigma 2-4 px; 4 px keeps the filtered beam, 5 px, well below the 15 px
#     scale and the 20 px block); V12 both (sensitivity only, not part of the pass rule).
#   Pass rule for a variant: (1) wing layer, at the S/N levels whose median rho_g brackets
#     the faint layer of the target cube (the range is passed to the summary script) and
#     pooled over S/N 2, 4, 8: paired dAM at 10' > 0 at >= 2 sigma and dAM at 120' not
#     negative beyond 1 sigma; (2) both layers at S/N 23: dAM >= -0.005 at both scales;
#     (3) the sign of (1) in at least 6 of the 8 cubes.
#   Amended 2026-11-17 after a smoke run and an S/N scan (--scan 1), before any variant was scored: at peak S/N
#     1-4 the 3 sigma mask leaves (almost) no pixel and rho_g is 0, so those levels cannot test anything. The scan put
#     the wing layer's median rho_g in the target faint range at peak S/N 12-23 (0.6-1.9 over the four runs, los 3)
#     and the target's bright layer at 23-32. New grid: 8, 12, 16, 23, 46. Rule (1) pools 8, 12, 16 instead of
#     2, 4, 8; rule (2) holds at 23 and 46; rule (3) at 8, 12, 16 and the matched levels.
#   rho_g (computed the same way on the target cube): per 20 px block, rho_g^2 = <|g|^2> / <2 sigma_g^2> - 1
#     over masked-in pixels of the window's channels, sigma_g from the injected noise.
#
#   julia -t 16 --project=. scripts/validate_faint_gas.jl --cats DIR --out DIR [--nseed 3]
#       [--runs b.1p.32,b.1p1,b1p.32,b1p1] [--los 2,3] [--snr 8,12,16,23,46] [--maps 1] [--scan 1]
#   --scan 1 prints rho_g and mask coverage per layer and S/N only (no VGT).

include(joinpath(@__DIR__, "synthetic_ppv_from_mhd.jl"))   # project_ppv, beam!, rd

const BEAM = 3.0
const PREFILTER = 4.0
const DN = 20
const KER = 4
const S10 = 15.0
const S120 = 180.0
const BLOCK = 32
const VARIANTS = ("plain", "amp2", "pre4", "amp2pre4")

function opts_from(args)
    o = Dict("runs" => "b.1p.32,b.1p1,b1p.32,b1p1", "los" => "2,3", "snr" => "8,12,16,23,46", "nseed" => "3",
             "nv" => "48", "maps" => "1", "scan" => "0")
    i = 1
    while i <= length(args)
        o[args[i][3:end]] = args[i + 1]
        i += 2
    end
    haskey(o, "cats") && haskey(o, "out") || throw(ArgumentError("--cats and --out are required"))
    return o
end

# Periodic Gaussian smoothing of a 2D field (the box is periodic), FWHM in px.
function psmooth(a::AbstractMatrix, fwhm)
    s = fwhm / sqrt(8log(2))
    kx, ky = fftfreq(size(a, 1)), fftfreq(size(a, 2))
    g = [exp(-2pi^2 * s^2 * (x^2 + y^2)) for x in kx, y in ky]
    r = ifft(fft(a) .* g)
    return eltype(a) <: Complex ? r : real.(r)
end

unitc(z) = abs(z) > 0 ? z / abs(z) : complex(NaN, NaN)

# Stokes average of exp(2i psi), weighted by w, normalised by the covered fraction; NaN
# under half coverage (as planck_fg_v2_lib.column_stokes / smooth_unit).
function stokes_at(psi, w, fwhm)
    ok = isfinite.(psi)
    e = ifelse.(ok, w .* cis.(2 .* ifelse.(ok, psi, 0.0)), 0.0im)
    cov = psmooth(Float64.(ok), fwhm)
    z = psmooth(e, fwhm) ./ max.(cov, 1e-9)
    z[cov .< 0.5] .= complex(NaN, NaN)
    return z
end

function window_psi(win_raw, win_grad, sigma, variant)
    amp = variant in ("amp2", "amp2pre4")
    q, u = VGT.VGT_QU_mv_2d(win_grad, DN, 1, size(win_raw, 3), 3sigma; min_valid=100,
                            gradient_mode=:central_checked, weight_cube=win_raw,
                            gradient_weight=amp ? :amplitude2 : :unit)
    I = dropdims(sum(win_raw; dims=3); dims=3)
    return VGT.finalize_vgt_products(q, u, Float64.(I), KER).psi, Float64.(I)
end

# rho_g per 20 px block (median over blocks); compute it the same way on the target cube.
function rho_g(win_raw, sigma, two_sg2)
    nx, ny, nc = size(win_raw)
    num = zeros(nx, ny)
    cnt = zeros(nx, ny)
    for k in 1:nc
        c = @view win_raw[:, :, k]
        for j in 2:ny-1, i in 2:nx-1
        all(c[ii, jj] >= 3sigma for (ii, jj) in ((i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1))) || continue
        gx = (c[i-1, j] - c[i+1, j]) / 2
        gy = (c[i, j-1] - c[i, j+1]) / 2
        num[i, j] += gx^2 + gy^2
        cnt[i, j] += 1
        end
    end
    vals = Float64[]
    for bj in 0:(ny ÷ DN - 1), bi in 0:(nx ÷ DN - 1)
        r = (bi*DN+1):(bi+1)*DN
        s = (bj*DN+1):(bj+1)*DN
        m = cnt[r, s] .> 0
        any(m) || continue
        push!(vals, sqrt(max(mean(num[r, s][m] ./ cnt[r, s][m]) / two_sg2 - 1, 0.0)))
    end
    return isempty(vals) ? NaN : median(vals)
end

function block_rows(io, key, layer, scale, am_maps, common)
    nx, ny = size(common)
    for (variant, cmap) in am_maps
        for bj in 0:(ny ÷ BLOCK - 1), bi in 0:(nx ÷ BLOCK - 1)
            r = (bi*BLOCK+1):(bi+1)*BLOCK
            s = (bj*BLOCK+1):(bj+1)*BLOCK
            m = common[r, s]
            n = count(m)
            am = n > 0 ? mean(cmap[r, s][m]) : NaN
            println(io, join((key..., layer, scale, variant, bj * (nx ÷ BLOCK) + bi, @sprintf("%.6f", am), n), ","))
        end
    end
end

function main(args)
    o = opts_from(args)
    out = o["out"]
    mkpath(joinpath(out, "maps"))
    runs = split(o["runs"], ",")
    loss = parse.(Int, split(o["los"], ","))
    snrs = parse.(Float64, split(o["snr"], ","))
    nseed = parse(Int, o["nseed"])
    nv = parse(Int, o["nv"])
    fb = open(joinpath(out, "blocks.csv"), "w")
    fw = open(joinpath(out, "windows.csv"), "w")
    println(fb, "run,los,snr,seed,layer,scale,variant,block,am,npix")
    println(fw, "run,los,snr,seed,variant,k0,vcen_std,rho_g,am10_window,cover")
    for run in runs, los in loss
        dir = joinpath(o["cats"], run)
        sky = Tuple(setdiff(1:3, los))
        perm = (sky..., los)
        rho = permutedims(rd(dir, "dens", "t800"), perm)
        v = permutedims(rd(dir, ("velx", "vely", "velz")[los], "t800"), perm)
        b1 = permutedims(rd(dir, ("magx", "magy", "magz")[sky[1]], "t800"), perm)
        b2 = permutedims(rd(dir, ("magx", "magy", "magz")[sky[2]], "t800"), perm)
        sv = std(v)
        vel = range(-4sv, 4sv; length=nv)
        sth = step(vel) / 2
        ppv = beam!(project_ppv(rho, v, vel, sth), BEAM)
        qt = beam!(project_ppv(rho .* (b1 .^ 2 .- b2 .^ 2), v, vel, sth), BEAM)
        ut = beam!(project_ppv(rho .* 2 .* b1 .* b2, v, vel, sth), BEAM)
        rho = v = b1 = b2 = nothing
        GC.gc()
        vstd = collect(vel) ./ sv
        k0s = [k for k in 1:nv-3 if all(abs.(vstd[k:k+3]) .<= 2.5)]
        layers = Dict("all" => k0s, "wing" => [k for k in k0s if mean(vstd[k:k+3]) < -1])
        peak = median(maximum(ppv; dims=3))
        nx, ny = size(ppv, 1), size(ppv, 2)
        # Truth per layer and scale: channels covered by the layer's windows, each once.
        truth = Dict{Tuple{String, Float64}, Matrix{Float64}}()
        for (ln, ks) in layers, sc in (S10, S120)
            chans = unique(vcat([k:k+3 for k in ks]...))
            z = dropdims(sum(qt[:, :, chans]; dims=3) .+ im .* sum(ut[:, :, chans]; dims=3); dims=3)
            truth[(ln, sc)] = 0.5 .* angle.(psmooth(z, sc))
        end
        for seed in 1:nseed
            noise = randn(MersenneTwister(1000 * seed + los), Float32, size(ppv))
            beam!(noise, BEAM)
            noise ./= std(noise)
            gn = [((noise[i-1, j, k] - noise[i+1, j, k]) / 2)^2 + ((noise[i, j-1, k] - noise[i, j+1, k]) / 2)^2
                  for i in 2:nx-1, j in 2:ny-1, k in 1:4]
            two_sg2_unit = mean(gn)
            for snr in snrs
                sigma = peak / snr
                cube = ppv .+ Float32(sigma) .* noise
                if o["scan"] == "1"   # gradient S/N and mask coverage only, to place a target cube on the S/N axis
                    for (ln, ks) in layers
                        rg = [rho_g(cube[:, :, k:k+3], sigma, two_sg2_unit * sigma^2) for k in ks]
                        cv = [mean(any(cube[:, :, k:k+3] .>= 3sigma; dims=3)) for k in ks]
                        @printf("scan %s los %d seed %d S/N %5.1f %-4s rho_g median %.2f [%.2f, %.2f]  mask cover %.2f [%.2f, %.2f]\n",
                                run, los, seed, snr, ln, median(rg), extrema(rg)..., median(cv), extrema(cv)...)
                    end
                    continue
                end
                cube_f = beam!(copy(cube), PREFILTER)
                key = (run, los, snr, seed)
                psis = Dict{Tuple{String, Int}, Matrix{Float64}}()
                cols = Dict{Int, Matrix{Float64}}()
                for variant in VARIANTS, k in k0s
                    raw = cube[:, :, k:k+3]
                    grad = variant in ("pre4", "amp2pre4") ? cube_f[:, :, k:k+3] : raw
                    psis[(variant, k)], cols[k] = window_psi(raw, grad, sigma, variant)
                end
                for k in k0s
                    rg = rho_g(cube[:, :, k:k+3], sigma, two_sg2_unit * sigma^2)
                    zt = dropdims(sum(qt[:, :, k:k+3]; dims=3) .+ im .* sum(ut[:, :, k:k+3]; dims=3); dims=3)
                    pt = 0.5 .* angle.(psmooth(zt, S10))
                    for variant in VARIANTS
                        z = stokes_at(psis[(variant, k)], cols[k], S10)
                        ok = isfinite.(z)
                        am = count(ok) > 0 ? mean(cos.(2 .* (0.5 .* angle.(z[ok]) .- pt[ok]))) : NaN
                        println(fw, join((key..., variant, k, @sprintf("%.3f", mean(vstd[k:k+3])),
                                          @sprintf("%.4f", rg), @sprintf("%.5f", am),
                                          @sprintf("%.4f", count(isfinite, psis[(variant, k)]) / (nx * ny))), ","))
                    end
                end
                for (ln, ks) in layers
                    isempty(ks) && continue
                    d10 = Dict{String, Matrix{Float64}}()
                    d120 = Dict{String, Matrix{Float64}}()
                    for variant in VARIANTS
                        zs = [stokes_at(psis[(variant, k)], cols[k], S10) for k in ks]
                        us = [unitc.(stokes_at(psis[(variant, k)], ones(nx, ny), S120)) for k in ks]
                        d10[variant] = 0.5 .* angle.(sum(zs))
                        d120[variant] = 0.5 .* angle.(sum(abs.(z) .* u for (z, u) in zip(zs, us)))
                    end
                    for (sc, d) in ((S10, d10), (S120, d120))
                        common = reduce(.&, [isfinite.(m) for m in values(d)])
                        t = truth[(ln, sc)]
                        am_maps = [(var, cos.(2 .* (d[var] .- t))) for var in VARIANTS]
                        block_rows(fb, key, ln, Int(sc), am_maps, common)
                    end
                    if o["maps"] == "1" && seed == 1 && ln == "wing" && snr in (12.0, 23.0)
                        p = joinpath(out, "maps", "$(run)_los$(los)_snr$(Int(snr))_wing.fits")
                        isfile(p) && rm(p)
                        FITS(p, "w") do f
                            write(f, Float32.(truth[(ln, S10)]); name="TRUE10")
                            for var in VARIANTS
                                write(f, Float32.(d10[var]); name=uppercase(var) * "10")
                            end
                            write(f, Float32.(dropdims(sum(cube[:, :, layers["wing"][1]:layers["wing"][end]+3]; dims=3); dims=3)); name="COLUMN")
                        end
                    end
                end
                flush(fb)
                flush(fw)
                @info "done" run los seed snr windows=length(k0s) wing=length(layers["wing"])
            end
        end
        ppv = qt = ut = nothing
        GC.gc()
    end
    close(fb)
    close(fw)
end

main(ARGS)
