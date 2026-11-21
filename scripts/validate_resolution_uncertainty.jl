using DrWatson
quickactivate(@__DIR__, "VGT")

using FFTW
using FITSIO
using Images
using Printf
using Random
using Statistics
using VGT

include(joinpath(@__DIR__, "synthetic_orientation_fields.jl"))

# Effective resolution, window choice and noise uncertainty of the VGT angle map, on synthetic cubes with a known
# orientation field. Declared rules and results: docs/validation-resolution-uncertainty.md.
#
# Why: VGT averages gradient orientations over a dn x dn window centred on every pixel (stride 1) and then smooths
# the Stokes maps with a Gaussian of sigma `ker` px, so the angle map is coarser than the cube. Comparisons with
# Planck or starlight must be made at a common resolution, which needs the effective resolution measured, not
# assumed; and the map needs a per-pixel noise error bar.
#
# Part tf (transfer function): orientation field t = t0 + a sin(2 pi x / lambda) (t0 30 deg, a 15 deg, so the
#   doubled angle stays in the linear range), lambda = 512/m for m in 2..40 (periodic on the 512 px grid). Recovered
#   amplitude of the sinusoid in the column-averaged (Stokes) VGT angle over the interior, divided by a: T(k).
#   The effective FWHM is the Gaussian FWHM with the same half-transfer wavelength: FWHM = lambda_half * 2 ln 2 / pi
#   (= lambda_half / 2.266), lambda_half where T crosses 0.5. A Gaussian exp(-k^2 s^2 / 2) fitted to T between 0.1
#   and 0.95 (the declared primary estimator) is kept as the secondary column fwhm_gaussfit_px: the box window has a
#   sinc response whose sidelobes beyond the first zero (T 0.1-0.2) fall inside the fit range and bias it low for
#   dn >= 20 at ker 1-2 and dn >= 32 (replaced after the full run, docs/validation-resolution-uncertainty.md).
#   The beam is 3 px; noise is beam-correlated (FWHM 3 px) at S/N r. Every dn and ker, plus rbar2_weight at dn 12,
#   20, 32.
# Part sweep (accuracy against effective resolution): random orientation fields (random_tmap, delta 0.8, power
#   k^-3) and, if --cats-dirs is given, CATS MHD cubes from scripts/synthetic_ppv_from_mhd.jl (ppv.fits,
#   psi_true.fits, beam 3 px, no noise; noise added here at S/N r with a 3 sigma mask). AM = <cos 2(psi - truth)>
#   against the pixel-scale truth and, for each target FWHM F (px), after Stokes-smoothing the VGT map from its
#   measured FWHM (part tf, r = Inf) to F and the truth (unit Stokes) to F. Configs whose FWHM exceeds F are
#   scored as they are and flagged.
# Part mc (noise uncertainty): psi_noise_mc (src/science/uncertainty.jl) on synthetic signal + one noise
#   realisation n0; error = psi(signal + n0) - psi(signal) (the noise-free VGT map, not the truth, so method error
#   is not counted as noise); coverage = fraction of pixels with |error| <= sigma_psi (68 % if calibrated) and
#   <= 2 sigma_psi (95 %), and rms(error / sigma_psi); a 32 px block bootstrap gives the 68 % interval of the
#   1-sigma coverage. scale 1 and 0.5 (added noise equal to / half of the cube's).
#
#   julia --project=. -t 16 scripts/validate_resolution_uncertainty.jl --part tf|sweep|mc --out DIR
#       [--seeds 3] [--snr Inf,2,1,0.5] [--dn 8,12,16,20,32,48] [--ker 1,2,4] [--targets 10.5,15,21,30]
#       [--cats-dirs DIR1,DIR2] [--tf-fit DIR/tf_fit.csv] [--realizations 16] [--quick 1]
# Writes DIR/<part>_runs.csv (and tf_fit.csv for part tf). --part fit refits DIR/tf_runs.csv into DIR/tf_fit.csv.

const N = 512
const BEAM = 3.0
const MARGIN = 40          # common interior: >= dn/2 + 2 ker + 8 for the largest window (48) and ker 4

function parse_args(args)
    o = Dict("part" => "tf", "out" => "validation_output/resolution", "seeds" => "3", "snr" => "Inf,2,1,0.5",
             "dn" => "8,12,16,20,32,48", "ker" => "1,2,4", "targets" => "10.5,15,21,30", "cats-dirs" => "",
             "tf-fit" => "", "realizations" => "16", "quick" => "0", "noise-elong" => "1.0", "noise-angle" => "0")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("bad argument $(args[i])"))
        o[args[i][3:end]] = args[i + 1]
        i += 2
    end
    return o
end
floats(s) = [x == "Inf" ? Inf : parse(Float64, x) for x in split(s, ",") if !isempty(x)]
ints(s) = [parse(Int, x) for x in split(s, ",") if !isempty(x)]
axial(a) = mod(a + pi / 2, pi) - pi / 2

function gsmooth(a::AbstractMatrix, s::Real)
    s <= 0 && return copy(a)
    l = 2 * ceil(Int, 3s) + 1
    return imfilter(a, Kernel.gaussian((s, s), (l, l)), NA())
end

"Effective std (px) of the truncated Gaussian that finalize_vgt_products uses for `ker`."
function ker_std(ker)
    k = Kernel.gaussian(ker)
    w = vec(sum(collect(k); dims=2))
    x = collect(axes(k, 1))
    return sqrt(sum(w .* x .^ 2) / sum(w))
end

noise_cube(rng, n, nv, o) = correlated_noise(rng, n, n, nv; fwhm_major_px=BEAM * sqrt(parse(Float64, o["noise-elong"])),
                                             fwhm_minor_px=BEAM / sqrt(parse(Float64, o["noise-elong"])),
                                             angle=deg2rad(parse(Float64, o["noise-angle"])))

# VGT on a cube: one kernel pass per (dn, rbar2_weight), finalised for every ker. Returns Dict((dn, ker, rb) => psi).
function vgt_maps(cube, threshold, dns, kers; rb_dns=Int[])
    nv = size(cube, 3)
    intensity = dropdims(sum(cube; dims=3); dims=3)
    out = Dict{Tuple{Int, Int, Bool}, Matrix{Float64}}()
    for dn in dns, rb in (false, true)
        rb && !(dn in rb_dns) && continue
        q, u = VGT_QU_mv_2d(cube, dn, 1, nv, threshold; min_valid=dn^2 ÷ 4, rbar2_weight=rb)
        for ker in kers
            out[(dn, ker, rb)] = finalize_vgt_products(copy(q), copy(u), copy(intensity), ker).psi
        end
    end
    return out
end

# Signal plus level so that a 3 sigma mask never acts (synthetic GRF part; CATS keeps its own zeros).
add_level(s) = s .- minimum(s) .+ 10.0

function interior(n, m=MARGIN)
    return (m + 1):(n - m)
end

# ---------------------------------------------------------------- part tf
function sinusoid_amplitude(psi, t0, lambda)
    rng = interior(size(psi, 1))
    xs = Float64[]
    ys = Float64[]
    for i in rng
        c = 0.0
        s = 0.0
        for j in rng
            p = psi[i, j]
            isfinite(p) || continue
            c += cos(2p)
            s += sin(2p)
        end
        (c == 0 && s == 0) && continue
        push!(xs, i)
        push!(ys, axial(0.5 * atan(s, c) - t0))
    end
    X = hcat(sin.(2pi .* xs ./ lambda), cos.(2pi .* xs ./ lambda), ones(length(xs)))
    β = X \ ys
    return β[1], β[2]
end

function part_tf(o)
    dns, kers, snrs = ints(o["dn"]), ints(o["ker"]), floats(o["snr"])
    seeds = parse(Int, o["seeds"])
    ms = o["quick"] == "1" ? [4, 16] : [2, 4, 6, 8, 12, 16, 24, 32, 40]
    t0, a = deg2rad(30.0), deg2rad(15.0)
    rb_dns = intersect(dns, [12, 20, 32])
    path = joinpath(o["out"], "tf_runs.csv")
    open(path, "w") do io
        println(io, "dn,ker,rb,snr,seed,lambda,k,A,B,T")
        for r in snrs, seed in 1:seeds, m in ms
            lambda = N / m
            rng = Xoshiro(1000 * seed + m)
            tmap = sinusoid_tmap(N, t0, a, lambda)
            sig = orientation_signal(rng, tmap; r=r, beam=BEAM)
            cube = add_level(sig)
            isinf(r) || (cube .+= noise_cube(rng, N, size(cube, 3), o))
            maps = vgt_maps(cube, 3.0, dns, kers; rb_dns)
            for ((dn, ker, rb), psi) in maps
                A, B = sinusoid_amplitude(psi, t0, lambda)
                @printf(io, "%d,%d,%d,%s,%d,%.3f,%.5f,%.5f,%.5f,%.5f\n", dn, ker, rb, r, seed, lambda, 2pi / lambda,
                        A, B, hypot(A, B) / a)
            end
            flush(io)
            @info "tf" r seed lambda
        end
    end
    fit_tf(path, joinpath(o["out"], "tf_fit.csv"))
end

"""
    half_transfer(ks, T) -> (lambda_half, fwhm)

Wavelength (px) where the transfer T(k), sampled at increasing k, first falls through 0.5 (log-linear interpolation
in wavelength), and the Gaussian FWHM with the same half-transfer point, lambda_half * 2 ln 2 / pi. Only the first
crossing counts, so sinc sidelobes at higher k do not move it. NaN if T never crosses 0.5.
"""
function half_transfer(ks::AbstractVector, T::AbstractVector)
    issorted(ks) || throw(ArgumentError("ks must be increasing"))
    lam = 2pi ./ ks
    for i in 1:(length(ks) - 1)
        if T[i] >= 0.5 && T[i + 1] < 0.5
            f = (0.5 - T[i]) / (T[i + 1] - T[i])
            lh = exp(log(lam[i]) + f * (log(lam[i + 1]) - log(lam[i])))
            return lh, lh * 2log(2) / pi
        end
    end
    return NaN, NaN
end

"Gaussian FWHM fitted to log T between 0.1 and 0.95 (secondary; biased by box sidelobes, see the header)."
function gaussfit_fwhm(ks::AbstractVector, T::AbstractVector)
    sel = (T .> 0.1) .& (T .< 0.95)
    count(sel) >= 2 || return NaN
    X = hcat(-ks[sel] .^ 2 ./ 2, ones(count(sel)))
    β = X \ log.(T[sel])
    return β[1] > 0 ? sqrt(8log(2)) * sqrt(β[1]) : NaN
end

function fit_tf(runs, out)
    lines = readlines(runs)[2:end]
    rows = [split(l, ",") for l in lines]
    groups = Dict{Tuple{Int, Int, Int, String}, Vector{Tuple{Float64, Float64}}}()
    for r in rows
        key = (parse(Int, r[1]), parse(Int, r[2]), parse(Int, r[3]), r[4])
        push!(get!(groups, key, Tuple{Float64, Float64}[]), (parse(Float64, r[7]), parse(Float64, r[10])))
    end
    open(out, "w") do io
        # Column 5 (fwhm_px, half-transfer) is what part sweep reads as the measured FWHM.
        println(io, "dn,ker,rb,snr,fwhm_px,fwhm_gaussfit_px,fwhm_formula_px,lambda_half_px,T_lowk,npoints")
        for key in sort(collect(keys(groups)))
            pts = groups[key]
            ks = sort(unique(first.(pts)))
            T = [mean(last.(filter(p -> p[1] == k, pts))) for k in ks]
            lh, fw = half_transfer(ks, T)
            gf = gaussfit_fwhm(ks, T)
            dn, ker = key[1], key[2]
            formula = sqrt(8log(2)) * sqrt(dn^2 / 12 + ker_std(ker)^2 + (BEAM / sqrt(8log(2)))^2)
            @printf(io, "%d,%d,%d,%s,%.3f,%.3f,%.3f,%.3f,%.4f,%d\n", dn, ker, key[3], key[4], fw, gf, formula, lh, T[1],
                    count((T .> 0.1) .& (T .< 0.95)))
        end
    end
end

# ---------------------------------------------------------------- part sweep
function read_tf_fit(path)
    d = Dict{Tuple{Int, Int, Int}, Float64}()
    isempty(path) && return d
    ls = readlines(path)
    # Old files (before the half-transfer change) carry the sidelobe-biased Gaussian fit in column 5; refuse them.
    split(ls[1], ",")[5] == "fwhm_px" ||
        throw(ArgumentError("$path has no half-transfer fwhm_px column; refit it with --part fit"))
    for l in ls[2:end]
        r = split(l, ",")
        r[4] == "Inf" || continue
        d[(parse(Int, r[1]), parse(Int, r[2]), parse(Int, r[3]))] = parse(Float64, r[5])
    end
    return d
end

function am_scores(psi, truth, fwhm, targets)
    rg = interior(size(psi, 1))
    p = psi[rg, rg]
    t = truth[rg, rg]
    ok = isfinite.(p) .& isfinite.(t)
    am0 = mean(cos.(2 .* (p[ok] .- t[ok])))
    out = Float64[am0]
    s2 = 1 / sqrt(8log(2))
    for F in targets
        add = F > fwhm ? sqrt(F^2 - fwhm^2) * s2 : 0.0
        pq, pu = gsmooth(map(x -> isfinite(x) ? cos(2x) : 0.0, p), add), gsmooth(map(x -> isfinite(x) ? sin(2x) : 0.0, p), add)
        tq, tu = gsmooth(cos.(2 .* t), F * s2), gsmooth(sin.(2 .* t), F * s2)
        m = ceil(Int, 2 * F * s2)
        inner = (m + 1):(size(p, 1) - m)
        a = 0.5 .* atan.(pu[inner, inner], pq[inner, inner])
        b = 0.5 .* atan.(tu[inner, inner], tq[inner, inner])
        push!(out, mean(cos.(2 .* (a .- b))))
    end
    return out
end

function load_cats(dir)
    cube = Float64.(FITS(f -> read(f[1]), joinpath(dir, "ppv.fits")))
    truth = Float64.(FITS(f -> read(f[1]), joinpath(dir, "psi_true.fits")))
    # Crop / tile to N x N if needed.
    n = min(size(cube, 1), size(cube, 2), N)
    return cube[1:n, 1:n, :], truth[1:n, 1:n]
end

function hp_rms(cube)
    n = size(cube, 1)
    k = kgrid(n)
    hptf = [highpass_tf(a, b) for a in k, b in k]
    vals = [std(real.(ifft(fft(cube[:, :, v]) .* hptf))) for v in 1:size(cube, 3)]
    # S/N of the bright channels (the top half by high-passed rms), as the narrow-line regime that dominates VGT.
    return median(sort(vals)[(end ÷ 2 + 1):end])
end

function part_sweep(o)
    dns, kers, snrs = ints(o["dn"]), ints(o["ker"]), floats(o["snr"])
    targets = floats(o["targets"])
    seeds = parse(Int, o["seeds"])
    fw = read_tf_fit(o["tf-fit"])
    rb_dns = intersect(dns, [12, 20, 32])
    fields = Any[("grf", nothing)]
    for d in split(o["cats-dirs"], ",")
        isempty(d) || push!(fields, ("cats:" * basename(rstrip(d, '/')), d))
    end
    path = joinpath(o["out"], "sweep_runs.csv")
    open(path, "w") do io
        println(io, "field,snr,seed,dn,ker,rb,fwhm_px,am_native," * join(["am_F$(F)" for F in targets], ","))
        for (name, dir) in fields, r in snrs, seed in 1:seeds
            (isinf(r) && seed > 1 && dir !== nothing) && continue  # noise-free CATS is deterministic
            rng = Xoshiro(7000 + 31 * seed + round(Int, isinf(r) ? 999 : 10r))
            if dir === nothing
                truth = random_tmap(rng, N)
                cube = add_level(orientation_signal(rng, truth; r=r, beam=BEAM))
                isinf(r) || (cube .+= noise_cube(rng, N, size(cube, 3), o))
                thr = 3.0
            else
                cube, truth = load_cats(dir)
                sigma = isinf(r) ? 0.0 : hp_rms(cube) / r
                isinf(r) || (cube = cube .+ sigma .* noise_cube(rng, size(cube, 1), size(cube, 3), o))
                thr = isinf(r) ? 1e-6 * maximum(cube) : 3sigma
            end
            maps = vgt_maps(cube, thr, dns, kers; rb_dns)
            for ((dn, ker, rb), psi) in maps
                f = get(fw, (dn, ker, Int(rb)), NaN)
                isfinite(f) || (f = sqrt(8log(2)) * sqrt(dn^2 / 12 + ker_std(ker)^2 + (BEAM / sqrt(8log(2)))^2))
                sc = am_scores(psi, truth, f, targets)
                println(io, join(Any[name, r, seed, dn, ker, Int(rb), round(f; digits=3),
                                     (round.(sc; digits=5))...], ","))
            end
            flush(io)
            @info "sweep" name r seed
        end
    end
end

# ---------------------------------------------------------------- part mc
function block_coverage(z, ok; B=32, nboot=400, rng=Xoshiro(3))
    n = size(z, 1)
    nb = n ÷ B
    blocks = [(count(ok[(bi - 1) * B + 1:bi * B, (bj - 1) * B + 1:bj * B] .& (abs.(z[(bi - 1) * B + 1:bi * B, (bj - 1) * B + 1:bj * B]) .<= 1)),
               count(ok[(bi - 1) * B + 1:bi * B, (bj - 1) * B + 1:bj * B])) for bi in 1:nb, bj in 1:nb]
    blocks = vec(blocks)
    vals = Float64[]
    for _ in 1:nboot
        s = rand(rng, blocks, length(blocks))
        push!(vals, sum(first.(s)) / max(sum(last.(s)), 1))
    end
    return quantile(vals, 0.16), quantile(vals, 0.84)
end

function part_mc(o)
    snrs = floats(o["snr"])
    dns = ints(o["dn"])
    seeds = parse(Int, o["seeds"])
    K = parse(Int, o["realizations"])
    fields = Any[("grf", nothing)]
    for d in split(o["cats-dirs"], ",")
        isempty(d) || push!(fields, ("cats:" * basename(rstrip(d, '/')), d))
    end
    path = joinpath(o["out"], "mc_runs.csv")
    open(path, "w") do io
        println(io, "field,snr,seed,dn,scale,npix,cov1,cov1_lo,cov1_hi,cov2,rms_z,median_sigma_deg,rms_err_deg," *
                    "mean_err_deg,median_offset_deg,cov1_tot,cov2_tot,rms_z_tot,bin_sigma_deg,bin_rms_err_deg")
        for (name, dir) in fields, r in snrs, seed in 1:seeds
            isinf(r) && continue
            rng = Xoshiro(9000 + 17 * seed + round(Int, 10r))
            if dir === nothing
                tmap = random_tmap(rng, N)
                sig = add_level(orientation_signal(rng, tmap; r=r, beam=BEAM))
                sigma, thr = 1.0, 3.0
            else
                sig, _ = load_cats(dir)
                sigma = hp_rms(sig) / r
                thr = 3sigma
            end
            n, nv = size(sig, 1), size(sig, 3)
            el = parse(Float64, o["noise-elong"])
            nkw = (; fwhm_major_px=BEAM * sqrt(el), fwhm_minor_px=BEAM / sqrt(el),
                   noise_angle=deg2rad(parse(Float64, o["noise-angle"])))
            obs = sig .+ correlated_noise(rng, n, n, nv; fwhm_major_px=nkw.fwhm_major_px,
                                          fwhm_minor_px=nkw.fwhm_minor_px, angle=nkw.noise_angle, sigma=sigma)
            for dn in dns
                q, u = VGT_QU_mv_2d(sig, dn, 1, nv, thr; min_valid=dn^2 ÷ 4)
                clean = finalize_vgt_products(q, u, dropdims(sum(sig; dims=3); dims=3), 4).psi
                for scale in (1.0, 0.5)
                    res = psi_noise_mc(obs, dn, 4, thr; realizations=K, noise_sigma=sigma, scale=scale,
                                       rng=Xoshiro(seed * 100 + dn), min_valid=dn^2 ÷ 4, nkw...)
                    rg = interior(n)
                    err = axial.(res.psi[rg, rg] .- clean[rg, rg])
                    sg = res.sigma_psi[rg, rg]
                    ok = isfinite.(err) .& isfinite.(sg) .& (sg .> 0)
                    z = err ./ sg
                    lo, hi = block_coverage(z, ok)
                    e, s = err[ok], sg[ok]
                    # Calibration curve: rms error in 5 quantile bins of sigma_psi.
                    qs = quantile(s, 0:0.2:1)
                    bs, be = String[], String[]
                    for b in 1:5
                        sel = (s .>= qs[b]) .& (s .<= qs[b + 1])
                        push!(bs, @sprintf("%.2f", rad2deg(median(s[sel]))))
                        push!(be, @sprintf("%.2f", rad2deg(sqrt(mean(e[sel] .^ 2)))))
                    end
                    off = res.offset[rg, rg][ok]
                    # Post hoc variant (added after the smoke test, before the full runs): sigma combined with the
                    # Monte Carlo mean shift, sqrt(sigma^2 + offset^2), as a total error including the noise bias.
                    zt = e ./ sqrt.(s .^ 2 .+ off .^ 2)
                    @printf(io, "%s,%s,%d,%d,%.2f,%d,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.3f,%s,%s\n", name, r,
                            seed, dn, scale, count(ok), mean(abs.(z[ok]) .<= 1), lo, hi, mean(abs.(z[ok]) .<= 2),
                            sqrt(mean(z[ok] .^ 2)), rad2deg(median(s)), rad2deg(sqrt(mean(e .^ 2))),
                            rad2deg(mean(e)), rad2deg(median(off)), mean(abs.(zt) .<= 1), mean(abs.(zt) .<= 2),
                            sqrt(mean(zt .^ 2)), join(bs, ";"), join(be, ";"))
                    flush(io)
                end
                @info "mc" name r seed dn
            end
        end
    end
end

function main(args)
    o = parse_args(args)
    mkpath(o["out"])
    o["part"] == "tf" && part_tf(o)
    o["part"] == "sweep" && part_sweep(o)
    o["part"] == "mc" && part_mc(o)
    o["part"] == "fit" && fit_tf(joinpath(o["out"], "tf_runs.csv"), joinpath(o["out"], "tf_fit.csv"))
end

# Guarded so tests can include the file for its functions without running a part.
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
