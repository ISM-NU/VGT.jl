using DrWatson
quickactivate(@__DIR__, "VGT")

using VGT
using FFTW
using FITSIO
using Images: imfilter, KernelFactors, distance_transform, feature_transform
using Printf
using Random
using Statistics

# VGT behaviour study: how the pipeline responds to noise, smoothing, thresholds, Stokes
# weights, NaNs, survey edges and mosaic seams, measured on MHD cubes with a known field.
# Results: docs/validation-behaviour-2026-11-22.md.
#
# Declared before running (do not change after seeing results):
# - Truth: plane-of-sky field direction of the MHD cube, psi_true = 0.5 atan2(U, Q) with
#   Q = sum rho (B1^2 - B2^2), U = sum rho 2 B1 B2, beam-smoothed and passed through the same
#   smoothing and binning as the cube ("T_beam"). A second truth "T_prod" smooths T_beam's
#   Q/U further to FWHM 16 px of the analysis grid (about the effective resolution of
#   dn 20 / ker 4); it is reported where settings change the resolution.
# - Metrics per map, on the evaluation pixels (interior, 16 px from the array border on the
#   analysis grid, finite in the map): AM = <cos 2(psi - psi_true)>; bias = the doubled-angle
#   mean of psi - psi_true (deg); scatter = rms of the headless difference about the bias
#   (deg); coverage = finite fraction of the evaluation pixels. Paired differences use the
#   pixels finite in both maps, the same field view and the same noise seed.
# - Errors: mean over field views x noise seeds; the standard error is over those units
#   (4 views from 2 sub-Alfvenic runs seen along 2 axes each, or 2 views of a super-Alfvenic
#   run). Views of the same run are not independent fields.
# - Noise: per-pixel rms fixed on the native grid by the peak S/N (median over the map of the
#   peak brightness / noise rms); the gradient S/N on the analysis grid is reported too.
#
# Grids: the native cube has a 3 px Gaussian beam (3 px per beam). The production chain
# ("2' chain" in project terms; here "chain2") smooths every channel to a 6 px beam with
# NaN-aware normalised convolution (output kept where the kernel weight on finite pixels
# exceeds `cut`, 0.5 as the project's smoothing scripts) and averages 2 x 2 blocks of finite
# pixels, so the analysis grid again has 3 px per beam. VGT on it: dn 20, ker 4,
# min_valid 100, central_checked, axial_vonmises, threshold = 1 x the noise rms of the
# analysis grid, smoothed_stokes_support, nan_before_smoothing. The 256^3 CATS boxes are
# periodic; each view is tiled 2 x 2 (512^2 native, 256^2 analysis grid) so windows near the
# middle see continuous gas.
#
# Experiments (--exp):
#   noise    white / spatially correlated grain (FWHM 3, 6 px native) / grain + channel
#            correlation (lag-1 0.25, 0.5) / anisotropic grain, at peak S/N 2..64
#   split    split-data reproducibility under channel-correlated noise: odd/even channels,
#            alternating 2-channel blocks, every 4th channel (gap) halves
#   chains   smoothing before VGT (chain2) against Q/U smoothing after (native grid, dn 40,
#            ker 8 = the same sky kernel, and ker 9 = about the same total resolution), smoothing without binning, and gradient whitening on both grids
#   thr      noise threshold multiplier (none, 0.5, 1, 2, 3) and min_valid (50..200)
#   weights  Stokes weight of a window: centre intensity (default), clipped at 0, squared,
#            square root, uniform (1 where the centre passes the threshold)
#   support  output support and zero-Stokes policies, and the error of filled pixels
#            against their distance from the nearest direct fit
#   dn       window size 12, 20, 32 (min_valid dn^2 / 4)
#   edges    survey edge (half plane) across or along the mean field: NaN outside with the
#            smoothing cut at 0.5 or 0.99, NaN applied after smoothing, zero outside, and a
#            noise rise toward the edge with a scalar or per-pixel threshold; edge_mask off/on
#   wsplit   split-data reproducibility (odd/even, gap halves) of the default, square-root and
#            uniform Stokes weights, against their accuracy
#   seams    one region where a block of channels (blue line wing holding 15 % of the flux,
#            or a quarter of the channels around the peak) fades by a gain g over
#            a ramp of L px, or is NaN; T_seam = <cos 2(psi - seam tangent)> minus the same
#            pixels without the seam
#
#   julia --project=. -t 32 scripts/validate_behaviour.jl --exp noise --cats-dir DIR \
#       --out validation_output/behaviour [--seeds 3] [--quick 1]
# Writes <out>/<exp>.csv in long format: exp,cond,snr,view,seed,metric,value.

const BEAM = 3.0                    # native beam FWHM (px)
const TARGET = 6.0                  # chain2 beam FWHM on the native grid (px)
const KER_FWHM = sqrt(TARGET^2 - BEAM^2)
const PROD_FWHM = 16.0              # T_prod smoothing on the analysis grid (px)
const BORDER = 16                   # evaluation border on the analysis grid (px)
const VIEWS_SUB = [("b1p1", 3), ("b1p1", 2), ("b1p.32", 3), ("b1p.32", 2)]   # sub-Alfvenic, sub- to transonic
const VIEWS_STR = [("b.1p.032", 3), ("b.1p.032", 2)]

function parse_args(args)
    o = Dict("seeds" => "3", "quick" => "0", "nv" => "64", "sigma-th" => "1.0", "out" => "validation_output/behaviour")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("Bad argument: $(args[i])"))
        o[args[i][3:end]] = args[i + 1]
        i += 2
    end
    haskey(o, "exp") && haskey(o, "cats-dir") || throw(ArgumentError("--exp and --cats-dir are required"))
    return o
end

# ------------------------------------------------------------------ smoothing helpers

# Periodic (elliptical) Gaussian smoothing of every slice along dims 1, 2. fwhm is the
# geometric-mean FWHM; ratio the axis ratio; angle (rad) the long-axis direction from axis 1.
function gsmooth_periodic!(a::AbstractArray, fwhm; ratio=1.0, angle=0.0)
    fwhm > 0 || return a
    s = fwhm / sqrt(8log(2))
    s1, s2 = s * sqrt(ratio), s / sqrt(ratio)
    c, sn = cos(angle), sin(angle)
    sxx = c^2 * s1^2 + sn^2 * s2^2
    syy = sn^2 * s1^2 + c^2 * s2^2
    sxy = c * sn * (s1^2 - s2^2)
    nx, ny = size(a, 1), size(a, 2)
    kx, ky = fftfreq(nx), fftfreq(ny)
    g = [exp(-2pi^2 * (sxx * x^2 + syy * y^2 + 2sxy * x * y)) for x in kx, y in ky]
    for k in axes(a, 3)
        v = view(a, :, :, k)
        v .= real.(ifft(fft(v) .* g))
    end
    return a
end

# NaN-aware normalised Gaussian convolution, kept where the weight on finite pixels exceeds
# `cut` (the project's smoothing scripts use 0.5); circular at the array border (tiled box).
function smooth_nan(cube::AbstractArray{<:Real, 3}, fwhm; cut=0.5)
    s = fwhm / sqrt(8log(2))
    k = KernelFactors.gaussian((s, s))
    out = similar(cube, Float64)
    for c in axes(cube, 3)
        p = view(cube, :, :, c)
        m = isfinite.(p)
        num = imfilter(ifelse.(m, p, 0.0), k, "circular")
        den = imfilter(Float64.(m), k, "circular")
        out[:, :, c] .= ifelse.(den .> cut, num ./ max.(den, 1e-12), NaN)
    end
    return out
end

# Mean of the finite pixels of each b x b block (NaN if none), every slice.
function bin_nan(a::AbstractArray, b::Int)
    nx, ny = size(a, 1) ÷ b, size(a, 2) ÷ b
    out = fill(NaN, nx, ny, size(a, 3))
    @inbounds for k in axes(a, 3), j in 1:ny, i in 1:nx
        t, n = 0.0, 0
        for jj in (j - 1) * b + 1:j * b, ii in (i - 1) * b + 1:i * b
            v = a[ii, jj, k]
            isfinite(v) || continue
            t += v
            n += 1
        end
        n > 0 && (out[i, j, k] = t / n)
    end
    return ndims(a) == 2 ? out[:, :, 1] : out
end

chain2(cube; cut=0.5) = bin_nan(smooth_nan(cube, KER_FWHM; cut), 2)

# ------------------------------------------------------------------ field views

function project_ppv(rho, v, vel, sth)
    nx, ny, nz = size(rho)
    nv = length(vel)
    dv = step(vel)
    vmin = first(vel)
    h = ceil(Int, 4sth / dv)
    ppv = zeros(nx, ny, nv)
    Threads.@threads for y in 1:ny
        for x in 1:nx, z in 1:nz
            r, w = rho[x, y, z], v[x, y, z]
            k0 = clamp(round(Int, (w - vmin) / dv) + 1, 1, nv)
            for k in max(1, k0 - h):min(nv, k0 + h)
                ppv[x, y, k] += r * exp(-0.5 * ((vel[k] - w) / sth)^2) / (sqrt(2pi) * sth) * dv
            end
        end
    end
    return ppv
end

# Noiseless native cube (tiled 2 x 2, beam BEAM) and the true Q, U on the native and the
# analysis grid. Same physics as scripts/synthetic_ppv_from_mhd.jl.
function load_view(catsdir, run, los; nv=64, sth=1.0)
    dir = joinpath(catsdir, run)
    snap = replace(first(filter(n -> startswith(n, "dens_"), readdir(dir))), "dens_" => "", ".fits" => "")
    rd(name) = Float64.(FITS(f -> read(f[1]), joinpath(dir, "$(name)_$(snap).fits")))
    sky = Tuple(setdiff(1:3, los))
    perm = (sky..., los)
    rho = permutedims(rd("dens"), perm)
    v = permutedims(rd(("velx", "vely", "velz")[los]), perm)
    b1 = permutedims(rd(("magx", "magy", "magz")[sky[1]]), perm)
    b2 = permutedims(rd(("magx", "magy", "magz")[sky[2]]), perm)
    vmax = 4std(v) + 4sth
    vel = range(-vmax, vmax; length=nv)
    ppv = project_ppv(rho, v, vel, max(sth, step(vel) / 2))
    q = sum(rho .* (b1 .^ 2 .- b2 .^ 2); dims=3)
    u = sum(rho .* 2 .* b1 .* b2; dims=3)
    rho = v = b1 = b2 = nothing
    gsmooth_periodic!(ppv, BEAM)
    gsmooth_periodic!(q, BEAM)
    gsmooth_periodic!(u, BEAM)
    sig = repeat(ppv, 2, 2, 1)
    tq, tu = repeat(q, 2, 2, 1), repeat(u, 2, 2, 1)
    tq2 = bin_nan(gsmooth_periodic!(copy(tq), KER_FWHM), 2)
    tu2 = bin_nan(gsmooth_periodic!(copy(tu), KER_FWHM), 2)
    pq = gsmooth_periodic!(reshape(copy(tq2), size(tq2)..., 1), sqrt(PROD_FWHM^2 - BEAM^2))[:, :, 1]
    pu = gsmooth_periodic!(reshape(copy(tu2), size(tu2)..., 1), sqrt(PROD_FWHM^2 - BEAM^2))[:, :, 1]
    peak = median(maximum(sig; dims=3))
    return (; name="$(run)_los$(los)", sig, tq=tq[:, :, 1], tu=tu[:, :, 1], tq2, tu2, pq, pu, peak)
end

# ------------------------------------------------------------------ noise

lag1_tap(rho) = rho <= 0 ? 0.0 : (1 - sqrt(1 - 2rho^2)) / (2rho)

# Unit-rms noise: white, convolved with an (elliptical) Gaussian grain of FWHM `grain`
# (periodic, on the full native tile), then a 3-tap [a, 1, a] filter along the channels
# with lag-1 correlation `lag1`.
function make_noise(rng, nx, ny, nv; grain=3.0, ratio=1.0, angle=0.0, lag1=0.0)
    n = randn(rng, nx, ny, nv + 2)
    grain > 0 && gsmooth_periodic!(n, grain; ratio, angle)
    a = lag1_tap(lag1)
    m = Array{Float64}(undef, nx, ny, nv)
    @inbounds for k in 1:nv
        @views m[:, :, k] .= a .* n[:, :, k] .+ n[:, :, k + 1] .+ a .* n[:, :, k + 2]
    end
    m ./= std(m)
    return m
end

# ------------------------------------------------------------------ VGT and scoring

function vgt(cube, thr; dn=20, ker=4, min_valid=dn^2 ÷ 4, gm=:central_checked, whiten=false, cov=Float64[],
             weight=nothing, edge_mask=false, support=:smoothed_stokes_support, zero=:nan_before_smoothing, bin=1)
    q, u = VGT_QU_mv_2d(cube, dn, 1, size(cube, 3), thr; min_valid, gradient_mode=gm, gradient_whiten=whiten,
                        gradient_noise_cov=cov, weight_cube=weight, edge_mask)
    I = dropdims(sum(cube; dims=3); dims=3)
    fin = finalize_vgt_products(copy(q), copy(u), I, ker, zero, support; edge_mask)
    qs, us = fin.q_smooth, fin.u_smooth
    if bin > 1
        ok = isfinite.(qs) .& isfinite.(us)
        qs, us = block_mean(qs, bin; valid=ok), block_mean(us, bin; valid=ok)
    end
    psi = 0.5 .* atan.(us, qs) .+ pi / 2
    return (; psi, q, u, rel=fin.reliability)
end

wrapd(d) = mod(d + pi / 2, pi) - pi / 2

function evalmask(n)
    m = falses(n, n)
    m[BORDER + 1:n - BORDER, BORDER + 1:n - BORDER] .= true
    return m
end

# AM, bias, scatter, coverage of psi against the truth Q/U on the pixels `sel`.
function score(psi, tq, tu, sel)
    n = count(sel)
    s = c = 0.0
    ds = Float64[]
    @inbounds for i in eachindex(psi)
        sel[i] || continue
        isfinite(psi[i]) || continue
        d = wrapd(psi[i] - 0.5 * atan(tu[i], tq[i]))
        push!(ds, d)
        s += sin(2d)
        c += cos(2d)
    end
    isempty(ds) && return (; AM=NaN, bias=NaN, scatter=NaN, cov=0.0)
    b = 0.5 * atan(s, c)
    return (; AM=c / length(ds), bias=rad2deg(b), scatter=rad2deg(sqrt(mean(wrapd.(ds .- b) .^ 2))), cov=length(ds) / n)
end

# AM of a against b (two maps) on the common finite pixels of sel.
function am_pair(a, b, sel)
    t, n = 0.0, 0
    @inbounds for i in eachindex(a)
        (sel[i] && isfinite(a[i]) && isfinite(b[i])) || continue
        t += cos(2wrapd(a[i] - b[i]))
        n += 1
    end
    return n == 0 ? NaN : t / n
end

# Paired AM difference of maps a and r against the truth on their common finite pixels.
function dam_pair(a, r, tq, tu, sel)
    common = sel .& isfinite.(a) .& isfinite.(r)
    return score(a, tq, tu, common).AM - score(r, tq, tu, common).AM
end

function emit(io, exp, cond, snr, view, seed, nt::NamedTuple)
    for (k, v) in pairs(nt)
        @printf(io, "%s,%s,%g,%s,%d,%s,%.6g\n", exp, cond, snr, view, seed, k, v)
    end
    flush(io)
end

# Gradient S/N: median |grad signal| over voxels above the threshold / rms of one noise
# gradient component (central differences, analysis grid).
function gradient_snr(s2, n2, thr)
    gs = Float64[]
    gn = Float64[]
    nx, ny, nv = size(s2)
    @inbounds for k in 1:nv, j in 2:ny - 1, i in 2:nx - 1
        gx, gy = (s2[i - 1, j, k] - s2[i + 1, j, k]) / 2, (s2[i, j - 1, k] - s2[i, j + 1, k]) / 2
        s2[i, j, k] >= thr && push!(gs, hypot(gx, gy))
        push!(gn, (n2[i - 1, j, k] - n2[i + 1, j, k]) / 2)
    end
    return isempty(gs) ? NaN : median(gs) / std(gn)
end

score_nt(r, v, sel) = begin
    a = score(r.psi, v.tq2, v.tu2, sel)
    p = score(r.psi, v.pq, v.pu, sel)
    (; AM=a.AM, bias=a.bias, scatter=a.scatter, cov=a.cov, AM_prod=p.AM)
end

# ------------------------------------------------------------------ experiments

const NOISE_TYPES = [
    ("white", (; grain=0.0)),
    ("grain3", (; grain=3.0)),
    ("grain6", (; grain=6.0)),
    ("grain3_lag0.25", (; grain=3.0, lag1=0.25)),
    ("grain3_lag0.5", (; grain=3.0, lag1=0.5)),
    ("grain3_aniso", (; grain=3.0, ratio=1.05, angle=pi / 4)),
]

function exp_noise(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    sel = evalmask(size(s2, 1))
    maps = Dict{Tuple{String, Float64}, Matrix{Float64}}()
    for (tname, kw) in NOISE_TYPES   # same seed for every type: paired draws
        n2 = chain2(make_noise(MersenneTwister(1000seed + 1), nx, ny, nv; kw...))
        for snr in snrs
            sigma = v.peak / snr
            thr = sigma * std(n2)
            r = vgt(s2 .+ sigma .* n2, thr)
            maps[(tname, snr)] = r.psi
            emit(io, "noise", tname, snr, v.name, seed, (; score_nt(r, v, sel)..., gsnr=gradient_snr(s2, sigma .* n2, thr)))
        end
    end
    for (tname, _) in NOISE_TYPES, snr in snrs
        tname == "grain3" && continue
        emit(io, "noise", tname, snr, v.name, seed,
             (; dAM_vs_grain3=dam_pair(maps[(tname, snr)], maps[("grain3", snr)], v.tq2, v.tu2, sel)))
    end
end

function exp_split(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    sel = evalmask(size(s2, 1))
    halves = Dict(
        "oddeven" => (1:2:nv, 2:2:nv),
        "block2" => (filter(k -> mod(k - 1, 4) < 2, 1:nv), filter(k -> mod(k - 1, 4) >= 2, 1:nv)),
        "gap" => (1:4:nv, 3:4:nv),
    )
    for lag in (0.0, 0.25, 0.5)
        n2 = chain2(make_noise(MersenneTwister(1000seed + 2), nx, ny, nv; grain=3.0, lag1=lag))
        for snr in snrs
            sigma = v.peak / snr
            cube = s2 .+ sigma .* n2
            thr = sigma * std(n2)
            full = vgt(cube, thr)
            emit(io, "split", "lag$(lag)_full", snr, v.name, seed, score_nt(full, v, sel))
            for (hname, (h1, h2)) in halves
                r1 = vgt(cube[:, :, h1], thr)
                r2 = vgt(cube[:, :, h2], thr)
                a1 = score(r1.psi, v.tq2, v.tu2, sel)
                a2 = score(r2.psi, v.tq2, v.tu2, sel)
                emit(io, "split", "lag$(lag)_$(hname)", snr, v.name, seed,
                     (; AM_h1=a1.AM, AM_h2=a2.AM, R12=am_pair(r1.psi, r2.psi, sel),
                      AM_h1_full=am_pair(r1.psi, full.psi, sel)))
            end
        end
    end
end

function exp_chains(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    kw = (; grain=3.0, ratio=1.05, angle=pi / 4)
    n = make_noise(MersenneTwister(1000seed + 3), nx, ny, nv; kw...)
    nref = make_noise(MersenneTwister(1000seed + 99), nx, ny, 16; kw...)   # independent noise-only channels
    s2 = chain2(v.sig)
    n2 = chain2(n)
    cov1 = measure_gradient_noise(nref).cov
    cov2 = measure_gradient_noise(chain2(nref)).cov
    sm = smooth_nan(v.sig, KER_FWHM)
    nsm = smooth_nan(n, KER_FWHM)
    sel = evalmask(size(s2, 1))
    emit(io, "chains", "noise_anisotropy", 0, v.name, seed, (; A_native=gradient_anisotropy(cov1)[1], A_chain2=gradient_anisotropy(cov2)[1]))
    for snr in snrs
        sigma = v.peak / snr
        c2 = s2 .+ sigma .* n2
        t2 = sigma * std(n2)
        rA = vgt(c2, t2)
        emit(io, "chains", "2p_dn20_k4", snr, v.name, seed, score_nt(rA, v, sel))
        emit(io, "chains", "2p_dn20_k4_whiten", snr, v.name, seed,
             score_nt(vgt(c2, t2; whiten=true, cov=cov2), v, sel))
        c1 = v.sig .+ sigma .* n
        for (name, ker, w, cv) in (("1p_dn40_k8", 8, false, cov1), ("1p_dn40_k9", 9, false, cov1),
                                   ("1p_dn40_k9_whiten", 9, true, cov1), ("1p_dn40_k9_whiten_identity", 9, true, [1.0, 1.0, 0.0]))
            r = vgt(c1, sigma; dn=40, ker, min_valid=400, whiten=w, cov=cv, bin=2)
            emit(io, "chains", name, snr, v.name, seed, (; score_nt(r, v, sel)..., dAM_vs_2p=dam_pair(r.psi, rA.psi, v.tq2, v.tu2, sel)))
        end
        r = vgt(sm .+ sigma .* nsm, sigma * std(nsm); dn=40, ker=8, min_valid=400, bin=2)
        emit(io, "chains", "2p_nobin_dn40_k8", snr, v.name, seed,
             (; score_nt(r, v, sel)..., dAM_vs_2p=dam_pair(r.psi, rA.psi, v.tq2, v.tu2, sel)))
    end
end

function exp_thr(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    n2 = chain2(make_noise(MersenneTwister(1000seed + 4), nx, ny, nv; grain=3.0))
    sel = evalmask(size(s2, 1))
    for snr in snrs
        sigma = v.peak / snr
        cube = s2 .+ sigma .* n2
        s2n = sigma * std(n2)
        ref = vgt(cube, s2n)
        for k in (0.0, 0.5, 1.0, 2.0, 3.0)
            r = k == 1.0 ? ref : vgt(cube, k == 0.0 ? -1e30 : k * s2n)
            emit(io, "thr", "k$(k)", snr, v.name, seed, (; score_nt(r, v, sel)..., dAM_vs_k1=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, sel)))
        end
        for mv in (50, 140, 200)
            r = vgt(cube, s2n; min_valid=mv)
            emit(io, "thr", "minvalid$(mv)", snr, v.name, seed, (; score_nt(r, v, sel)..., dAM_vs_k1=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, sel)))
        end
    end
end

function exp_weights(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    n2 = chain2(make_noise(MersenneTwister(1000seed + 5), nx, ny, nv; grain=3.0))
    sel = evalmask(size(s2, 1))
    for snr in snrs
        sigma = v.peak / snr
        cube = s2 .+ sigma .* n2
        t = sigma * std(n2)
        ref = vgt(cube, t)
        # Fraction of the Stokes weight that is negative: centre pixels below 0 in windows that
        # the default run fits (approximated by pixels with a finite fit in that channel's mask).
        negfrac = count(<(0), cube) / length(cube)
        emit(io, "weights", "I", snr, v.name, seed, (; score_nt(ref, v, sel)..., negvox=negfrac))
        for (name, w, tt) in (("Ipos", max.(cube, 0.0), t), ("I2", max.(cube, 0.0) .^ 2, t^2),
                              ("sqrtI", sqrt.(max.(cube, 0.0)), sqrt(t)), ("unif", Float64.(cube .>= t), 0.5))
            r = vgt(cube, tt; weight=w)
            emit(io, "weights", name, snr, v.name, seed, (; score_nt(r, v, sel)..., dAM_vs_I=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, sel)))
        end
    end
end

# Does split-data reproducibility rank the Stokes weightings as their accuracy does?
function exp_wsplit(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    n2 = chain2(make_noise(MersenneTwister(1000seed + 10), nx, ny, nv; grain=3.0, lag1=0.25))
    sel = evalmask(size(s2, 1))
    halves = Dict("oddeven" => (1:2:nv, 2:2:nv), "gap" => (1:4:nv, 3:4:nv))
    for snr in snrs
        sigma = v.peak / snr
        cube = s2 .+ sigma .* n2
        t = sigma * std(n2)
        wts = Dict("I" => (nothing, t), "sqrtI" => (sqrt.(max.(cube, 0.0)), sqrt(t)), "unif" => (Float64.(cube .>= t), 0.5))
        for (wname, (w, tt)) in wts
            full = vgt(cube, tt; weight=w)
            emit(io, "wsplit", "$(wname)_full", snr, v.name, seed, score_nt(full, v, sel))
            for (hname, (h1, h2)) in halves
                sub(a, h) = a === nothing ? nothing : a[:, :, h]
                r1 = vgt(cube[:, :, h1], tt; weight=sub(w, h1))
                r2 = vgt(cube[:, :, h2], tt; weight=sub(w, h2))
                emit(io, "wsplit", "$(wname)_$(hname)", snr, v.name, seed,
                     (; AM_h1=score(r1.psi, v.tq2, v.tu2, sel).AM, AM_h2=score(r2.psi, v.tq2, v.tu2, sel).AM,
                      R12=am_pair(r1.psi, r2.psi, sel)))
            end
        end
    end
end

function exp_support(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    n2 = chain2(make_noise(MersenneTwister(1000seed + 6), nx, ny, nv; grain=3.0))
    sel = evalmask(size(s2, 1))
    for snr in snrs
        sigma = v.peak / snr
        cube = s2 .+ sigma .* n2
        t = sigma * std(n2)
        ref = vgt(cube, t)
        # distance of each pixel from the nearest pixel with a direct fit
        fit = (ref.q .!= 0) .| (ref.u .!= 0)
        dist = distance_transform(feature_transform(fit))
        a = score(ref.psi, v.tq2, v.tu2, sel)
        emit(io, "support", "smoothed_nan", snr, v.name, seed, (; a..., fitfrac=count(fit .& sel) / count(sel)))
        for (lo, hi) in ((0, 0), (1, 2), (3, 5), (6, 10), (11, 1000))
            m = sel .& (dist .>= lo) .& (dist .<= hi)
            count(m) > 50 || continue
            b = score(ref.psi, v.tq2, v.tu2, m)
            emit(io, "support", "smoothed_nan_d$(lo)-$(hi)", snr, v.name, seed, (; AM=b.AM, cov=b.cov, npix=count(m)))
        end
        for (name, sup, z) in (("direct_nan", :direct_fit_support, :nan_before_smoothing),
                               ("smoothed_keepzero", :smoothed_stokes_support, :keep_zero_before_smoothing))
            r = vgt(cube, t; support=sup, zero=z)
            emit(io, "support", name, snr, v.name, seed, (; score(r.psi, v.tq2, v.tu2, sel)..., dAM_vs_default=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, sel)))
        end
    end
end

function exp_dn(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    s2 = chain2(v.sig)
    n2 = chain2(make_noise(MersenneTwister(1000seed + 7), nx, ny, nv; grain=3.0))
    sel = evalmask(size(s2, 1))
    for snr in snrs
        sigma = v.peak / snr
        cube = s2 .+ sigma .* n2
        for dn in (12, 20, 32)
            r = vgt(cube, sigma * std(n2); dn)
            emit(io, "dn", "dn$(dn)", snr, v.name, seed, score_nt(r, v, sel))
        end
    end
end

# Distance bins from an edge or seam on the analysis grid (px).
const DBINS = ((0, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 50), (50, 80))

# T = <cos 2(psi - theta)> on pixels m.
tcos(psi, theta, m) = (vals = [cos(2(psi[i] - theta)) for i in eachindex(psi) if m[i] && isfinite(psi[i])];
                       isempty(vals) ? NaN : mean(vals))

function exp_edges(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    n = make_noise(MersenneTwister(1000seed + 8), nx, ny, nv; grain=3.0)
    s2full = chain2(v.sig)
    n2full = chain2(n)
    na = size(s2full, 1)
    x0 = 300                                  # native edge position; live side 1:x0
    for normal in (1, 2)                      # edge normal along map axis 1 or 2
        theta_n = normal == 1 ? 0.0 : pi / 2  # psi along the normal = across the edge
        live = [((normal == 1 ? i : j) <= x0) for i in 1:nx, j in 1:ny]
        dlive = [(x0 - (normal == 1 ? i : j) + 0.5) for i in 1:nx, j in 1:ny]   # native px inside
        d2 = [((x0 ÷ 2) - (normal == 1 ? i : j) + 0.5) for i in 1:na, j in 1:na] # analysis px inside
        interior = evalmask(na)
        for snr in snrs
            sigma = v.peak / snr
            t2 = sigma * std(n2full)
            ref = vgt(s2full .+ sigma .* n2full, t2)
            ramp = 1 .+ 12 .* exp.(-max.(dlive, 0) ./ 25)               # noise rise toward the edge
            cases = Dict{String, Any}()
            masked(a) = (b = copy(a); for k in 1:nv; @views b[:, :, k][.!live] .= NaN; end; b)
            zeroed(a) = (b = copy(a); for k in 1:nv; @views b[:, :, k][.!live] .= 0.0; end; b)
            cases["nan_cut0.5"] = (chain2(masked(v.sig .+ sigma .* n)), t2)
            cases["nan_cut0.99"] = (chain2(masked(v.sig .+ sigma .* n); cut=0.99), t2)
            post = s2full .+ sigma .* n2full
            live2 = d2 .> 0
            for k in 1:nv
                @views post[:, :, k][.!live2] .= NaN
            end
            cases["nan_after_smoothing"] = (post, t2)
            cases["zero_cut0.5"] = (chain2(zeroed(v.sig .+ sigma .* n)), t2)
            cr = chain2(masked(v.sig .+ sigma .* n .* ramp))
            cases["ramp_scalar"] = (cr, t2)
            # per-pixel threshold: rms map of the same chain applied to noise alone
            nr = chain2(masked(sigma .* n .* ramp))
            rmsmap = dropdims(std(replace(nr, NaN => 0.0); dims=3); dims=3)
            cases["ramp_pixel"] = (cr, repeat(rmsmap, 1, 1, nv))
            for (cname, (cube, thr)) in cases, em in (false, true)
                r = vgt(cube, thr; edge_mask=em)
                for (lo, hi) in DBINS
                    m = interior .& (d2 .> lo) .& (d2 .<= hi)
                    count(m) > 0 || continue
                    a = score(r.psi, v.tq2, v.tu2, m)
                    common = m .& isfinite.(r.psi) .& isfinite.(ref.psi)
                    emit(io, "edges", "n$(normal)_$(cname)_em$(Int(em))_d$(lo)-$(hi)", snr, v.name, seed,
                         (; AM=a.AM, cov=a.cov, dAM=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, m),
                          T=tcos(r.psi, theta_n, common), T_ref=tcos(ref.psi, theta_n, common),
                          T_true=tcos(0.5 .* atan.(v.tu2, v.tq2), theta_n, common)))
                end
            end
        end
    end
end

function exp_seams(io, v, seed, snrs)
    nx, ny, nv = size(v.sig)
    n = make_noise(MersenneTwister(1000seed + 9), nx, ny, nv; grain=3.0)
    s2full = chain2(v.sig)
    n2full = chain2(n)
    na = size(s2full, 1)
    xs = nx ÷ 2
    prof = dropdims(sum(v.sig; dims=(1, 2)); dims=(1, 2))
    kp = argmax(prof)
    cum = cumsum(prof) ./ sum(prof)
    # wing: the blue channels up to 15 % of the flux (a block whose band misses one line wing);
    # core: a quarter of the channels around the peak
    chsets = Dict("wing" => 1:findfirst(>=(0.15), cum), "core" => max(1, kp - nv ÷ 8):min(nv, kp + nv ÷ 8 - 1))
    for (csname, chs) in chsets
        emit(io, "seams", "flux_share_$(csname)", 0, v.name, seed, (; share=sum(prof[chs]) / sum(prof), nch=length(chs)))
    end
    for normal in (1, 2)
        theta_t = normal == 1 ? pi / 2 : 0.0   # seam tangent
        coord = [(normal == 1 ? i : j) for i in 1:nx, j in 1:ny]
        d2 = [abs((normal == 1 ? i : j) - na ÷ 2 - 0.5) for i in 1:na, j in 1:na]
        interior = evalmask(na)
        for snr in snrs
            sigma = v.peak / snr
            t2 = sigma * std(n2full)
            ref = vgt(s2full .+ sigma .* n2full, t2)
            for (csname, chs) in chsets
                for g in (0.5, 0.0), L in (0, 175)
                    gain = [L == 0 ? (c > xs ? g : 1.0) : 1 - (1 - g) * clamp((c - xs + L / 2) / L, 0, 1) for c in coord]
                    sig = copy(v.sig)
                    for k in chs
                        @views sig[:, :, k] .*= gain
                    end
                    r = vgt(chain2(sig) .+ sigma .* n2full, t2)
                    seam_rows(io, r, ref, v, seed, snr, "n$(normal)_$(csname)_g$(g)_L$(L)", d2, interior, theta_t)
                end
                c = v.sig .+ sigma .* n
                for k in chs
                    @views c[:, :, k][coord .> xs] .= NaN
                end
                cc = chain2(c)
                r = vgt(cc, t2)
                seam_rows(io, r, ref, v, seed, snr, "n$(normal)_$(csname)_nan", d2, interior, theta_t)
                # NaN in a few channels only, replaced by 0 after smoothing (partial spectra)
                anyfin = dropdims(any(isfinite, cc; dims=3); dims=3)
                cz = copy(cc)
                for k in 1:nv
                    @views cz[:, :, k][anyfin .& .!isfinite.(cz[:, :, k])] .= 0.0
                end
                r = vgt(cz, t2)
                seam_rows(io, r, ref, v, seed, snr, "n$(normal)_$(csname)_nan0", d2, interior, theta_t)
            end
        end
    end
end

function seam_rows(io, r, ref, v, seed, snr, cname, d2, interior, theta_t)
    for (lo, hi) in DBINS
        m = interior .& (d2 .> lo) .& (d2 .<= hi)
        count(m) > 0 || continue
        a = score(r.psi, v.tq2, v.tu2, m)
        common = m .& isfinite.(r.psi) .& isfinite.(ref.psi)
        emit(io, "seams", "$(cname)_d$(lo)-$(hi)", snr, v.name, seed,
             (; AM=a.AM, cov=a.cov, dAM=dam_pair(r.psi, ref.psi, v.tq2, v.tu2, m),
              T=tcos(r.psi, theta_t, common), T_ref=tcos(ref.psi, theta_t, common)))
    end
end

const EXPS = Dict(
    "noise" => (exp_noise, (2.0, 4.0, 8.0, 16.0, 64.0), VIEWS_SUB),
    "split" => (exp_split, (4.0, 16.0, 64.0), VIEWS_SUB),
    "chains" => (exp_chains, (4.0, 16.0, 64.0), vcat(VIEWS_SUB, VIEWS_STR)),
    "thr" => (exp_thr, (4.0, 16.0), VIEWS_SUB),
    "weights" => (exp_weights, (4.0, 8.0, 16.0, 64.0), vcat(VIEWS_SUB, VIEWS_STR)),
    "support" => (exp_support, (2.0, 4.0), VIEWS_SUB),
    "dn" => (exp_dn, (4.0, 16.0), VIEWS_SUB),
    "edges" => (exp_edges, (8.0, 64.0), VIEWS_SUB),
    "seams" => (exp_seams, (8.0, 64.0), VIEWS_SUB),
    "wsplit" => (exp_wsplit, (4.0, 16.0, 64.0), vcat(VIEWS_SUB, VIEWS_STR)),
)

function main(args)
    o = parse_args(args)
    f, snrs, views = EXPS[o["exp"]]
    seeds = parse(Int, o["seeds"])
    nv = parse(Int, o["nv"])
    if o["quick"] == "1"
        views, seeds, snrs, nv = views[1:1], 1, snrs[end:end], 32
    end
    mkpath(o["out"])
    path = joinpath(o["out"], "$(o["exp"]).csv")
    open(path, "w") do io
        println(io, "exp,cond,snr,view,seed,metric,value")
        for (run, los) in views
            t0 = time()
            v = load_view(o["cats-dir"], run, los; nv, sth=parse(Float64, o["sigma-th"]))
            emit(io, o["exp"], "noiseless", Inf, v.name, 0, score_nt(vgt(chain2(v.sig), 1e-6 * v.peak), v, evalmask(size(v.tq2, 1))))
            for seed in 1:seeds
                f(io, v, seed, snrs)
            end
            @printf("%s %s done in %.0f s\n", o["exp"], v.name, time() - t0)
        end
    end
    println("wrote ", path)
end

main(ARGS)
