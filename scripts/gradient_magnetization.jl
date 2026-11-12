using DrWatson
quickactivate(@__DIR__, "VGT")

using FFTW
using FITSIO
using Printf
using Random
using Statistics

# Magnetization from the anisotropy of velocity-centroid gradients (Lazarian, Yuen &
# Pogosyan 2024, ApJ 974, 237, eqs. 22 and 40). In every B x B block the gradient
# covariance s_ij = <g_i g_j> gives
#   J2 = ((s_xx - s_yy)^2 + 4 s_xy^2) / (s_xx + s_yy)^2,
# 1 for perfectly aligned gradients, 0 for isotropic ones. For a field in the plane of
# the sky the theory predicts J2 ~ exp(-2 M_A^2), so M_A,mes = sqrt(-ln(J2) / 2) (the
# measured value is M_A / sin(gamma) for a field inclined by gamma to the line of sight).
# With the velocity dispersion this gives the plane-of-sky field (their eqs. 66, 70):
#   B ~ sqrt(4 pi rho) dv / M_A      (magnetic driving, classical DCF form)
#   B ~ sqrt(4 pi rho) dv / sqrt(M_A) (velocity driving)
# so the script writes the inputs (J2, M_A,mes, dv) and leaves rho and the form to the user.
#
#   julia --project=. scripts/gradient_magnetization.jl --input cube.fits --out mag.fits \
#       [--block 32] [--noise 0] [--crop x0,y0,size] [--smooth 0] [--source centroid] [--sigma-th 0] \
#       [--null 0] [--highpass 0]
#
# Range (tested on the ten CATS Cho runs, M_A = dB/B0 as in the paper's eq. 6):
# - Strong and weak field only separate in transonic/subsonic gas (sonic Mach <~ 1.5); in
#   supersonic gas every run gives M_A,mes ~ 0.8-1, whatever the field.
# - Subsonic, dB/B0 = 0.5-0.6: M_A,mes 0.5-0.6 clean, 0.6-0.7 at a 3 px beam and peak S/N 23.
#   dB/B0 ~ 5: M_A,mes only ~1.1 (J2 floors at ~0.1 per 64 px block), so read it as "M_A >~ 1".
# - Noise drives J2 down (towards weak field): at S/N 23 use --noise (3 sigma threshold) and
#   --smooth at the beam (sigma = FWHM / 2.35), which restores the beam-only J2; at S/N 8 J2
#   is lost. Block 64 separates the fields better than 32.
# - These absolute values hold only for maps whose --null level is near the simulations'
#   (~0.05 at 64 px); the LMC's is 0.004 at 64 px (more eddies per block), so on other maps
#   compare J2 / null between runs, not M_A,mes.
# - DV is not usable on noisy cubes: the threshold clips the line wings (half the clean value
#   at S/N 23), so do not use it for a field strength yet.
#
# --source centroid takes the gradients of the velocity centroid (as in the paper);
# channels sums the gradient products of every channel map instead (the gradients VGT
# uses). --smooth s smooths each map with a Gaussian of sigma s pixels (FWHM 2.35 s) first.
# Centroid C = sum v T / sum T over channels with T > 3 * --noise (all channels if 0).
# dv per block: sqrt(var of C over the block + mean of the line's second moment - sigma_th^2),
# so thermal broadening (--sigma-th, same units as the velocity axis) does not count as turbulence.
# --highpass h (centroid only): subtract the centroid smoothed with a Gaussian of sigma h px
# first, so a large-scale velocity field (galactic rotation) does not count as anisotropy.
# --null N (centroid only): N isotropic surrogates of the centroid map (same azimuthally
# averaged power spectrum and mask, random phases), run through the same smoothing and
# blocks; their median J2 is printed as the level a map without any preferred direction
# gives at this block size (the J2 floor). Phase randomization alone would keep the
# anisotropy, so the amplitudes are averaged over direction first.
# Output: a FITS with HDUs J2, MA_MES, DV, CENTROID_STD (one value per block), CENTROID
# (the full-resolution map the gradients came from, before --smooth) and a
# printed summary (medians), for scoring against simulations with known M_A.

function parse_args(args)
    opts = Dict{String,String}("block" => "32", "noise" => "0", "crop" => "", "smooth" => "0", "source" => "centroid", "sigma-th" => "0", "null" => "0", "highpass" => "0")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("Bad argument: $(args[i])"))
        opts[args[i][3:end]] = args[i + 1]
        i += 2
    end
    for key in ("input", "out")
        haskey(opts, key) || throw(ArgumentError("--$key is required"))
    end
    return opts
end

function velocity_axis(h, nv)
    cdelt = Float64(h["CDELT3"])
    crval = Float64(h["CRVAL3"])
    crpix = haskey(h, "CRPIX3") ? Float64(h["CRPIX3"]) : 1.0
    return [(crval + (k - crpix) * cdelt) / 1000 for k in 1:nv]   # CDELT3 in m/s -> km/s (simulations: code units)
end

function moments(cube, vel, thr)
    nx, ny, nv = size(cube)
    c = fill(NaN, nx, ny)
    w2 = fill(NaN, nx, ny)
    Threads.@threads for y in 1:ny
        for x in 1:nx
            s0 = s1 = s2 = 0.0
            for k in 1:nv
                t = cube[x, y, k]
                (isfinite(t) && t > thr) || continue
                s0 += t
                s1 += t * vel[k]
                s2 += t * vel[k]^2
            end
            s0 > 0 || continue
            c[x, y] = s1 / s0
            w2[x, y] = max(s2 / s0 - (s1 / s0)^2, 0.0)
        end
    end
    return c, w2
end

# NaN-aware separable Gaussian smoothing (sigma in pixels; 0 = none).
function smooth2(a, s)
    s > 0 || return a
    r = ceil(Int, 3s)
    w = [exp(-0.5 * (i / s)^2) for i in -r:r]
    out = a
    for dim in 1:2
        src = out
        out = fill(NaN, size(a))
        n1, n2 = size(a)
        for y in 1:n2, x in 1:n1
            num = den = 0.0
            for (o, wt) in zip(-r:r, w)
                xx, yy = dim == 1 ? (x + o, y) : (x, y + o)
                (1 <= xx <= n1 && 1 <= yy <= n2) || continue
                v = src[xx, yy]
                isfinite(v) || continue
                num += wt * v
                den += wt
            end
            isfinite(a[x, y]) && den > 0 && (out[x, y] = num / den)
        end
    end
    return out
end

# Add one map's central-difference gradient products to the per-block sums
# (s_xx, s_yy, s_xy, count); pixels with a missing neighbour are skipped.
function block_cov!(acc, m, b)
    n1, n2 = size(m)
    for y in 2:(n2 - 1), x in 2:(n1 - 1)
        i, j = (x - 1) ÷ b + 1, (y - 1) ÷ b + 1
        (i <= size(acc, 1) && j <= size(acc, 2)) || continue
        gx = (m[x + 1, y] - m[x - 1, y]) / 2
        gy = (m[x, y + 1] - m[x, y - 1]) / 2
        (isfinite(gx) && isfinite(gy)) || continue
        acc[i, j, 1] += gx^2
        acc[i, j, 2] += gy^2
        acc[i, j, 3] += gx * gy
        acc[i, j, 4] += 1
    end
end

# Isotropic Gaussian surrogate of a map: |FFT| averaged over direction in rings of |k|, random
# phases, same NaN mask.
function isotropic_surrogate(m, rng)
    ok = isfinite.(m)
    a = zeros(size(m))
    a[ok] .= m[ok] .- mean(m[ok])
    amp = abs.(fft(a))
    n1, n2 = size(m)
    kr = [round(Int, hypot(k1 * n2, k2 * n1) / max(n1, n2)) for k1 in fftfreq(n1, n1), k2 in fftfreq(n2, n2)]
    ring = zeros(maximum(kr) + 1)
    cnt = zeros(maximum(kr) + 1)
    for i in eachindex(kr)
        ring[kr[i] + 1] += amp[i]
        cnt[kr[i] + 1] += 1
    end
    s = real.(ifft([ring[k + 1] / cnt[k + 1] for k in kr] .* cis.(2pi .* rand(rng, n1, n2))))
    s[.!ok] .= NaN
    return s
end

function block_j2(acc, b)
    j2 = fill(NaN, size(acc, 1), size(acc, 2))
    for j in axes(acc, 2), i in axes(acc, 1)
        sxx, syy, sxy, n = acc[i, j, :]
        n >= b^2 ÷ 2 && (j2[i, j] = ((sxx - syy)^2 + 4sxy^2) / (sxx + syy)^2)
    end
    return j2
end

function main(args)
    opts = parse_args(args)
    cube, h = FITS(opts["input"]) do f
        h = read_header(f[1])
        if isempty(opts["crop"])
            Float64.(read(f[1])), h
        else
            x0, y0, sz = parse.(Int, split(opts["crop"], ","))
            Float64.(read(f[1], x0:(x0 + sz - 1), y0:(y0 + sz - 1), :)), h
        end
    end
    nx, ny, nv = size(cube)
    vel = velocity_axis(h, nv)
    thr = 3 * parse(Float64, opts["noise"])
    c, w2 = moments(cube, vel, thr)
    sm = parse(Float64, opts["smooth"])
    b = parse(Int, opts["block"])
    sth = parse(Float64, opts["sigma-th"])
    mx, my = nx ÷ b, ny ÷ b
    acc = zeros(mx, my, 4)
    hp = parse(Float64, opts["highpass"])
    (hp > 0 && opts["source"] != "centroid") && throw(ArgumentError("--highpass needs --source centroid"))
    cg = hp > 0 ? c .- smooth2(c, hp) : c                # map the gradients are taken from
    if opts["source"] == "centroid"
        block_cov!(acc, smooth2(cg, sm), b)
    elseif opts["source"] == "channels"
        for k in 1:nv
            m = cube[:, :, k]
            thr > 0 && (m[m .<= thr] .= NaN)
            block_cov!(acc, smooth2(m, sm), b)
        end
    else
        throw(ArgumentError("--source must be centroid or channels"))
    end
    j2 = block_j2(acc, b)
    dv = fill(NaN, mx, my)
    cstd = fill(NaN, mx, my)
    for j in 1:my, i in 1:mx
        xs, ys = ((i - 1) * b + 1):(i * b), ((j - 1) * b + 1):(j * b)
        acc[i, j, 4] >= b^2 ÷ 2 || continue
        cb = filter(isfinite, vec(c[xs, ys]))
        wb = filter(isfinite, vec(w2[xs, ys]))
        cstd[i, j] = std(cb)
        dv[i, j] = sqrt(var(cb) + max(mean(wb) - sth^2, 0.0))
    end
    ma = sqrt.(max.(-log.(j2) ./ 2, 0.0))
    nnull = parse(Int, opts["null"])
    null_j2 = Float64[]
    if nnull > 0
        opts["source"] == "centroid" || throw(ArgumentError("--null needs --source centroid"))
        rng = MersenneTwister(1)
        for _ in 1:nnull
            nacc = zeros(mx, my, 4)
            block_cov!(nacc, smooth2(isotropic_surrogate(cg, rng), sm), b)
            push!(null_j2, median(filter(isfinite, vec(block_j2(nacc, b)))))
        end
    end
    out = opts["out"]
    mkpath(dirname(abspath(out)))
    isfile(out) && rm(out)
    FITS(out, "w") do f
        for (name, arr) in (("J2", j2), ("MA_MES", ma), ("DV", dv), ("CENTROID_STD", cstd), ("CENTROID", c))
            write(f, Float32.(arr); name=name)
        end
    end
    fin(a) = filter(isfinite, vec(a))
    @printf("%s %s smooth %s highpass %s block %d: blocks %d  median J2 %.3f  median M_A,mes %.2f  median dv %.3f  (centroid std %.3f)\n",
            opts["input"], opts["source"], opts["smooth"], opts["highpass"], b, length(fin(j2)), median(fin(j2)), median(fin(ma)), median(fin(dv)), median(fin(cstd)))
    isempty(null_j2) || @printf("  isotropic null (%d surrogates): median J2 %.3f (range %.3f-%.3f)\n",
                                nnull, median(null_j2), minimum(null_j2), maximum(null_j2))
end

main(ARGS)
