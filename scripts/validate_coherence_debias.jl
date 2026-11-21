using DrWatson
quickactivate(@__DIR__, "VGT")

# Synthetic validation of the noise-debiased beam coherence (coherence_phasor,
# coherence_floor2, debias_coherence, cross_coherence). Synthetic data only.
#
# Question: the beam coherence R = |<exp(2i psi)>| of a VGT angle map is biased
# up by noise (pure noise gives a floor, not 0). Does subtracting the floor
# measured on a pure-noise run (R_db), or the cross product of two independent
# noise realizations (R_x), recover the coherence of the noiseless VGT map
# (R_true, same settings) better than the raw R?
#
# Signal: baseline 1000 (the intensity mask never acts) plus, in every channel
# (envelope 1 + 0.4 sin(2 pi k / nv)), a mixture a(x) ramp + (1 - a(x)) GRF:
# a(x) falls linearly from 1 at the left edge to 0 at the right, so the true
# coherence runs from 1 (one gradient direction) to the low value of a random
# field (power-law spectrum, slope -3). Both parts have the same rms gradient
# before mixing; the mixture is scaled so its median central-difference
# gradient magnitude is S/N * sigma_g.
# Noise: white noise convolved with a Gaussian (FWHM 3 px), isotropic (axis
# ratio 1) and anisotropic (1.05), as validate_gradient_whitening.jl.
# VGT as in production: dn = 20, min_valid = 100, central_checked, Q/U
# smoothed with ker = 4, direct-fit support. Coherence beam FWHM 30 px (a 10'
# beam on 20" pixels). Interior pixels only (beam and window edges cut).
# Pass rule (written before running): at S/N <= 1 the debiased estimator's
# median error |R_db - R_true| is at most half the raw |R - R_true| where
# R_true < 0.5, and at S/N >= 2 it does not add more than 0.02 of bias where
# R_true > 0.8.
#
# Environment: VGT_CD_SIZE (256), VGT_CD_NV (8), VGT_CD_SEEDS (3).
# Run (cluster): julia --project=. --threads=16 scripts/validate_coherence_debias.jl

import VGT
using FFTW
using Printf
using Random
using Statistics
using VGT: coherence_phasor, coherence_floor2, debias_coherence, cross_coherence

const GW = Module(:GW)
Base.include(GW, joinpath(@__DIR__, "validate_gradient_whitening.jl"))   # correlated_noise, vgt_psi, BASE, DN, KER

const SIZE = parse(Int, get(ENV, "VGT_CD_SIZE", "256"))
const NV = parse(Int, get(ENV, "VGT_CD_NV", "8"))
const NSEEDS = parse(Int, get(ENV, "VGT_CD_SEEDS", "3"))
const BEAM = 30.0
const LEVELS = (0.5, 1.0, 2.0, 4.0)

function grf(n, rng; slope=3.0)
    k = [hypot(min(i - 1, n - i + 1), min(j - 1, n - j + 1)) for i in 1:n, j in 1:n]
    amp = @. ifelse(k == 0, 0.0, k^(-slope / 2))
    f = amp .* cis.(2π .* rand(rng, n, n)) .* randn(rng, n, n)
    return real.(ifft(f))
end

function cdiff_mag(img)
    n1, n2 = size(img)
    return [hypot((img[i + 1, j] - img[i - 1, j]) / 2, (img[i, j + 1] - img[i, j - 1]) / 2) for i in 2:(n1 - 1), j in 2:(n2 - 1)]
end

function signal_image(rng)
    ic, jc = (SIZE + 1) / 2, (SIZE + 1) / 2
    θ = deg2rad(20.0)
    ramp = [cos(θ) * (i - ic) + sin(θ) * (j - jc) for i in 1:SIZE, j in 1:SIZE]
    g = grf(SIZE, rng)
    g .*= 1 / sqrt(mean(abs2, cdiff_mag(g)))     # unit rms gradient, as the unit-slope ramp
    a = [clamp(1 - (i - 1) / (SIZE - 1) * 1.25, 0.0, 1.0) for i in 1:SIZE, j in 1:SIZE]
    return a .* ramp .+ (1 .- a) .* g
end

cube_from(img, scale) = [GW.BASE + scale * img[i, j] * (1 + 0.4 * sin(2π * k / NV)) for i in 1:SIZE, j in 1:SIZE, k in 1:NV]

function interior()
    h = GW.DN ÷ 2 + round(Int, BEAM)
    m = falses(SIZE, SIZE)
    m[h:(SIZE - h), h:(SIZE - h)] .= true
    return m
end

med(v) = (x = filter(isfinite, v); isempty(x) ? NaN : median(x))

function main()
    mask = interior()
    println("### noise-debiased VGT coherence, size $SIZE, nv $NV, beam FWHM $(BEAM) px, $NSEEDS seeds")
    nind = 4π * (BEAM / (2 * sqrt(2 * log(2))))^2
    @printf("independent pixels per beam %.0f -> floor of independent angles f^2 = %.4f\n", nind, 1 / nind)
    println("errors are median (estimate - R_true) over interior pixels in bins of R_true; f^2 from a pure-noise run (same noise type, other seed)")
    println("| noise | S/N | f^2 | N_eff | bin R_true | n px | raw R - R_true | R_db - R_true | R_x - R_true |")
    println("|---|---|---|---|---|---|---|---|---|")
    for (nlabel, ratio) in (("isotropic", 1.0), ("anisotropic 1.05", GW.AXIS_RATIO))
        # sigma_g of this noise (rms of the two central-difference gradient components)
        ref = GW.correlated_noise(SIZE, SIZE, 2, MersenneTwister(900); ratio=ratio)
        σg = sqrt(mean(abs2, cdiff_mag(ref[:, :, 1])) / 2)
        rows = Dict{Tuple{Float64, String}, Vector{NTuple{4, Float64}}}()
        f2s = Float64[]
        for seed in 1:NSEEDS
            rng = MersenneTwister(100 + seed)
            img = signal_image(rng)
            med_g = median(cdiff_mag(img)[mask[2:end-1, 2:end-1]])
            psin = GW.vgt_psi(GW.BASE .+ GW.correlated_noise(SIZE, SIZE, NV, MersenneTwister(500 + seed); ratio=ratio))
            f2 = coherence_floor2(psin, BEAM; mask=mask)
            push!(f2s, f2)
            for snr in LEVELS
                scale = snr * σg / med_g
                psit = GW.vgt_psi(cube_from(img, scale))
                na = GW.correlated_noise(SIZE, SIZE, NV, MersenneTwister(1000 * seed + 1); ratio=ratio)
                nb = GW.correlated_noise(SIZE, SIZE, NV, MersenneTwister(1000 * seed + 2); ratio=ratio)
                c0 = cube_from(img, scale)
                za = coherence_phasor(GW.vgt_psi(c0 .+ na), BEAM)
                zb = coherence_phasor(GW.vgt_psi(c0 .+ nb), BEAM)
                rt = abs.(coherence_phasor(psit, BEAM))
                rr = abs.(za)
                rd = debias_coherence(rr, f2)
                rx = cross_coherence(za, zb)
                for (bl, lo, hi) in (("< 0.5", 0.0, 0.5), ("0.5-0.8", 0.5, 0.8), ("> 0.8", 0.8, 1.01))
                    sel = mask .& isfinite.(rt) .& (rt .>= lo) .& (rt .< hi) .& isfinite.(rr)
                    count(sel) < 50 && continue
                    push!(get!(rows, (snr, bl), NTuple{4, Float64}[]),
                          (count(sel), med(rr[sel] .- rt[sel]), med(rd[sel] .- rt[sel]), med(rx[sel] .- rt[sel])))
                end
            end
        end
        f2 = mean(f2s)
        for snr in LEVELS, bl in ("< 0.5", "0.5-0.8", "> 0.8")
            haskey(rows, (snr, bl)) || continue
            v = rows[(snr, bl)]
            se(k) = length(v) > 1 ? std([x[k] for x in v]) / sqrt(length(v)) : NaN
            m(k) = mean([x[k] for x in v])
            @printf("| %s | %.1f | %.3f | %.1f | %s | %d | %+.3f ± %.3f | %+.3f ± %.3f | %+.3f ± %.3f |\n",
                    nlabel, snr, f2, 1 / f2, bl, round(Int, sum(x[1] for x in v)), m(2), se(2), m(3), se(3), m(4), se(4))
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
