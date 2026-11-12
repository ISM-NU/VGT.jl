using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Printf
using Statistics

# Score VGT angle maps against a known field direction (psi_true.fits from
# synthetic_ppv_from_mhd.jl). Prints the alignment measure AM = <cos 2(psi - psi_true)>
# (1 = aligned, 0 = random, -1 = perpendicular) and the median |difference| over
# pixels where both maps are finite.
#
#   julia --project=. scripts/score_vgt_vs_true.jl --true psi_true.fits psi_a.fits psi_b.fits ...

function main(args)
    length(args) >= 3 && args[1] == "--true" || throw(ArgumentError("usage: --true psi_true.fits psi1.fits ..."))
    t = Float64.(FITS(f -> read(f[1]), args[2]))
    for p in args[3:end]
        s = Float64.(FITS(f -> read(f[1]), p))
        size(s) == size(t) || throw(ArgumentError("$p is $(size(s)), true map is $(size(t))"))
        d = mod.(s .- t .+ pi / 2, pi) .- pi / 2
        ok = isfinite.(d)
        @printf("%-70s AM %.3f  median |dpsi| %5.1f deg  pixels %d\n", p, mean(cos.(2 .* d[ok])),
                rad2deg(median(abs.(d[ok]))), count(ok))
    end
end

main(ARGS)
