using Printf, Statistics
module MW end
for e in Meta.parseall(read("scripts/validate_mw_like_synthetic.jl", String)).args
    e isa LineNumberNode && continue
    s = string(e)
    (occursin("quickactivate", s) || startswith(s, "using DrWatson") || startswith(s, "if abspath")) && continue
    Core.eval(MW, e)
end
# Fraction of fitted channel-windows (>= min_valid central_checked gradients) whose
# centre intensity is <= 0 or below the noise threshold; their Stokes weight is the
# raw centre intensity (vgt.jl:381-395), so a negative centre flips the orientation.
for (label, b) in (("many_filaments", MW.build_many_filaments), ("velocity_components", MW.build_velocity_components),
                   ("sparse_co_like", MW.build_sparse_co_like), ("planck_misaligned", MW.build_planck_misaligned))
    cube, truth, noise = b(; nx=96, ny=96)
    nx, ny, nv = size(cube); dn = 12; h = 6
    nfit = 0; nneg = 0; nbelow = 0; wneg = 0.0; wpos = 0.0
    for k in 1:nv
        ok = falses(nx, ny)
        for j in 2:ny-1, i in 2:nx-1
            ok[i, j] = all(cube[a, c, k] >= noise for (a, c) in ((i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1)))
        end
        for jj in h:ny-h, ii in h:nx-h
            count(@view ok[ii-h+1:ii+h, jj-h+1:jj+h]) >= 50 || continue
            nfit += 1
            w = cube[ii, jj, k]
            w < noise && (nbelow += 1)
            w <= 0 && (nneg += 1; wneg += -w)
            w > 0 && (wpos += w)
        end
    end
    @printf("%-20s fitted windows %7d  centre<noise %.4f  centre<=0 %.5f  |neg weight|/pos weight %.2e\n",
            label, nfit, nbelow / nfit, nneg / nfit, wneg / wpos)
end
