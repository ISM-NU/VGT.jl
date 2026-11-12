using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Statistics

# Velocity Decomposition Algorithm (Yuen, Ho & Lazarian 2021; eq. B33 of Lazarian et al.
# 2024, arXiv:2406.01545): split every channel map p into the part that follows the
# integrated intensity I (density fluctuations) and the rest, the velocity caustics
#   p_v = p - (<pI> - <p><I>) (I - <I>) / sigma_I^2,   p_d = p - p_v.
# VGT traces the field through velocity caustics, so running it on p_v removes the
# density structures that can lie perpendicular to the field.
#
#   julia --project=. scripts/vda_velocity_caustics.jl --input cube.fits --out pv.fits \
#       [--crop x0,y0,size] [--block 0] [--part velocity]
#
# --part velocity (p_v, default), density (p_d), or none (the crop unchanged, as the
# baseline for comparing VGT with and without the decomposition on the same pixels).
# --block 0 takes the averages over the whole map, as in the papers. --block B uses
# B x B pixel blocks instead (the regression coefficient is computed per block and
# applied to that block), for maps where large-scale gradients such as galactic
# rotation would otherwise dominate the p-I correlation. Use the raw (unclipped) cube;
# NaN pixels stay NaN and are left out of the averages.

function parse_args(args)
    opts = Dict{String,String}("crop" => "", "block" => "0", "part" => "velocity")
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

# p_v (or p_d with density = true) for one block of pixels (all channels), in place.
function caustics!(c, xs, ys, density)
    nv = size(c, 3)
    I = dropdims(sum(view(c, xs, ys, :); dims=3); dims=3)    # NaN wherever any channel is NaN
    ok = isfinite.(I)
    Iok = I[ok]
    vI = length(Iok) < 16 ? 0.0 : var(Iok; corrected=false)
    if !(vI > 0)                     # no regression possible: blank the block rather than pass p through
        view(c, xs, ys, :) .= NaN32
        return
    end
    mI = mean(Iok)
    dI = I .- mI
    for k in 1:nv
        p = view(c, xs, ys, k)
        pk = p[ok]
        beta = (mean(pk .* Iok) - mean(pk) * mI) / vI
        if density
            @. p = beta * dI
        else
            @. p = p - beta * dI
        end
    end
end

function main(args)
    opts = parse_args(args)
    cube, header = FITS(opts["input"]) do f
        h = read_header(f[1])
        if isempty(opts["crop"])
            Float32.(read(f[1])), h
        else
            x0, y0, sz = parse.(Int, split(opts["crop"], ","))
            for (key, off) in (("CRPIX1", x0), ("CRPIX2", y0))
                haskey(h, key) && (h[key] = Float64(h[key]) - (off - 1))
            end
            Float32.(read(f[1], x0:(x0 + sz - 1), y0:(y0 + sz - 1), :)), h
        end
    end
    nx, ny, _ = size(cube)
    b = parse(Int, opts["block"])
    part = opts["part"]
    part in ("velocity", "density", "none") || throw(ArgumentError("--part must be velocity, density or none"))
    density = part == "density"
    if part == "none"
    elseif b == 0
        caustics!(cube, 1:nx, 1:ny, density)
    else
        blocks = [(bx:min(bx + b - 1, nx), by:min(by + b - 1, ny)) for by in 1:b:ny for bx in 1:b:nx]
        Threads.@threads for (xs, ys) in blocks
            caustics!(cube, xs, ys, density)
        end
    end
    out = opts["out"]
    mkpath(dirname(out))
    isfile(out) && rm(out)
    FITS(o -> write(o, cube; header=header), out, "w")
    println("wrote ", out, " (", part, ", block ", b, ")")
end

main(ARGS)
