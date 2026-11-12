using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Printf

# Smooth a cube from one Gaussian beam to a larger one and resample it to a new
# pixel size, keeping the same sky window. Used to put the native-resolution
# Milky Way foreground (about 30" beam, 12" pixels) on the sampling of the LMC
# cube (60" beam, 20" pixels) so ROHSA results can be compared at equal angular
# resolution and equal pixels per beam.
#
#   julia --project=. scripts/smooth_cube_to_beam.jl --input crop.fits --output smooth.fits \
#       --fwhm-in-arcsec 30 --fwhm-out-arcsec 60 [--pixel-out-arcsec 20] [--crop x0,y0,size]
#
# Use --crop on the raw input cube rather than a ROHSA region FITS (those are
# clipped to >= 0).
# NaN pixels are excluded with normalised convolution. Without --pixel-out-arcsec
# the original pixel grid is kept.

function parse_args(args)
    opts = Dict{String,Any}(
        "input" => nothing,
        "output" => nothing,
        "fwhm-in-arcsec" => nothing,
        "fwhm-out-arcsec" => nothing,
        "pixel-out-arcsec" => "",
        "crop" => "",
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("Unexpected positional argument: $arg"))
        key = arg[3:end]
        haskey(opts, key) || throw(ArgumentError("Unknown option --$key"))
        i == length(args) && throw(ArgumentError("Missing value for --$key"))
        opts[key] = args[i + 1]
        i += 2
    end
    for key in ("input", "output", "fwhm-in-arcsec", "fwhm-out-arcsec")
        opts[key] === nothing && throw(ArgumentError("--$key is required"))
    end
    return opts
end

function gaussian_kernel(sigma_px::Float64)
    half = ceil(Int, 4sigma_px)
    k = [exp(-(i^2) / (2sigma_px^2)) for i in -half:half]
    return k ./ sum(k), half
end

# Separable normalised convolution of one channel; `w` is 1 where data are
# finite, 0 elsewhere.
function smooth_plane(plane::AbstractMatrix{Float64}, w::AbstractMatrix{Float64}, kernel, half)
    nx, ny = size(plane)
    num = zeros(nx, ny)
    den = zeros(nx, ny)
    tmpn = zeros(nx, ny)
    tmpd = zeros(nx, ny)
    @inbounds for y in 1:ny, x in 1:nx
        sn = 0.0
        sd = 0.0
        for d in -half:half
            xx = x + d
            (1 <= xx <= nx) || continue
            kv = kernel[d + half + 1]
            sn += kv * plane[xx, y]
            sd += kv * w[xx, y]
        end
        tmpn[x, y] = sn
        tmpd[x, y] = sd
    end
    @inbounds for y in 1:ny, x in 1:nx
        sn = 0.0
        sd = 0.0
        for d in -half:half
            yy = y + d
            (1 <= yy <= ny) || continue
            kv = kernel[d + half + 1]
            sn += kv * tmpn[x, yy]
            sd += kv * tmpd[x, yy]
        end
        num[x, y] = sn
        den[x, y] = sd
    end
    return num, den
end

function bilinear(plane, fx::Float64, fy::Float64)
    nx, ny = size(plane)
    x0 = clamp(floor(Int, fx), 1, nx - 1)
    y0 = clamp(floor(Int, fy), 1, ny - 1)
    tx = clamp(fx - x0, 0.0, 1.0)
    ty = clamp(fy - y0, 0.0, 1.0)
    return (1 - tx) * (1 - ty) * plane[x0, y0] + tx * (1 - ty) * plane[x0 + 1, y0] +
           (1 - tx) * ty * plane[x0, y0 + 1] + tx * ty * plane[x0 + 1, y0 + 1]
end

function main(args)
    opts = parse_args(args)
    fwhm_in = parse(Float64, opts["fwhm-in-arcsec"])
    fwhm_out = parse(Float64, opts["fwhm-out-arcsec"])
    fwhm_out > fwhm_in || throw(ArgumentError("--fwhm-out-arcsec must exceed --fwhm-in-arcsec"))
    cube, header = FITS(opts["input"]) do f
        hdr = read_header(f[1])
        if isempty(opts["crop"])
            Float64.(read(f[1])), hdr
        else
            x0, y0, sz = parse.(Int, split(opts["crop"], ","))
            for (key, off) in (("CRPIX1", x0), ("CRPIX2", y0))
                haskey(hdr, key) && (hdr[key] = Float64(hdr[key]) - (off - 1))
            end
            Float64.(read(f[1], x0:(x0 + sz - 1), y0:(y0 + sz - 1), :)), hdr
        end
    end
    nx, ny, nv = size(cube)
    pix_in = abs(Float64(header["CDELT2"])) * 3600
    sigma_px = sqrt(fwhm_out^2 - fwhm_in^2) / (2sqrt(2log(2))) / pix_in
    kernel, half = gaussian_kernel(sigma_px)
    pix_out = isempty(opts["pixel-out-arcsec"]) ? pix_in : parse(Float64, opts["pixel-out-arcsec"])
    scale = pix_out / pix_in
    nxo, nyo = floor(Int, nx / scale), floor(Int, ny / scale)
    out = fill(Float32(NaN), nxo, nyo, nv)
    valid = Float64.(isfinite.(cube[:, :, 1]))
    for k in 1:nv
        plane = map(v -> isfinite(v) ? v : 0.0, @view cube[:, :, k])
        num, den = smooth_plane(plane, valid, kernel, half)
        sm = map((n, d) -> d > 0.5 ? n / d : NaN, num, den)
        for yo in 1:nyo, xo in 1:nxo
            # Output pixel centre in input pixel coordinates (same window edges).
            fx = (xo - 0.5) * scale + 0.5
            fy = (yo - 0.5) * scale + 0.5
            out[xo, yo, k] = Float32(bilinear(sm, fx, fy))
        end
    end
    for (ax, n) in ((1, nx), (2, ny))
        key = "CRPIX$ax"
        haskey(header, key) && (header[key] = (Float64(header[key]) - 0.5) / scale + 0.5)
        dkey = "CDELT$ax"
        haskey(header, dkey) && (header[dkey] = Float64(header[dkey]) * scale)
    end
    header["BMAJ"] = fwhm_out / 3600
    header["BMIN"] = fwhm_out / 3600
    mkpath(dirname(abspath(opts["output"])))
    FITS(opts["output"], "w") do f
        write(f, out; header=header)
    end
    @printf("Smoothed %d x %d x %d from %.0f\" to %.0f\" (kernel sigma %.2f px) and resampled to %d x %d at %.1f\" pixels\n",
            nx, ny, nv, fwhm_in, fwhm_out, sigma_px, nxo, nyo, pix_out)
end

main(ARGS)
