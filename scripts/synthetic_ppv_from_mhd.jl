using DrWatson
quickactivate(@__DIR__, "VGT")

using FFTW
using FITSIO
using Random
using Statistics

# Synthetic optically thin 21 cm cube from an MHD snapshot (CATS, Burkhart et al. 2020:
# dens, velx/y/z, magx/y/z FITS cubes in code units), plus the true plane-of-sky field
# direction, so VGT variants (cube as is, VDA velocity / density parts) can be scored
# against a known field instead of starlight.
#
#   julia --project=. scripts/synthetic_ppv_from_mhd.jl --dir data/external/cats/b1p.032 \
#       --snap t800 --los 3 --out validation_output/vda_sim/b1p.032 [--nv 64] [--sigma-th 0]
#
# Brightness: T(x, y, v) = sum over the line of sight of rho * G(v - v_los; sigma_th),
# with G a unit-area Gaussian per cell (thermal broadening sigma_th, code units; 0 =
# one channel). The velocity axis spans +-(4 std(v_los) + 4 sigma_th) in --nv channels.
# True angle: psi_true = 0.5 atan2(U, Q), Q = sum rho (B1^2 - B2^2), U = sum rho 2 B1 B2,
# in the same array-axis convention as the VGT psi output (radians; axis 1 then axis 2
# of the map). Writes ppv.fits (CDELT3 in m/s for the pipeline) and psi_true.fits.
#
# To match observed data: --beam-px F smooths every channel (and Q, U of the true angle)
# with a Gaussian beam of FWHM F pixels (periodic, as the simulation box is); --snr S adds
# white noise so that the median over the map of the peak brightness is S times the noise
# (fixed seed, --seed). The noise level is printed, for the VGT noise threshold.
# --noise-beam 1 smooths the noise with the same beam and rescales it to the same per-pixel
# level, as in a real radio map, where noise is correlated over the beam (default 0: white,
# which hurts gradient methods more than real data does).

function parse_args(args)
    opts = Dict{String,String}("los" => "3", "nv" => "64", "sigma-th" => "0", "snap" => "t800", "beam-px" => "0",
                              "snr" => "0", "seed" => "1", "noise-beam" => "0")
    i = 1
    while i <= length(args)
        startswith(args[i], "--") && i < length(args) || throw(ArgumentError("Bad argument: $(args[i])"))
        opts[args[i][3:end]] = args[i + 1]
        i += 2
    end
    for key in ("dir", "out")
        haskey(opts, key) || throw(ArgumentError("--$key is required"))
    end
    return opts
end

# Periodic Gaussian smoothing of every 2D slice along dims 1, 2 (FWHM in pixels).
function beam!(a, fwhm)
    fwhm > 0 || return a
    s = fwhm / sqrt(8log(2))
    nx, ny = size(a, 1), size(a, 2)
    kx, ky = fftfreq(nx), fftfreq(ny)
    g = [exp(-2pi^2 * s^2 * (x^2 + y^2)) for x in kx, y in ky]
    for k in axes(a, 3)
        v = view(a, :, :, k)
        v .= real.(ifft(fft(v) .* g))
    end
    return a
end

# Optically thin PPV cube: every cell adds rho * G(v - v_los; sth) * dv. Kept in its own
# function so the threaded loop's locals cannot be captured from (and shared with) main.
function project_ppv(rho, v, vel, sth)
    nx, ny, nz = size(rho)
    nv = length(vel)
    dv = step(vel)
    vmin = first(vel)
    h = ceil(Int, 4sth / dv)
    ppv = zeros(Float32, nx, ny, nv)
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

rd(dir, name, snap) = Float64.(FITS(f -> read(f[1]), joinpath(dir, "$(name)_$(snap).fits")))

function main(args)
    opts = parse_args(args)
    dir, snap = opts["dir"], opts["snap"]
    los = parse(Int, opts["los"])
    sky = Tuple(setdiff(1:3, los))                       # map axes, in order
    perm = (sky..., los)
    rho = permutedims(rd(dir, "dens", snap), perm)
    v = permutedims(rd(dir, ("velx", "vely", "velz")[los], snap), perm)
    b1 = permutedims(rd(dir, ("magx", "magy", "magz")[sky[1]], snap), perm)
    b2 = permutedims(rd(dir, ("magx", "magy", "magz")[sky[2]], snap), perm)
    nv = parse(Int, opts["nv"])
    sth_in = parse(Float64, opts["sigma-th"])
    vmax = 4 * std(v) + 4 * sth_in                       # room for the thermal wings, or the edge clamp drops flux
    vel = range(-vmax, vmax; length=nv)
    dv = step(vel)
    sth = max(sth_in, dv / 2)                            # at least half a channel, so every cell lands
    ppv = project_ppv(rho, v, vel, sth)
    fwhm = parse(Float64, opts["beam-px"])
    beam!(ppv, fwhm)
    sigma_n = 0.0
    snr = parse(Float64, opts["snr"])
    if snr > 0
        sigma_n = median(maximum(ppv; dims=3)) / snr
        noise = randn(MersenneTwister(parse(Int, opts["seed"])), Float32, size(ppv))
        if opts["noise-beam"] == "1"
            beam!(noise, fwhm)
            noise ./= std(noise)
        end
        ppv .+= Float32(sigma_n) .* noise
    end
    q = sum(rho .* (b1 .^ 2 .- b2 .^ 2); dims=3)
    u = sum(rho .* 2 .* b1 .* b2; dims=3)
    beam!(q, fwhm)
    beam!(u, fwhm)
    psi = 0.5 .* atan.(dropdims(u; dims=3), dropdims(q; dims=3))
    out = opts["out"]
    mkpath(out)
    for (name, arr, hdr) in (("ppv.fits", ppv, FITSHeader(["CDELT3", "CRVAL3", "CRPIX3"], [dv * 1000, vel[1] * 1000, 1.0], ["", "", ""])),
                             ("psi_true.fits", Float32.(psi), FITSHeader(String[], [], String[])))
        p = joinpath(out, name)
        isfile(p) && rm(p)
        FITS(o -> write(o, arr; header=hdr), p, "w")
    end
    println("map axes ", sky, ", line of sight ", los, "; mean B along map axes (", round(mean(b1); digits=3), ", ",
            round(mean(b2); digits=3), "), along los ", round(mean(permutedims(rd(dir, ("magx", "magy", "magz")[los], snap), perm)); digits=3),
            "; v_los std ", round(std(v); digits=3), ", channel ", round(dv; digits=4), ", sigma_th ", round(sth; digits=4),
            "; beam ", fwhm, " px, noise ", sigma_n, "; wrote ", out)
end

# Run only as a script, so validate_faint_gas.jl can include the generator functions.
abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
