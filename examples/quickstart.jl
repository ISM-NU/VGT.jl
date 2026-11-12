# Quick start: run VGT end to end on a small synthetic cube with a known field.
# Background of the map panels: the velocity centroid, which VGT traces.
#
#   julia --project=. examples/quickstart.jl
#
# The cube's velocity centroid varies smoothly across the map, with a curved
# gradient, so the true gradient direction is known at every pixel; the
# magnetic-field proxy is that direction rotated by 90 degrees. The script
# writes the cube to a temporary FITS file, runs the full pipeline (the same
# code path as main.jl), compares the recovered angle with the truth and saves
# a figure to docs/images/quickstart.png. It takes about a minute.

using VGT
using CairoMakie
using FITSIO
using Random
using Statistics

const NX, NY, NV = 96, 96, 40
const ANGLE = deg2rad(37.0)    # gradient direction at the map centre
const SLOPE = 0.15             # centroid shift per pixel, channels
const CURVATURE = 1.2e-3       # bends the gradient across the map

# Gradient direction of the centroid field c(x, y) = SLOPE (x cos a + y sin a) + CURVATURE x y.
function true_gradient_angle(i, j)
    x, y = i - (NX + 1) / 2, j - (NY + 1) / 2
    return atan(SLOPE * sin(ANGLE) + CURVATURE * x, SLOPE * cos(ANGLE) + CURVATURE * y)
end

function synthetic_cube(; noise=0.05, rng=MersenneTwister(1))
    cube = Array{Float32}(undef, NX, NY, NV)
    v0 = (NV + 1) / 2
    for j in 1:NY, i in 1:NX
        x, y = i - (NX + 1) / 2, j - (NY + 1) / 2
        centroid = v0 + SLOPE * (x * cos(ANGLE) + y * sin(ANGLE)) + CURVATURE * x * y
        amplitude = 20.0  # constant, so every gradient in the cube comes from velocity
        for k in 1:NV
            cube[i, j, k] = amplitude * exp(-0.5 * ((k - centroid) / 6.0)^2) + noise * randn(rng)
        end
    end
    return cube
end

headless(d) = mod(d + pi / 2, pi) - pi / 2

function main()
    dir = mktempdir()
    input = joinpath(dir, "synthetic_cube.fits")
    cube = synthetic_cube()
    FITS(input, "w") do f
        write(f, cube)
    end

    cfg = VGTConfig(
        input_file=input, output_dir=joinpath(dir, "out"), keep_intermediate_chunks=false,
        noise=1.0, dn=10, ker=2, ny_chunks=1, min_valid=50, fitter_mode=:axial_vonmises,
    )
    run = run_pipeline(cfg)
    outputs = write_pipeline_outputs(run)

    truth = [true_gradient_angle(i, j) + pi / 2 for i in 1:NX, j in 1:NY]
    err = rad2deg.(abs.(headless.(run.psi .- truth)))
    inner = err[cfg.dn+1:NX-cfg.dn, cfg.dn+1:NY-cfg.dn]   # skip the window edge
    ok = filter(isfinite, inner)
    println("psi map: ", outputs.psi_path)
    println("median |psi - truth| away from the edges: ", round(median(ok); digits=2), " deg")

    save_figure(cube, truth, run.psi, ok)
end

function segments(angle; step=8, len=3.2)
    pts = Point2f[]
    for j in step÷2:step:NY, i in step÷2:step:NX
        a = angle[i, j]
        isfinite(a) || continue
        d = Point2f(len * cos(a), len * sin(a))
        push!(pts, Point2f(i, j) - d, Point2f(i, j) + d)
    end
    return pts
end

function save_figure(cube, truth, psi, errors)
    blue, orange = "#2a78d6", "#eb6834"
    ink, muted = "#2b2b29", "#6b6b66"
    m0 = dropdims(sum(cube; dims=3); dims=3)
    centroid = dropdims(sum(cube .* reshape(1:NV, 1, 1, NV); dims=3); dims=3) ./ m0
    fig = Figure(size=(1050, 360), fontsize=13, backgroundcolor=:white)
    for (col, (angle, colour, title)) in enumerate([(truth, blue, "True field"), (psi, orange, "VGT field (psi)")])
        ax = Axis(fig[1, col]; title=title, titlecolor=ink, aspect=DataAspect(),
                  xlabel="x (pixel)", ylabel="y (pixel)", xlabelcolor=muted, ylabelcolor=muted,
                  xgridvisible=false, ygridvisible=false)
        heatmap!(ax, centroid; colormap=[:white, "#c3c2b7"])  # velocity centroid (channels)
        linesegments!(ax, segments(angle); color=colour, linewidth=2)
        limits!(ax, 0.5, NX + 0.5, 0.5, NY + 0.5)
    end
    ax = Axis(fig[1, 3]; title="Recovered minus true angle", titlecolor=ink,
              xlabel="|psi - truth| (deg)", ylabel="pixels", xlabelcolor=muted, ylabelcolor=muted,
              xgridvisible=false, ygridcolor=(:black, 0.06))
    hist!(ax, errors; bins=range(0, 1.1 * quantile(errors, 0.999); length=30), color=blue,
          strokecolor=:white, strokewidth=1)
    out = joinpath(dirname(@__DIR__), "docs", "images", "quickstart.png")
    mkpath(dirname(out))
    save(out, fig; px_per_unit=2)
    println("figure: ", out)
end

main()
