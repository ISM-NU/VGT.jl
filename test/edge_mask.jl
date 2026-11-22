using VGT
using FITSIO
using Test

# Smooth deterministic cube with enough texture for every window to fit.
function edge_mask_cube(nx=60, ny=44, nv=3)
    cube = zeros(Float64, nx, ny, nv)
    for k in 1:nv, j in 1:ny, i in 1:nx
        cube[i, j, k] = 20 + 0.3 * (cos(0.5) * i + sin(0.5) * j) +
                        0.4 * sin(1.7i + 2.3j + 3k) + 0.3 * cos(0.9i * j / 7 + k)
    end
    return cube
end

# Offsets of the pixels one gradient reads, written out independently of the package.
const EDGE_MASK_STENCILS = Dict(
    :central_checked => [(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)],
    :sobel_checked => [(a, b) for a in -1:1 for b in -1:1],
    :isotropic_checked => [(a, b) for a in -1:1 for b in -1:1],
    :nan_edge_backward => [(0, 0), (-1, 0), (0, -1)],
)

# Brute force: does the footprint of the window at (i, j) reach a blank or leave the array?
function footprint_touches_blank(cube, i, j, dn, mode)
    nx, ny, _ = size(cube)
    half = div(dn, 2)
    for b in (j - half + 1):(j + half), a in (i - half + 1):(i + half), (di, dj) in EDGE_MASK_STENCILS[mode]
        x, y = a + di, b + dj
        (1 <= x <= nx && 1 <= y <= ny) || return true
        all(isfinite, @view cube[x, y, :]) || return true
    end
    return false
end

@testset "edge_mask is off by default and leaves output unchanged" begin
    cube = edge_mask_cube()
    # Values from the kernel before the option existed (origin/main): Q and U at
    # one pixel and their sums. rtol allows for FMA differences between CPUs.
    reference = Dict(
        :central_checked => (41.6773132782186, 83.04135271219344, 81183.65312313646, 175404.24625816356),
        :sobel_checked => (58.75134899016493, 72.32687429569425, 112451.04824709806, 163303.23586601944),
        :isotropic_checked => (64.66746903706935, 67.04935510351345, 114693.53319797202, 159659.61769228347),
        :nan_edge_backward => (0.12821949704532457, 92.98738493210904, -14239.781494888335, 196255.4980367299),
    )
    for (mode, (q_ref, u_ref, q_sum, u_sum)) in reference
        q, u = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_mode=mode)
        @test q[30, 22] ≈ q_ref rtol = 1e-10
        @test u[30, 22] ≈ u_ref rtol = 1e-10
        @test sum(q) ≈ q_sum rtol = 1e-10
        @test sum(u) ≈ u_sum rtol = 1e-10
        q0, u0 = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_mode=mode, edge_mask=false)
        @test q0 == q && u0 == u
    end
    @test !VGTConfig(input_file="x.fits", dn=6, ker=1, noise=1.0, min_valid=4).edge_mask
    @test !haskey(VGT.run_parameters(VGTConfig(input_file="x.fits", dn=6, ker=1, noise=1.0, min_valid=4)), "edge_mask")
end

@testset "edge_mask masks exactly the pixels whose footprint reaches a blank or the border" begin
    dn = 6
    half = dn ÷ 2
    cube = edge_mask_cube()
    nx, ny, _ = size(cube)
    # Blank block x = 25:34, y = 15:24 in every channel, plus one blank voxel in one channel.
    cube[25:34, 15:24, :] .= NaN
    cube[50, 38, 2] = Inf
    for mode in (:central_checked, :sobel_checked, :isotropic_checked, :nan_edge_backward)
        q0, u0 = VGT_QU_mv_2d(cube, dn, 1, 3, 1.0; min_valid=10, gradient_mode=mode)
        q, u = VGT_QU_mv_2d(cube, dn, 1, 3, 1.0; min_valid=10, gradient_mode=mode, edge_mask=true)
        expected = [footprint_touches_blank(cube, i, j, dn, mode) for i in 1:nx, j in 1:ny]
        masked = isnan.(q)
        @test masked == expected
        @test isnan.(u) == expected
        # Unmasked pixels are untouched by the option.
        @test q[.!masked] == q0[.!masked]
        @test u[.!masked] == u0[.!masked]
        @test count(.!masked) > 0
        # The option also removes partial-window fits that the default makes next to a blank.
        @test any(masked .& .!isnan.(q0) .& (q0 .!= 0))
    end

    # Mask width on a row through the blank block (y = 20), and at the array border.
    q, _ = VGT_QU_mv_2d(cube, dn, 1, 3, 1.0; min_valid=10, edge_mask=true)
    row = isnan.(q[:, 20])
    # central_checked footprint is i-half : i+half+1: half+1 pixels masked left of the block, half on the right.
    @test findlast(!, row[1:24]) == 24 - (half + 1)
    @test findfirst(!, row[35:end]) + 34 == 35 + half
    # Border: unmasked from half+1 to nx-half-1 on a blank-free map.
    clean = edge_mask_cube()
    qc, _ = VGT_QU_mv_2d(clean, dn, 1, 3, 1.0; min_valid=10, edge_mask=true)
    cols = findall(!isnan, qc[:, 22])
    @test first(cols) == half + 1 && last(cols) == nx - half - 1
    rows = findall(!isnan, qc[30, :])
    @test first(rows) == half + 1 && last(rows) == ny - half - 1
    # The default keeps the pixels at the outer ring of the computed range.
    qd, _ = VGT_QU_mv_2d(clean, dn, 1, 3, 1.0; min_valid=10)
    @test qd[half, 22] != 0 && qd[nx - half, 22] != 0
end

@testset "edge_mask ignores blank channels and checks weight_cube" begin
    cube = edge_mask_cube()
    cube[:, :, 3] .= NaN
    q, _ = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, edge_mask=true)
    @test count(!isnan, q) > 0
    cube = edge_mask_cube()
    weight = copy(cube)
    weight[20:22, 20:22, 1] .= NaN
    q, _ = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, weight_cube=weight, edge_mask=true)
    @test isnan(q[21, 21]) && isnan(q[21 + 3, 21]) && !isnan(q[40, 21])
end

@testset "edge_mask through config, TOML and the pipeline" begin
    keys = ["CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3", "CRPIX1", "CRPIX2", "CRPIX3",
            "CDELT1", "CDELT2", "CDELT3", "CUNIT3"]
    h = FITSHeader(keys, Any["RA---TAN", "DEC--TAN", "VRAD", 80.0, -69.0, 0.0, 12.0, 20.0, 1.0,
                             -0.5, 0.5, 1000.0, "m/s"], fill("", length(keys)))
    mktempdir() do dir
        cube = Float32.(edge_mask_cube(30, 40, 3))
        cube[1:6, :, :] .= NaN32          # survey edge written as blanks
        cube[20:23, 10:14, :] .= NaN32
        path = joinpath(dir, "c.fits")
        FITS(f -> write(f, cube; header=h), path, "w")
        base = (; input_file=path, keep_intermediate_chunks=false, noise=1.0, dn=6, ker=1, min_valid=4)

        off = run_pipeline(VGTConfig(; base..., ny_chunks=1, output_dir=joinpath(dir, "off")))
        on1 = run_pipeline(VGTConfig(; base..., edge_mask=true, ny_chunks=1, output_dir=joinpath(dir, "on1")))
        on3 = run_pipeline(VGTConfig(; base..., edge_mask=true, ny_chunks=3, output_dir=joinpath(dir, "on3")))
        # Chunked and one-piece runs agree with the mask on.
        for f in (:psi, :Q, :U, :I, :fit_support)
            @test isequal(getfield(on3, f), getfield(on1, f))
        end
        # Final psi is blank wherever the kernel masked, which smoothing would otherwise fill.
        q, _ = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=4, edge_mask=true)
        @test all(isnan, on1.psi[isnan.(q)])
        @test count(isfinite, on1.psi) > 0
        @test count(isfinite, on1.psi) < count(isfinite, off.psi)
        # No direction within the footprint radius of the blank edge: the first
        # finite column is the first with a clean window (x >= 7 + half = 10).
        @test findfirst(isfinite, on1.psi[:, 25]) >= 7 + 3

        toml = joinpath(dir, "config.toml")
        write(toml, """
        [input]
        input_file = "$(path)"
        output_dir = ""
        [noise]
        noise = 1.0
        [vgt]
        dn = 6
        ker = 1
        min_valid = 4
        edge_mask = true
        [execution]
        ny_chunks = 1
        """)
        loaded = load_config(toml)
        @test loaded.edge_mask
        @test VGT.run_parameters(loaded)["edge_mask"] == "true"
        @test VGT.run_parameter_hash(loaded) != VGT.run_parameter_hash(VGT.config_with(loaded; edge_mask=false))
    end
end
