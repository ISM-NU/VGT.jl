using VGT
using FITSIO
using Random
using Test

# velocity_range_kms must give exactly what a run on the cube cut to those channels gives,
# including estimated noise (which reads the cube from disk) and chunk padding.
@testset "velocity_range_kms: same products as a cube cut to the window" begin
    keys = ["CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3", "CRPIX1", "CRPIX2", "CRPIX3",
            "CDELT1", "CDELT2", "CDELT3", "CUNIT3"]
    hdr(crval3, crpix3, cdelt3) = FITSHeader(keys, Any["RA---TAN", "DEC--TAN", "VRAD", 80.0, -69.0, crval3,
                                                     12.0, 10.0, crpix3, -0.5, 0.5, cdelt3, "m/s"],
                                             fill("", length(keys)))

    # Channels at -2, -1, ..., 5 km/s.
    h = hdr(-2000.0, 1.0, 1000.0)
    cfg(r) = VGTConfig(input_file="c.fits", noise=1.0, dn=4, ker=1, min_valid=4, velocity_range_kms=r)
    @test selected_channels(cfg(Float64[]), h, 8) == 1:8
    @test selected_channels(cfg([0.5, 3.0]), h, 8) == 4:6
    @test selected_channels(cfg([1.0, 3.0]), h, 8) == 4:6           # edges inclusive
    @test selected_channels(cfg([-9.0, 9.0]), h, 8) == 1:8
    @test selected_channels(cfg([0.5, 3.0]), hdr(5000.0, 1.0, -1000.0), 8) == 3:5   # descending axis
    @test_throws ArgumentError selected_channels(cfg([10.0, 12.0]), h, 8)
    @test_throws ArgumentError cfg([3.0, 1.0])
    @test_throws ArgumentError cfg([1.0, 2.0, 3.0])
    @test_throws ArgumentError cfg([NaN, 2.0])

    mktempdir() do dir
        rng = MersenneTwister(7)
        cube = zeros(Float32, 24, 20, 8)
        @inbounds for k in 1:8, j in 1:20, i in 1:24
            cube[i, j, k] = 2f0 * sin(0.3f0 * i * (1 + 0.1f0 * k) + 0.2f0 * j) + 0.5f0 * randn(rng, Float32)
        end
        full = joinpath(dir, "full.fits")
        cut = joinpath(dir, "cut.fits")
        FITS(f -> write(f, cube; header=h), full, "w")
        FITS(f -> write(f, cube[:, :, 4:6]; header=hdr(-2000.0, -2.0, 1000.0)), cut, "w")

        for noise_mode in (:scalar, :auto, :separable)
            base = (; keep_intermediate_chunks=false, noise=0.3, dn=4, ker=1, ny_chunks=2, min_valid=4,
                    noise_mode=noise_mode, noise_tile_size=8)
            a = run_pipeline(VGTConfig(; base..., input_file=full, velocity_range_kms=[0.5, 3.0],
                                       output_dir=joinpath(dir, "win_$(noise_mode)")))
            b = run_pipeline(VGTConfig(; base..., input_file=cut, output_dir=joinpath(dir, "cut_$(noise_mode)")))
            @test a.channels == 4:6 && b.channels == 1:3
            @test count(isfinite, a.psi) > 0
            for f in (:psi, :reliability, :fit_rbar2, :fit_support, :Q, :U, :I)
                @test isequal(getfield(a, f), getfield(b, f))
            end
            @test a.diagnostics == b.diagnostics
            @test a.header["CRPIX3"] == b.header["CRPIX3"] == -2.0
        end

        win = run_pipeline(VGTConfig(; keep_intermediate_chunks=false, noise=0.3, dn=4, ker=1, ny_chunks=1,
                                     min_valid=4, input_file=full, velocity_range_kms=[0.5, 3.0],
                                     output_dir=joinpath(dir, "out")))
        outputs = write_pipeline_outputs(win)
        ph = read_fits_header(outputs.psi_path)
        @test ph["VGTVLO"] == 0.5 && ph["VGTVHI"] == 3.0
        meta = read(outputs.metadata_path, String)
        @test occursin("channels_used = [4, 6]", meta)
        @test run_parameters(win.config)["vrange"] == "0.5_3.0"

        toml = joinpath(dir, "cfg.toml")
        write(toml, """
            [input]
            input_file = "$(full)"
            output_dir = "$(joinpath(dir, "toml"))"
            velocity_range_kms = [3, 7]
            [noise]
            noise = 0.3
            [vgt]
            dn = 4
            ker = 1
            min_valid = 4
            [execution]
            ny_chunks = 1
            """)
        @test load_config(toml).velocity_range_kms == [3.0, 7.0]
        # TOML parses `[]` as Vector{Any}; this once recursed in the positional constructor until the stack overflowed.
        write(toml, replace(read(toml, String), "velocity_range_kms = [3, 7]" => "velocity_range_kms = []"))
        @test load_config(toml).velocity_range_kms == Float64[]
        @test VGTConfig(input_file="c.fits", noise=1.0, dn=4, ker=1, min_valid=4, velocity_range_kms=Any[]).velocity_range_kms ==
              Float64[]
    end
end
