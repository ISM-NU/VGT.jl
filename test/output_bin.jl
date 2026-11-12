using VGT
using FITSIO
using Test
using WCS: pix_to_world

@testset "output_bin: Q/U averaged onto a coarser grid with a matching WCS" begin
    a = [1.0 2.0 NaN 7.0 9.0;
         3.0 4.0 NaN 8.0 9.0;
         5.0 6.0 NaN NaN 9.0]
    m = block_mean(a, 2)
    @test size(m) == (1, 2)                       # partial edge blocks dropped
    @test m[1, 1] == 2.5
    @test m[1, 2] == 7.5                          # NaN entries skipped
    @test isnan(block_mean([NaN NaN; NaN NaN], 2)[1, 1])
    @test block_mean(a, 1; valid=isfinite.(a)) isa Matrix{Float64}

    # Each binned pixel centre must land on the centre of its input block.
    for (ctype, crval, cdelt) in ((("RA---TAN", "DEC--TAN"), (80.0, -69.0), (-0.1, 0.1)),
                                  (("GLON-CAR", "GLAT-CAR"), (279.0, -33.0), (-0.0953889974, 0.08)))
        keys = ["NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2"]
        h = FITSHeader(keys, Any[101, 100, ctype..., crval..., 51.0, 50.0, cdelt...], fill("", length(keys)))
        for b in (2, 3)
            hb = bin_header(h, b)
            @test hb["NAXIS1"] == 101 ÷ b && hb["NAXIS2"] == 100 ÷ b
            @test hb["CDELT1"] == b * cdelt[1]
            pix_b = Float64[1 7 hb["NAXIS1"]; 1 11 hb["NAXIS2"]]
            pix_in = (pix_b .- 0.5) .* b .+ 0.5
            wb = pix_to_world(VGT._celestial_wcs(hb), pix_b)
            win = pix_to_world(VGT._celestial_wcs(h), pix_in)
            @test maximum(abs.(wb .- win)) < 1e-9
        end
        @test h["CRPIX1"] == 51.0                 # input header untouched
    end
    cd_keys = ["CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"]
    hcd = bin_header(FITSHeader(cd_keys, Any["RA---TAN", "DEC--TAN", 10.0, 10.0, -0.1, 0.01, 0.02, 0.1],
                                fill("", 8)), 2)
    @test (hcd["CD1_1"], hcd["CD1_2"], hcd["CD2_1"], hcd["CD2_2"]) == (-0.2, 0.02, 0.04, 0.2)
    @test bin_header(nothing, 2) === nothing

    @test_throws ArgumentError VGTConfig(input_file="c.fits", noise=1.0, dn=4, ker=1, min_valid=4, output_bin=0)

    mktempdir() do dir
        cube = zeros(Float32, 16, 14, 6)
        @inbounds for k in 1:6, j in 1:14, i in 1:16
            cube[i, j, k] = 3f0 + 0.2f0 * i + 0.1f0 * j + 0.05f0 * k + 0.01f0 * i * j
        end
        keys = ["CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3", "CRPIX1", "CRPIX2", "CRPIX3",
                "CDELT1", "CDELT2", "CDELT3"]
        values = Any["RA---TAN", "DEC--TAN", "VRAD", 80.0, -69.0, 0.0, 8.0, 7.0, 1.0, -0.5, 0.5, 1000.0]
        sky = joinpath(dir, "sky.fits")
        FITS(sky, "w") do f
            write(f, cube; header=FITSHeader(keys, values, fill("", length(keys))))
        end
        base = (; input_file=sky, keep_intermediate_chunks=false, noise=1e-6, dn=4, ker=1,
                ny_chunks=1, min_valid=4, noise_mode=:auto)
        full = run_pipeline(VGTConfig(; base..., output_dir=joinpath(dir, "full")))
        binned = run_pipeline(VGTConfig(; base..., output_dir=joinpath(dir, "bin"), output_bin=3))
        @test size(binned.psi) == size(binned.reliability) == size(binned.fit_rbar2) == (5, 4)
        @test size(binned.Q) == size(full.Q)      # raw accumulators stay on the input grid
        @test count(isfinite, binned.psi) > 0
        @test binned.header["CDELT2"] == 1.5 && binned.header["CRPIX1"] == (8.0 - 0.5) / 3 + 0.5

        # psi comes from the averaged Q/U, which here equals the axial mean of the full-grid psi.
        fin = finalize_vgt_products(copy(full.Q), copy(full.U), copy(full.I), 1)
        @test isequal(binned.psi, bin_products(fin, full.fit_rbar2, full.fit_support, full.header, 3).psi)

        gal = run_pipeline(VGTConfig(; base..., output_dir=joinpath(dir, "gal"), output_bin=3, angle_frame=:galactic))
        @test isequal(gal.psi, pixel_to_position_angle(binned.psi, binned.header; frame=:galactic))

        outputs = write_pipeline_outputs(binned)
        h = read_fits_header(outputs.psi_path)
        @test h["NAXIS1"] == 5 && h["NAXIS2"] == 4 && h["VGTBIN"] == 3 && h["CDELT1"] == -1.5
        @test run_parameters(binned.config)["output_bin"] == 3
        @test !haskey(run_parameters(full.config), "output_bin")
        @test_throws ArgumentError run_pipeline(VGTConfig(; base..., output_dir=joinpath(dir, "big"), output_bin=20))

        toml = joinpath(dir, "cfg.toml")
        write(toml, """
            [input]
            input_file = "$(sky)"
            output_dir = "$(joinpath(dir, "toml"))"
            [noise]
            noise = 1e-6
            [vgt]
            dn = 4
            ker = 1
            min_valid = 4
            output_bin = 4
            [execution]
            ny_chunks = 1
            """)
        @test load_config(toml).output_bin == 4
    end
end
