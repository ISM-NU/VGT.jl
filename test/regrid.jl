using FITSIO

@testset "regrid helpers" begin
    @test VGT.source_axis_position(10.0, 100.0, 3.0, 2.0) ≈ -42.0
    @test VGT.source_axis_position(96.0, 100.0, 3.0, -2.0) ≈ 5.0

    bracket = VGT.linear_bracket(2.25, 4)
    @test bracket.lo == 2
    @test bracket.hi == 3
    @test bracket.w ≈ 0.25
    @test VGT.linear_bracket(0.0, 4) === nothing
    @test VGT.linear_bracket(5.0, 4) === nothing
    @test VGT.linear_bracket(4.0, 4) == (lo=4, hi=4, w=0.0)

    image = Float32[1 3; 5 7]
    @test VGT.bilinear_sample(image, 1.5, 1.5) ≈ 4.0f0
    @test VGT.bilinear_sample(image, 1.0, 2.0) ≈ 3.0f0
    image[2, 2] = NaN32
    @test isnan(VGT.bilinear_sample(image, 1.5, 1.5))

    target_header = Dict("CRVAL3" => 0.0, "CRPIX3" => 1.0, "CDELT3" => 1.0)
    source_header = Dict("CRVAL3" => -0.5, "CRPIX3" => 1.0, "CDELT3" => 2.0)
    @test VGT.velocity_overlap_indices(target_header, 5, source_header, 3) == 1:4
end

@testset "regrid script smoke test" begin
    regrid = Module(:RegridCubeSmoke)
    Core.eval(regrid, :(using FITSIO, Mmap, Printf, WCS, VGT))
    Base.include(regrid, joinpath(@__DIR__, "..", "scripts", "regrid_cube_to_reference.jl"))

    mktempdir() do dir
        function fixture_header(nx, ny, nv, crval3, cdelt3)
            keys = [
                "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
                "CTYPE1", "CTYPE2", "CTYPE3", "CUNIT3",
                "CRPIX1", "CRPIX2", "CRPIX3", "CRVAL1", "CRVAL2", "CRVAL3",
                "CDELT1", "CDELT2", "CDELT3",
            ]
            values = [
                true, -32, 3, nx, ny, nv,
                "LINEAR", "LINEAR", "LINEAR", "m/s",
                1.0, 1.0, 1.0, 10.0, 20.0, crval3,
                1.0, 1.0, cdelt3,
            ]
            return FITSHeader(keys, values, fill("", length(keys)))
        end

        source = zeros(Float32, 5, 4, 3)
        @inbounds for k in axes(source, 3), j in axes(source, 2), i in axes(source, 1)
            source[i, j, k] = 10f0 * i + 100f0 * j + 1000f0 * k
        end
        target = zeros(Float32, 3, 2, 3)
        source_path = joinpath(dir, "source.fits")
        target_path = joinpath(dir, "target.fits")
        output_path = joinpath(dir, "output.fits")
        FITS(source_path, "w") do f
            write(f, source; header=fixture_header(5, 4, 3, 0.0, 2.0))
        end
        FITS(target_path, "w") do f
            write(f, target; header=fixture_header(3, 2, 3, 1.0, 1.0))
        end

        regrid.regrid_cube(regrid.parse_args([
            "--source", source_path,
            "--reference", target_path,
            "--output", output_path,
            "--row-block", "1",
        ]))
        FITS(output_path) do f
            result = read(f[1])
            @test size(result) == (3, 2, 3)
            @test result[1, 1, 1] ≈ 1610f0
            @test result[3, 2, 3] ≈ 2730f0
            header = read_header(f[1])
            @test header["CRVAL3"] ≈ 1.0
            @test header["CDELT3"] ≈ 1.0
        end
    end
end
