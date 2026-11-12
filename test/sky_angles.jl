using VGT
using FITSIO
using Test

function _sky_test_header(ctype1, ctype2, crval, cdelt; naxis3=false)
    keys = ["CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "CUNIT1", "CUNIT2"]
    values = Any[ctype1, ctype2, crval[1], crval[2], 51.0, 51.0, cdelt[1], cdelt[2], "deg", "deg"]
    if ctype1 == "RA---TAN"
        push!(keys, "RADESYS"); push!(values, "ICRS")
    end
    if naxis3
        append!(keys, ["NAXIS3", "CTYPE3", "CRVAL3", "CRPIX3", "CDELT3", "CUNIT3", "PC3_3", "SPECSYS",
                       "BUNIT", "BMAJ", "BMIN", "BPA"])
        append!(values, Any[10, "VRAD", 0.0, 1.0, 1000.0, "m/s", 1.0, "LSRK", "K", 0.01, 0.01, 0.0])
    end
    return FITSHeader(keys, values, fill("", length(keys)))
end

const SKY_TEST_HEADERS = Dict(
    # LMC-like equatorial TAN grid, 10 deg across.
    "tan" => _sky_test_header("RA---TAN", "DEC--TAN", (80.0, -69.0), (-0.1, 0.1)),
    # Oblique Galactic CAR grid with non-square pixels, like the GASS LMC cube.
    "car" => _sky_test_header("GLON-CAR", "GLAT-CAR", (279.0, -33.0), (-0.0953889974, 0.08)),
)

# (grid, FITS x, FITS y, frame, IAU position angle chi, pixel-frame psi). psi was
# computed with astropy by stepping along chi on the sky and projecting through
# the WCS, the method of `galactic_angle_to_pixel` in
# scripts/compare_vgt_polarization.py in ISM-NU/LMC-foreground (1e-3 deg step).
const SKY_ANGLE_FIXTURES = [
    ("tan", 51, 51, :galactic, 0.7, -2.249655743906),
    ("tan", 51, 51, :icrs, 0.7, 2.270819539521),
    ("tan", 1, 51, :galactic, 0.7, -2.261210270065),
    ("tan", 1, 51, :icrs, 0.7, 2.495362530866),
    ("tan", 101, 51, :galactic, 0.7, -2.239141837802),
    ("tan", 101, 51, :icrs, 0.7, 2.049649351371),
    ("tan", 51, 1, :galactic, 0.7, -2.193071017771),
    ("tan", 51, 1, :icrs, 0.7, 2.268958912788),
    ("tan", 101, 101, :galactic, 0.7, -2.302333165146),
    ("tan", 101, 101, :icrs, 0.7, 2.089278648134),
    ("tan", 101, 101, :galactic, 2.0, -1.003516303873),
    ("tan", 101, 101, :icrs, 2.0, -2.888677396658),
    ("car", 51, 51, :galactic, 0.7, 2.185812202268),
    ("car", 51, 51, :icrs, 0.7, 0.613716050922),
    ("car", 1, 51, :galactic, 0.7, 2.237876556719),
    ("car", 1, 51, :icrs, 0.7, 0.564878485670),
    ("car", 101, 51, :galactic, 0.7, 2.134674948030),
    ("car", 101, 51, :icrs, 0.7, 0.644279883290),
    ("car", 51, 1, :galactic, 0.7, 2.186962302678),
    ("car", 51, 1, :icrs, 0.7, 0.427423227294),
    ("car", 101, 101, :galactic, 0.7, 2.137852346998),
    ("car", 101, 101, :icrs, 0.7, 0.791250784355),
    ("car", 101, 101, :galactic, 2.0, -2.700921654821),
    ("car", 101, 101, :icrs, 2.0, 1.942399566580),
]

@testset "pixel-frame psi converts to IAU position angles (astropy reference)" begin
    for (grid, x, y, frame, chi, psi_pix) in SKY_ANGLE_FIXTURES
        psi = fill(NaN, 101, 101)
        psi[x, y] = psi_pix
        out = pixel_to_position_angle(psi, SKY_TEST_HEADERS[grid]; frame=frame)
        @test abs(headless_angle_difference(out[x, y], chi)) < 1e-4
        @test count(isfinite, out) == 1
    end
end

@testset "celestial_header keeps the sky WCS and drops cube-only cards" begin
    cube = _sky_test_header("RA---TAN", "DEC--TAN", (80.0, -69.0), (-0.1, 0.1); naxis3=true)
    h = celestial_header(cube)
    @test h["WCSAXES"] == 2
    for key in ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CDELT1", "CDELT2", "RADESYS")
        @test h[key] == cube[key]
    end
    for key in ("NAXIS3", "CTYPE3", "CRVAL3", "CRPIX3", "CDELT3", "CUNIT3", "PC3_3", "SPECSYS",
                "BUNIT", "BMAJ", "BMIN", "BPA")
        @test !haskey(h, key)
    end
    @test h["INBMAJ"] == 0.01
    # The cube header gives the same position angles as its 2D product header.
    psi = fill(0.3, 101, 101)
    @test pixel_to_position_angle(psi, cube) == pixel_to_position_angle(psi, h)
end

@testset "VGT product headers record their own unit, convention and parameters" begin
    cfg = VGT.VGTConfig(
        input_file="cube.fits", output_dir="", keep_intermediate_chunks=false,
        noise=0.16, dn=20, ker=4, ny_chunks=1, min_valid=100, fitter_mode=:axial_vonmises,
    )
    cube = _sky_test_header("RA---TAN", "DEC--TAN", (80.0, -69.0), (-0.1, 0.1); naxis3=true)
    psi_header = vgt_product_header(cube, cfg, :psi)
    @test psi_header["BUNIT"] == "rad"
    @test psi_header["PSIFRAME"] == "PIXEL"
    @test psi_header["VGTDN"] == 20
    @test psi_header["VGTKER"] == 4.0
    @test psi_header["VGTKFWHM"] ≈ 4 * 2.3548200450309493
    @test !haskey(psi_header, "CTYPE3")
    support_header = vgt_product_header(cube, cfg, :support)
    @test support_header["BUNIT"] == ""
    @test !haskey(support_header, "PSIFRAME")

    mktempdir() do dir
        path = joinpath(dir, "psi.fits")
        write_fits(path, rand(101, 101); header=psi_header)
        written = read_fits_header(path)
        @test written["NAXIS"] == 2
        @test written["WCSAXES"] == 2
        @test written["CTYPE1"] == "RA---TAN"
        @test !haskey(written, "CTYPE3")
        @test get_comment(written, "PSIFRAME") != ""
    end
end

@testset "pipeline products: unchanged data, 2D sky headers, commentary cards kept" begin
    mktempdir() do dir
        cube = zeros(Float32, 12, 12, 3)
        @inbounds for k in 1:3, j in 1:12, i in 1:12
            cube[i, j, k] = 3f0 + 0.2f0 * i + 0.1f0 * j + 0.05f0 * k
        end
        keys = ["CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3", "CRPIX1", "CRPIX2", "CRPIX3",
                "CDELT1", "CDELT2", "CDELT3", "BUNIT", "BMAJ", "COMMENT", "COMMENT", "HISTORY"]
        values = Any["RA---TAN", "DEC--TAN", "VRAD", 80.0, -69.0, 0.0, 6.0, 6.0, 1.0,
                     -0.01, 0.01, 1000.0, "K", 0.02, nothing, nothing, nothing]
        comments = vcat(fill("", 14), ["first note", "second note", "made by a test"])
        input = joinpath(dir, "cube.fits")
        FITS(input, "w") do f
            write(f, cube; header=FITSHeader(keys, values, comments))
        end
        cfg = VGT.VGTConfig(
            input_file=input, output_dir=joinpath(dir, "out"), keep_intermediate_chunks=false,
            noise=1e-6, dn=4, ker=1, ny_chunks=1, min_valid=4, noise_mode=:auto,
        )
        run = VGT.run_pipeline(cfg)
        outputs = write_pipeline_outputs(run)
        for (path, data) in ((outputs.psi_path, run.psi), (outputs.reliability_path, run.reliability),
                             (outputs.fit_rbar2_path, run.fit_rbar2),
                             (outputs.support_path, Float32.(run.fit_support)))
            @test isequal(read_fits_image(path), data)
            h = read_fits_header(path)
            @test h["NAXIS"] == 2 && h["CTYPE1"] == "RA---TAN" && h["VGTDN"] == 4
            @test !haskey(h, "CTYPE3") && !haskey(h, "BMAJ")
            texts = [h.comments[i] for i in eachindex(h.keys) if h.keys[i] in ("COMMENT", "HISTORY")]
            @test all(t -> any(c -> occursin(t, c), texts), ("first note", "second note", "made by a test"))
        end
        @test read_fits_header(outputs.psi_path)["PSIFRAME"] == "PIXEL"
    end
end
