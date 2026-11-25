using VGT
using FITSIO
using Random
using Test

# A run split into y chunks must give exactly the one-piece result: each chunk reads enough rows around its own rows
# for every window and every gradient stencil in it (chunk_padding_pixels).
@testset "chunked runs equal the one-piece run" begin
    keys = ["CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3", "CRPIX1", "CRPIX2", "CRPIX3",
            "CDELT1", "CDELT2", "CDELT3", "CUNIT3"]
    h = FITSHeader(keys, Any["RA---TAN", "DEC--TAN", "VRAD", 80.0, -69.0, 0.0, 12.0, 20.0, 1.0,
                             -0.5, 0.5, 1000.0, "m/s"], fill("", length(keys)))
    mktempdir() do dir
        rng = MersenneTwister(11)
        cube = zeros(Float32, 24, 40, 5)
        @inbounds for k in 1:5, j in 1:40, i in 1:24
            cube[i, j, k] = 3f0 + 2f0 * sin(0.3f0 * i + 0.25f0 * j * (1 + 0.1f0 * k)) + 0.5f0 * randn(rng, Float32)
        end
        path = joinpath(dir, "c.fits")
        FITS(f -> write(f, cube; header=h), path, "w")
        for mode in (:central_checked, :sobel_checked, :isotropic_checked, :nan_edge_backward)
            base = (; input_file=path, keep_intermediate_chunks=false, noise=1.0, dn=6, ker=1, min_valid=4,
                    gradient_mode=mode)
            one = run_pipeline(VGTConfig(; base..., ny_chunks=1, output_dir=joinpath(dir, "one_$mode")))
            @test count(isfinite, one.psi) > 0
            for n in (3, 4)
                split = run_pipeline(VGTConfig(; base..., ny_chunks=n, output_dir=joinpath(dir, "split_$(mode)_$n")))
                for f in (:psi, :reliability, :Q, :U, :I, :W, :fit_rbar2, :fit_support)
                    @test isequal(getfield(split, f), getfield(one, f))
                end
            end
        end
    end
end
