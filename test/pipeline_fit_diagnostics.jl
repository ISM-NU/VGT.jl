using VGT
using FITSIO
using Test
using TOML

@testset "stitch and write fitter diagnostics outputs" begin
    mktempdir() do dir
        cube = zeros(Float32, 12, 10, 3)
        @inbounds for k in 1:3, j in 1:10, i in 1:12
            cube[i, j, k] = 3f0 + 0.2f0 * i + 0.1f0 * j + 0.05f0 * k
        end

        input = joinpath(dir, "cube.fits")
        FITS(input, "w") do f
            write(f, cube)
        end

        cfg = VGT.VGTConfig(
            input_file=input,
            output_dir=joinpath(dir, "out"),
            keep_intermediate_chunks=true,
            noise=1e-6,
            dn=4,
            ker=1,
            ny_chunks=1,
            min_valid=4,
            fitter_mode=:axial_vonmises,
        )
        mkpath(cfg.output_dir)

        chunk_path = VGT.process_chunk(1, cfg, size(cube, 2), size(cube, 3))
        @test isfile(chunk_path)

        q_stitched, u_stitched, intensity_stitched, diagnostics, fit_maps = VGT.stitch_chunk_outputs(cfg, size(cube, 1), size(cube, 2))
        @test size(fit_maps.fit_rbar2) == (size(cube, 1), size(cube, 2))
        @test size(fit_maps.fit_support) == (size(cube, 1), size(cube, 2))
        @test all(0 .<= fit_maps.fit_support .<= 1)
        @test any(isfinite, fit_maps.fit_rbar2)
        @test diagnostics.fit_attempts > 0
        @test diagnostics.fit_successes > 0
        @test diagnostics.fit_attempts == diagnostics.fit_successes + diagnostics.fit_failures

        finalized = VGT.finalize_vgt_products(q_stitched, u_stitched, intensity_stitched, cfg.ker)
        run = VGT.PipelineRun(
            config=cfg,
            psi=finalized.psi,
            reliability=Float32.(finalized.reliability),
            Q=q_stitched,
            U=u_stitched,
            I=intensity_stitched,
            fit_rbar2=fit_maps.fit_rbar2,
            fit_support=fit_maps.fit_support,
            header=nothing,
            diagnostics=diagnostics,
        )

        outputs = VGT.write_pipeline_outputs(run)
        @test isfile(outputs.psi_path)
        @test isfile(outputs.reliability_path)
        @test isfile(outputs.fit_rbar2_path)
        @test isfile(outputs.support_path)
        @test isfile(outputs.metadata_path)

        metadata = TOML.parsefile(outputs.metadata_path)
        @test metadata["config"]["input_file"] == input
        @test metadata["config"]["fitter_mode"] == "axial_vonmises"
        @test metadata["config"]["gradient_mode"] == "central_checked"
        @test !haskey(metadata["config"], "n_threads")
        @test metadata["performance"]["julia_threads"] >= 1
        @test metadata["diagnostics"]["fit_attempts"] == diagnostics.fit_attempts
        @test metadata["outputs"]["psi_path"] == outputs.psi_path
        @test metadata["outputs"]["support_path"] == outputs.support_path
        @test haskey(metadata, "gitcommit")
    end
end
