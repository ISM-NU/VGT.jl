using VGT
using DrWatson
using FITSIO
using JLD2
using Test
using TOML

@testset "next available output path" begin
    mktempdir() do dir
        @test VGT.next_available_output_path(dir, "psi_vgt_results.fits") ==
            joinpath(dir, "psi_vgt_results.fits")

        write(joinpath(dir, "psi_vgt_results.fits"), "first")
        @test VGT.next_available_output_path(dir, "psi_vgt_results.fits") ==
            joinpath(dir, "psi_vgt_results_1.fits")

        write(joinpath(dir, "psi_vgt_results_1.fits"), "second")
        @test VGT.next_available_output_path(dir, "psi_vgt_results.fits") ==
            joinpath(dir, "psi_vgt_results_2.fits")
    end
end

@testset "process_chunk uses exact path and records diagnostics" begin
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
            noise_mode=:auto,
            noise_tile_size=4,
            noise_sigma_multiplier=2.0,
        )
        mkpath(cfg.output_dir)

        expected_noise = VGT.build_noise_cube(
            cube,
            cfg.noise,
            cfg.noise_mode;
            tile_size=cfg.noise_tile_size,
            multiplier=cfg.noise_sigma_multiplier,
        )
        expected_q, expected_u, profile = VGT.VGT_QU_mv_2d_profiled(
            cube, cfg.dn, 1, size(cube, 3), expected_noise;
            min_valid=cfg.min_valid,
            fitter_mode=cfg.fitter_mode,
        )

        chunk_path = VGT.process_chunk(1, cfg, size(cube, 2), size(cube, 3))
        @test isfile(chunk_path)

        jldopen(chunk_path, "r") do file
            rbar_mask = isfinite.(profile[:fit_rbar2])
            @test maximum(abs.(file["q"] .- expected_q)) < 1e-10
            @test maximum(abs.(file["u"] .- expected_u)) < 1e-10
            @test file["fit_attempts"] == profile[:fit_attempts]
            @test file["fit_successes"] == profile[:fit_successes]
            @test file["fit_failures"] == profile[:fit_failures]
            @test file["valid_windows"] == profile[:valid_windows]
            @test file["fit_low_confidence_windows"] == profile[:fit_low_confidence_windows]
            @test maximum(abs.(file["fit_rbar2"][rbar_mask] .- profile[:fit_rbar2][rbar_mask])) < 1e-10
            @test maximum(abs.(file["fit_attempt_fraction"] .- profile[:fit_attempt_fraction])) < 1e-10
            @test maximum(abs.(file["fit_support"] .- profile[:fit_support])) < 1e-10
            @test !haskey(file, "fit_kappa")
            @test !haskey(file, "fit_ci95")
            @test !haskey(file, "fit_angle_entropy")
            @test !haskey(file, "fit_multimodal_score")
        end
    end
end

@testset "prefixed input filename" begin
    @test VGT.prefixed_input_filename("data/exp_raw/CO_data_reproject.fits", "psi_", ".fits") ==
        "psi_CO_data_reproject.fits"
    @test VGT.prefixed_input_filename("data/exp_raw/CO_data_reproject.fits", "psi_reliability_", ".fits") ==
        "psi_reliability_CO_data_reproject.fits"
    @test VGT.prefixed_input_filename("data/exp_raw/CO_data_reproject.fits", "psi_plot_", ".png") ==
        "psi_plot_CO_data_reproject.png"
    @test VGT.prefixed_input_filename("/tmp/galaxy.cube.fits", "psi_", ".fits") ==
        "psi_galaxy.cube.fits"
end

@testset "DrWatson run paths are parameter named" begin
    cfg = VGT.VGTConfig(
        input_file="data/exp_raw/example_cube.fits",
        output_dir="",
        keep_intermediate_chunks=false,
        noise=0.16,
        dn=20,
        ker=4,
        ny_chunks=1,
        min_valid=100,
        fitter_mode=:axial_vonmises,
    )

    name = VGT.run_savename(cfg)
    dir = VGT.default_output_dir(cfg)
    @test occursin("dn=20", name)
    @test occursin("fitter=axial_vonmises", name)
    @test occursin("noise_mode=scalar", name)
    @test occursin("cfg=$(VGT.run_parameter_hash(cfg))", name)
    @test length(name) <= 160
    @test occursin("input=example_cube", name)
    @test occursin(joinpath("data", "sims", "vgt"), dir)
    params = VGT.run_parameters(cfg)
    @test params["gradient"] == "central_checked"
    @test params["output_support"] == "smoothed_stokes_support"
    @test params["rbar2_weight"] == "false"
    for removed in ("finalize_mask", "histogram_edges", "gradient_sample",
                    "noise_spatial_exponent", "noise_channel_exponent")
        @test !haskey(params, removed)
    end
    @test VGT.resolved_output_dir(cfg) == dir
end

@testset "DrWatson project path helpers" begin
    @test VGT.project_config_path() == DrWatson.projectdir("config.toml")
    @test VGT.rawdatadir("cube.fits") == DrWatson.datadir("exp_raw", "cube.fits")
    @test VGT.processeddatadir("map.fits") == DrWatson.datadir("exp_pro", "map.fits")
    @test VGT.simulationdir("vgt", "run") == DrWatson.datadir("sims", "vgt", "run")
    @test VGT.experiment_default_dir("planck_report") == DrWatson.datadir("sims", "planck_report")
    params = Dict("dn" => 20, "fitter" => "axial_vonmises")
    @test VGT.parameterized_output_dir("root", params) == joinpath("root", DrWatson.savename(params))
    @test VGT.experiment_output_dir("scan", params) == DrWatson.datadir("sims", "scan", DrWatson.savename(params))
    @test VGT.experiment_summary_dir("out") == joinpath("out", "summaries")
    @test startswith(
        VGT.experiment_summary_path("out", params; prefix="summary"),
        joinpath("out", "summaries", "summary_"),
    )
end

@testset "DrWatson summary collector writes compact aggregate" begin
    mktempdir() do dir
        params = Dict("case" => "a")
        path = VGT.experiment_summary_path(dir, params; prefix="example")
        mkpath(dirname(path))
        jldsave(path; case="a", value=1.5)

        collected = VGT.collect_experiment_summaries(dir)
        @test collected == joinpath(dir, "collected_results.jld2")
        @test isfile(collected)
        loaded = JLD2.load(collected)
        @test haskey(loaded, "results")
        @test length(loaded["results"]) == 1
        @test loaded["results"][1]["case"] == "a"
        @test loaded["results"][1]["value"] == 1.5
    end
end

@testset "FITS wrapper writes and reads simple image products" begin
    mktempdir() do dir
        path = joinpath(dir, "wrapped.fits")
        data = Float32[1 2; 3 4]
        VGT.write_fits(path, data)
        loaded = VGT.read_fits_image(path)
        @test loaded == data
        header = VGT.read_fits_header(path)
        @test haskey(header, "NAXIS")
        @test header["NAXIS"] == 2
    end
end

@testset "FITS extension normalization matches write_fits compression trigger" begin
    @test VGT.normalize_fits_extension(".fits") == ".fits"
    @test VGT.normalize_fits_extension("fits") == ".fits"
    @test VGT.normalize_fits_extension(".fit") == ".fit"
    @test VGT.normalize_fits_extension(".fz") == ".fits.fz"
    @test VGT.normalize_fits_extension("fits.fz") == ".fits.fz"
    @test_throws ArgumentError VGT.normalize_fits_extension(".fits.gz")
end

@testset "FITS wrapper writes and reads compressed image products" begin
    if Sys.which("fpack") === nothing
        @test_skip "fpack not available"
    else
        mktempdir() do dir
            path = joinpath(dir, "wrapped.fits.fz")
            data = reshape(Float32.(1:100), 10, 10)
            VGT.write_fits(path, data)
            @test isfile(path)

            f = FITS(path)
            try
                @test length(f) == 2
                header = read_header(f[2])
                @test header["ZIMAGE"] == true
                @test header["ZCMPTYPE"] == "GZIP_1"
            finally
                close(f)
            end

            loaded = VGT.read_fits_image(path)
            @test loaded == data
            @test VGT.read_fits_image(path, 2:4, 3:5) == data[2:4, 3:5]

            header = VGT.read_fits_header(path)
            @test header["ZIMAGE"] == true
            @test header["ZNAXIS"] == 2
        end
    end
end

@testset "getcoor uses WCS-compatible FITS coordinates" begin
    header = FITSHeader(
        [
            "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
            "CTYPE1", "CTYPE2", "CTYPE3",
            "CRPIX1", "CRPIX2", "CRPIX3",
            "CRVAL1", "CRVAL2", "CRVAL3",
            "CDELT1", "CDELT2", "CDELT3",
        ],
        [
            true, -64, 3, 4, 3, 2,
            "LINEAR", "LINEAR", "LINEAR",
            1.0, 1.0, 1.0,
            10.0, 20.0, 30.0,
            0.5, 0.25, 2.0,
        ],
        fill("", 18),
    )

    ra, dec, vz = VGT.getcoor(header, 4, 3, 2)
    @test ra ≈ [10.0, 10.5, 11.0, 11.5]
    @test dec ≈ [20.0, 20.25, 20.5]
    @test vz ≈ [30.0, 32.0]
end

@testset "run_pipeline passes global Noise model through chunk jobs" begin
    mktempdir() do dir
        cube = zeros(Float32, 12, 12, 3)
        @inbounds for k in 1:3, j in 1:12, i in 1:12
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
            noise_mode=:auto,
        )

        run = VGT.run_pipeline(cfg)
        @test size(run.psi) == (12, 12)
        @test run.config.noise_mode == :auto
        @test run.diagnostics.fit_successes > 0
        @test run.performance.julia_threads >= 1
        @test run.performance.ny_chunks == 1
        @test run.performance.elapsed_seconds >= 0
    end
end

@testset "all-NaN angle map: outputs and metadata are written, plot skipped" begin
    mktempdir() do dir
        input = joinpath(dir, "empty.fits")
        cfg = VGT.VGTConfig(input_file=input, output_dir=joinpath(dir, "out"), noise=1.0, dn=4, ker=1, ny_chunks=1, min_valid=4)
        n = 6
        run = VGT.PipelineRun(
            config=cfg,
            psi=fill(NaN, n, n),
            reliability=fill(NaN32, n, n),
            Q=zeros(n, n),
            U=zeros(n, n),
            I=zeros(n, n),
            W=zeros(n, n),
            fit_rbar2=fill(NaN, n, n),
            fit_support=zeros(n, n),
            header=nothing,
            diagnostics=VGT.VGTDiagnostics(),
        )
        outputs = @test_logs (:warn, r"no finite angles; plot skipped") VGT.write_pipeline_outputs(run)
        @test isfile(outputs.psi_path)
        @test isfile(outputs.support_path)
        @test outputs.plot_path === nothing
        @test isfile(outputs.metadata_path)
        metadata = TOML.parsefile(outputs.metadata_path)
        @test !haskey(metadata["outputs"], "plot_path")
        @test metadata["outputs"]["psi_path"] == outputs.psi_path
    end
end
