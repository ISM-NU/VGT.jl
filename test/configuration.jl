using VGT
using Test

@testset "auto execution config" begin
    small_auto = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=:auto,
        min_valid=100,
        fitter_mode=:axial_vonmises,
    )
    resolved_small = VGT.resolve_execution_config(small_auto, 300, 300, 50)
    @test resolved_small.ny_chunks >= 1
    @test resolved_small.fitter_mode == :axial_vonmises
    @test resolved_small.gradient_mode == :central_checked

    explicit_execution = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=8,
        min_valid=100,
    )
    resolved_explicit = VGT.resolve_execution_config(explicit_execution, 300, 300, 50)
    @test resolved_explicit.ny_chunks == 8

    circular_candidate_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=4,
        min_valid=100,
        fitter_mode=:axial_vonmises,
    )
    resolved_circular_candidate = VGT.resolve_execution_config(circular_candidate_cfg, 300, 300, 50)
    @test resolved_circular_candidate.fitter_mode == :axial_vonmises

    weighted_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=4,
        min_valid=100,
        rbar2_weight=true,
    )
    resolved_weighted = VGT.resolve_execution_config(weighted_cfg, 300, 300, 50)
    @test resolved_weighted.rbar2_weight

    dynamic_noise_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=4,
        min_valid=100,
        fitter_mode=:axial_vonmises,
        noise_mode=:auto,
        noise_tile_size=64,
        noise_sigma_multiplier=2.5,
    )
    resolved_dynamic_noise = VGT.resolve_execution_config(dynamic_noise_cfg, 300, 300, 50)
    @test resolved_dynamic_noise.noise_mode == :auto
    @test resolved_dynamic_noise.noise_tile_size == 64
    @test resolved_dynamic_noise.noise_sigma_multiplier == 2.5

    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=1,
        min_valid=100,
        fitter_mode=:bad_mode,
    )

    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        noise_mode=:bad_noise_mode,
    )
    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        noise_sigma_multiplier=0.0,
    )

    @test_throws ArgumentError VGT.VGTConfig(
        "dummy.fits",
        "out",
        false,
        1.0,
        20,
        4,
        1,
        100,
        :bad_mode,
        false,
        :bad_noise_mode,
        128,
        :nan_before_smoothing,
        :smoothed_stokes_support,
        :nan_edge_backward,
    )

    alternative_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        zero_stokes_policy=:keep_zero_before_smoothing,
        output_support_mode=:direct_fit_support,
        gradient_mode=:nan_edge_backward,
    )
    @test alternative_cfg.zero_stokes_policy == :keep_zero_before_smoothing
    @test alternative_cfg.output_support_mode == :direct_fit_support
    @test alternative_cfg.gradient_mode == :nan_edge_backward
    @test VGT.chunk_padding_pixels(alternative_cfg) == 10

    sobel_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        gradient_mode=:sobel_checked,
    )
    @test sobel_cfg.gradient_mode == :sobel_checked

    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        zero_stokes_policy=:bad_zero_policy,
    )
    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        output_support_mode=:bad_output_support,
    )
    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        gradient_mode=:bad_gradient_mode,
    )
end

@testset "TOML config loading" begin
    mktempdir() do dir
        path = joinpath(dir, "config.toml")
        write(path, """
        [input]
        input_file = "data/exp_raw/example.fits"
        output_dir = ""
        keep_intermediate_chunks = false

        [noise]
        noise = 0.17
        noise_mode = "scalar"
        noise_tile_size = 64
        noise_sigma_multiplier = 2.0

        [vgt]
        dn = 20
        ker = 4
        min_valid = 100
        fitter_mode = "axial_vonmises"

        [execution]
        ny_chunks = "auto"
        """)

        cfg = VGT.load_config(path)
        @test cfg.input_file == "data/exp_raw/example.fits"
        @test cfg.output_dir == ""
        @test cfg.keep_intermediate_chunks == false
        @test cfg.noise == 0.17
        @test cfg.noise_mode === :scalar
        @test cfg.noise_tile_size == 64
        @test cfg.noise_sigma_multiplier == 2.0
        @test cfg.dn == 20
        @test cfg.ker == 4
        @test cfg.ny_chunks === :auto
        @test cfg.min_valid == 100
        @test cfg.fitter_mode === :axial_vonmises
        @test cfg.weight_file == ""
        @test VGT.weight_source(cfg) == cfg.input_file
        @test cfg.gradient_mode === :central_checked

        write(path, replace(read(path, String), "keep_intermediate_chunks = false" =>
                            "keep_intermediate_chunks = false\nweight_file = \"w.fits\""))
        wcfg = VGT.load_config(path)
        @test VGT.weight_source(wcfg) == "w.fits"
        @test VGT.resolve_execution_config(wcfg, 64, 64, 4).weight_file == "w.fits"
    end

    mktempdir() do dir
        path = joinpath(dir, "bad.toml")
        write(path, """
        [input]
        input_file = "data/exp_raw/example.fits"
        output_dir = ""

        [noise]
        noise = 0.17
        nois_mode = "scalar"

        [vgt]
        dn = 20
        ker = 4
        min_valid = 100

        [execution]
        ny_chunks = 1
        """)

        @test_throws ArgumentError VGT.load_config(path)
    end
end

const REMOVED_BASE = (input_file="dummy.fits", output_dir="out", noise=1.0, dn=20, ker=4)

@testset "removed options keep a clear failure" begin
    # Modes removed in 0.1.0 are rejected with a pointer to the archived repository.
    for kw in (
        (fitter_mode=:lsqfit_reference,),
        (fitter_mode=:legacy_lsqfit_reference,),
        (fitter_mode="lsqfit_reference",),
        (gradient_mode=:legacy_circshift,),
        (noise_mode=:separable_shrinkage,),
    )
        err = try
            VGT.VGTConfig(; REMOVED_BASE..., kw...)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("removed in VGT 0.1.0", err.msg)
        @test occursin("archived private repository", err.msg)
    end
    @test_throws ArgumentError VGT.VGTConfig("dummy.fits", "out", false, 1.0, 20, 4,
                                             1, 100, :lsqfit_reference)
    # An unknown (never existing) value keeps the generic message.
    err = try VGT.VGTConfig(; REMOVED_BASE..., fitter_mode=:bad_mode); nothing catch e e end
    @test err isa ArgumentError && !occursin("removed in", err.msg)

    # Removed keys: the old default is tolerated with a warning ...
    for kw in (
        (finalize_mask_mode=:strict_nonpositive,),
        (finalize_mask_mode="strict_nonpositive",),
        (histogram_edge_mode=:modern_range,),
        (gradient_sample_mode=:paired_finite,),
        (noise_spatial_exponent=1.0,),
        (noise_channel_exponent=1,),
    )
        cfg = @test_logs (:warn, r"obsolete") VGT.VGTConfig(; REMOVED_BASE..., kw...)
        @test cfg.fitter_mode === :axial_vonmises
    end
    # ... anything else is rejected.
    for kw in (
        (finalize_mask_mode=:legacy_nan_only,),
        (histogram_edge_mode=:legacy_linspace,),
        (gradient_sample_mode=:legacy_independent_nan_filter,),
        (noise_spatial_exponent=0.5,),
        (noise_channel_exponent=0.25,),
    )
        err = try
            VGT.VGTConfig(; REMOVED_BASE..., kw...)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("removed in VGT 0.1.0", err.msg)
    end
end

@testset "removed options in TOML config" begin
    mktempdir() do dir
        path = joinpath(dir, "config.toml")
        body(extra_vgt, extra_noise="") = """
        [input]
        input_file = "x.fits"
        output_dir = ""
        [noise]
        noise = 0.17
        $extra_noise
        [vgt]
        dn = 20
        ker = 4
        min_valid = 100
        $extra_vgt
        [execution]
        ny_chunks = 1
        """
        write(path, body("finalize_mask_mode = \"strict_nonpositive\"\nhistogram_edge_mode = \"modern_range\"\ngradient_sample_mode = \"paired_finite\"",
                         "noise_spatial_exponent = 1.0\nnoise_channel_exponent = 1.0"))
        cfg = @test_logs (:warn, r"obsolete") (:warn, r"obsolete") (:warn, r"obsolete") (:warn, r"obsolete") (:warn, r"obsolete") VGT.load_config(path)
        @test cfg.fitter_mode === :axial_vonmises

        for line in ("finalize_mask_mode = \"legacy_nan_only\"",
                     "histogram_edge_mode = \"legacy_linspace\"",
                     "gradient_sample_mode = \"legacy_independent_nan_filter\"",
                     "fitter_mode = \"lsqfit_reference\"",
                     "fitter_mode = \"legacy_lsqfit_reference\"",
                     "gradient_mode = \"legacy_circshift\"")
            write(path, body(line))
            @test_throws ArgumentError VGT.load_config(path)
        end
        for line in ("noise_spatial_exponent = 0.5", "noise_mode = \"separable_shrinkage\"")
            write(path, body("", line))
            @test_throws ArgumentError VGT.load_config(path)
        end
    end
end

@testset "auto execution formulas" begin
    @test VGT.choose_auto_threads(1) == 1
    @test VGT.choose_auto_threads(6) == 6
    @test VGT.choose_auto_threads(64) == 64

    @test VGT.choose_auto_chunks(300, 20, 4) == 8
    @test VGT.choose_auto_chunks(300, 10, 4) == 8
    @test VGT.choose_auto_chunks(12, 20, 4) == 1

    @test VGT.estimate_peak_chunk_bytes(100, 10, 20) ==
          100 * 20 * 10 * 4 + 5 * 100 * 20 * 8

end

@testset "auto execution resolver choices" begin
    cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=:auto,
        min_valid=100,
        fitter_mode=:axial_vonmises,
    )

    generous_budget = typemax(Int) ÷ 4
    resolved = VGT.choose_auto_execution(cfg, 300, 300, 50, 8, generous_budget)
    @test resolved.julia_threads == 8
    @test resolved.ny_chunks == 10

    small = VGT.choose_auto_execution(cfg, 30, 30, 10, 8, generous_budget)
    @test small.ny_chunks == 1

    peak_one = VGT.estimate_peak_chunk_bytes(300, 50, cld(300, 10) + 2 * VGT.chunk_padding_pixels(cfg))
    memory_limited = VGT.choose_auto_execution(cfg, 300, 300, 50, 8, peak_one)
    @test memory_limited.ny_chunks == 10

    explicit_chunks_cfg = VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=3,
        min_valid=100,
        fitter_mode=:axial_vonmises,
    )
    explicit = VGT.choose_auto_execution(explicit_chunks_cfg, 300, 300, 50, 8, generous_budget)
    @test explicit.julia_threads == 8
    @test explicit.ny_chunks == 3
end

@testset "memory budget follows the Slurm allocation" begin
    node_half = fld(Int(Sys.total_memory()), 2)
    @test VGT.memory_budget_bytes(Dict{String, String}()) == node_half
    @test VGT.memory_budget_bytes(Dict("SLURM_MEM_PER_NODE" => "1024")) == fld(1024 * 1024^2, 2)
    @test VGT.memory_budget_bytes(Dict("SLURM_MEM_PER_NODE" => "not-a-number")) == node_half
end
