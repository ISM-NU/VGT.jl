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

    # Every key, in the table the example config documents, reaches its field.
    mktempdir() do dir
        path = joinpath(dir, "all.toml")
        write(path, """
        [input]
        input_file = "c.fits"
        output_dir = "o"
        keep_intermediate_chunks = false
        weight_file = "w.fits"
        velocity_range_kms = [1.0, 5.0]

        [noise]
        noise = 0.3
        noise_mode = "channel"
        noise_tile_size = 32
        noise_sigma_multiplier = 2.0
        gradient_whiten = true
        gradient_noise_cov = [0.3, 0.2, 0.05]

        [vgt]
        dn = 12
        ker = 2
        min_valid = 10
        fitter_mode = "axial_vonmises"
        zero_stokes_policy = "keep_zero_before_smoothing"
        output_support_mode = "direct_fit_support"
        gradient_mode = "sobel_checked"
        rbar2_weight = true
        rbar2_gate = 0.3
        angle_frame = "galactic"
        output_bin = 2
        edge_mask = true

        [execution]
        ny_chunks = 4
        """)
        expected = VGT.VGTConfig(
            input_file="c.fits", output_dir="o", keep_intermediate_chunks=false,
            weight_file="w.fits", velocity_range_kms=[1.0, 5.0], noise=0.3, noise_mode=:channel,
            noise_tile_size=32, noise_sigma_multiplier=2.0, gradient_whiten=true,
            gradient_noise_cov=[0.3, 0.2, 0.05], dn=12, ker=2, min_valid=10,
            fitter_mode=:axial_vonmises, zero_stokes_policy=:keep_zero_before_smoothing,
            output_support_mode=:direct_fit_support, gradient_mode=:sobel_checked,
            rbar2_weight=true, rbar2_gate=0.3, angle_frame=:galactic, output_bin=2,
            edge_mask=true, ny_chunks=4)
        loaded = VGT.load_config(path)
        for f in fieldnames(VGT.VGTConfig)
            @test getfield(loaded, f) == getfield(expected, f)
        end
        # No field is left at its default, so a key that is dropped shows up here.
        defaults = VGT.VGTConfig(input_file="c.fits", noise=0.3, dn=12, ker=2)
        # gradient_noise_cov_file is exclusive with gradient_noise_cov, so it gets its own load.
        @test all(f -> getfield(loaded, f) != getfield(defaults, f) ||
                       f in (:input_file, :noise, :dn, :ker, :fitter_mode, :gradient_noise_cov_file),
                  fieldnames(VGT.VGTConfig))
        write(path, replace(read(path, String), "gradient_noise_cov = [0.3, 0.2, 0.05]" =>
                            "gradient_noise_cov_file = \"cov.fits\""))
        @test VGT.load_config(path).gradient_noise_cov_file == "cov.fits"
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

@testset "tuple constructor validates like keywords" begin
    cfg = VGT.VGTConfig(input_file="dummy.fits", output_dir="out", noise=1.0, dn=20, ker=4)
    fields = collect(fieldnames(VGT.VGTConfig))
    values(; changes...) = Tuple(haskey(changes, f) ? changes[f] : getfield(cfg, f) for f in fields)
    @test VGT.VGTConfig(values()).dn == 20
    @test VGT.VGTConfig(values(fitter_mode="axial_vonmises")).fitter_mode === :axial_vonmises
    @test_throws ArgumentError VGT.VGTConfig(values(fitter_mode=:bogus))
    @test_throws ArgumentError VGT.VGTConfig(values(ny_chunks=0))
    @test_throws ArgumentError VGT.VGTConfig(values(ny_chunks=:foo))
end

@testset "config defaults and required TOML keys" begin
    cfg = VGT.VGTConfig(input_file="in.fits", noise=1.0, dn=20, ker=4)
    defaults = (
        input_file = "in.fits", output_dir = "subcube", keep_intermediate_chunks = true,
        noise = 1.0, dn = 20, ker = 4, ny_chunks = 1, min_valid = 100,
        fitter_mode = :axial_vonmises, rbar2_weight = false, noise_mode = :scalar,
        noise_tile_size = 128, zero_stokes_policy = :nan_before_smoothing,
        output_support_mode = :smoothed_stokes_support, gradient_mode = :central_checked,
        noise_sigma_multiplier = 1.0, rbar2_gate = 0.0, weight_file = "",
        angle_frame = :pixel, output_bin = 1, velocity_range_kms = Float64[],
        gradient_whiten = false, gradient_noise_cov = Float64[],
        gradient_noise_cov_file = "", edge_mask = false)
    @test NamedTuple{fieldnames(VGT.VGTConfig)}(Tuple(getfield(cfg, f) for f in fieldnames(VGT.VGTConfig))) == defaults

    # A TOML file must set these even though some have keyword defaults: a run
    # should not silently pick up a default for its output folder, chunking or window threshold.
    required = (input = ("input_file" => "\"in.fits\"", "output_dir" => "\"out\""),
                noise = ("noise" => "1.0",),
                vgt = ("dn" => "20", "ker" => "4", "min_valid" => "100"),
                execution = ("ny_chunks" => "1",))
    function toml(; skip="")
        io = IOBuffer()
        for (table, entries) in pairs(required)
            println(io, "[$table]")
            for (key, val) in entries
                key == skip || println(io, "$key = $val")
            end
        end
        return String(take!(io))
    end

    mktempdir() do dir
        path = joinpath(dir, "config.toml")
        write(path, toml())
        @test VGT.load_config(path).min_valid == 100
        for entries in required, (key, _) in entries
            write(path, toml(skip=key))
            @test_throws ArgumentError VGT.load_config(path)
        end
    end
end

@testset "auto execution trades chunks for memory" begin
    cfg = VGT.VGTConfig(input_file="dummy.fits", output_dir="out", noise=1.0, dn=20, ker=4,
                        ny_chunks=:auto)
    nx, ny, nv = 300, 300, 50
    pad = VGT.chunk_padding_pixels(cfg)
    chunks_for(budget) = VGT.choose_auto_execution(cfg, nx, ny, nv, 8, budget).ny_chunks
    peak(chunks) = VGT.estimate_peak_chunk_bytes(nx, nv, cld(ny, chunks) + 2 * pad)

    budgets = [typemax(Int) ÷ 4, 12_000_000, 6_000_000, 4_000_000, 3_000_000]
    chunks = chunks_for.(budgets)
    @test all(c -> 1 <= c <= ny, chunks)
    @test all(peak.(chunks) .<= budgets)
    @test issorted(chunks)   # a tighter budget never means fewer chunks
    @test chunks[end] > chunks[1]
    @test VGT.estimate_peak_chunk_bytes(100, 10, 20) == 160_000   # 100*20*10*4 + 5*100*20*8

    # Too little memory for even one-row chunks is an error, not a silent overrun.
    @test_throws ArgumentError chunks_for(1_000)

    # An explicit chunk count that fits the budget is kept.
    fixed = VGT.config_with(cfg; ny_chunks=3)
    @test VGT.choose_auto_execution(fixed, nx, ny, nv, 8, typemax(Int) ÷ 4).ny_chunks == 3
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
