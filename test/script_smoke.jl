using FITSIO
using Test
using Random

@testset "analyze_inputs noise estimator is one-sigma for negative samples" begin
    analyzer = Module(:AnalyzeInputsSmoke)
    Core.eval(analyzer, :(using Test))
    Base.include(analyzer, joinpath(@__DIR__, "..", "scripts", "analyze_inputs.jl"))

    sigma = 0.42
    normal_mad_constant = 0.6744897501960817
    negative_half_sample = fill(-sigma * normal_mad_constant, 10, 10)

    @test analyzer.estimate_noise(negative_half_sample) ≈ sigma
    @test isnan(analyzer.estimate_noise(fill(-1.0, 9, 10)))
end

@testset "channel_quality_report distinguishes empty and noise-only channels" begin
    reporter = Module(:ChannelQualityReportSmoke)
    Core.eval(reporter, :(using Test))
    Base.include(reporter, joinpath(@__DIR__, "..", "scripts", "channel_quality_report.jl"))

    mktempdir() do dir
        cube = zeros(Float32, 4, 3, 3)
        cube[:, :, 2] .= reshape(Float32[0.05, -0.05, 0.04, -0.04, 0.03, -0.03,
                                         0.02, -0.02, 0.01, -0.01, 0.02, -0.02], 4, 3)
        cube[:, :, 3] .= 1.0f0
        header = FITSHeader(
            [
                "SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "NAXIS3",
                "CRPIX1", "CRPIX2", "CRPIX3",
                "CRVAL1", "CRVAL2", "CRVAL3",
                "CDELT1", "CDELT2", "CDELT3",
            ],
            [true, -32, 3, 4, 3, 3, 1.0, 1.0, 1.0, 0.0, 0.0, -1000.0, 1.0, 1.0, 1000.0],
            fill("", 15),
        )
        input = joinpath(dir, "channels.fits")
        FITS(input, "w") do f
            write(f, cube; header=header)
        end
        outdir = joinpath(dir, "out")
        csv_path = reporter.channel_quality_report(Dict(
            "input" => input,
            "noise" => "0.1",
            "output-dir" => outdir,
            "min-fill-3sigma" => "0.001",
            "min-positive-negative-ratio" => "1.25",
        ))
        text = read(csv_path, String)
        @test occursin("no_signal", text)
        @test occursin("true,false", split(text, '\n')[2])
        @test occursin("false,true", split(text, '\n')[3])
    end
end

@testset "analyze_inputs fallback handles no-negative samples" begin
    analyzer = Module(:AnalyzeInputsFallbackSmoke)
    Core.eval(analyzer, :(using Test))
    Base.include(analyzer, joinpath(@__DIR__, "..", "scripts", "analyze_inputs.jl"))

    sigma = 0.27
    normal_mad_constant = analyzer.NORMAL_MAD
    positive_half_sample = fill(sigma * normal_mad_constant, 10, 10)
    offset_symmetric_sample = vcat(
        fill(10.0 - sigma * normal_mad_constant, 5, 10),
        fill(10.0 + sigma * normal_mad_constant, 5, 10),
    )

    @test isnan(analyzer.estimate_noise(positive_half_sample))
    @test analyzer.estimate_fallback_noise(positive_half_sample) ≈ sigma
    @test analyzer.estimate_fallback_noise(offset_symmetric_sample) ≈ sigma
    @test isnan(analyzer.estimate_fallback_noise(zeros(10, 10)))
end

@testset "noise threshold evaluator helper functions are importable" begin
    evaluator = Module(:EvaluateNoiseThresholdsSmoke)
    Core.eval(evaluator, :(using Test))
    Base.include(evaluator, joinpath(@__DIR__, "..", "scripts", "evaluate_noise_thresholds.jl"))

    sigma = 0.31
    negative_half_sample = fill(-sigma * evaluator.MAD_NORMAL, 10, 10)

    @test evaluator.robust_negative_noise(negative_half_sample) ≈ sigma
    @test evaluator.sigma_value(:stationary, 1, 1, 1, 4, 4, 4) ≈ 0.17
    @test evaluator.sigma_value(:channel, 1, 1, 4, 4, 4, 4) >
          evaluator.sigma_value(:channel, 1, 1, 1, 4, 4, 4)
    @test evaluator.sigma_value(:spatial, 4, 4, 1, 4, 4, 4) >
          evaluator.sigma_value(:spatial, 1, 1, 1, 4, 4, 4)
end

@testset "validate_vgt_parameters parser and candidate expansion" begin
    validator = Module(:ValidateVGTParametersSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    opts = validator.parse_args([
        "--input", "data/exp_raw/example.fits",
        "--noise", "0.16",
        "--dn", "10,14,20",
        "--min-valid", "fixed100,frac0.50,80",
        "--fitter", "axial_vonmises",
        "--gradient-mode", "nan_edge_backward",
        "--output-dir", "validation_output/example",
    ])

    @test opts["input"] == "data/exp_raw/example.fits"
    @test opts["noise"] == "0.16"
    @test opts["dn"] == "10,14,20"
    @test opts["gradient-mode"] == "nan_edge_backward"
    @test opts["noise-sigma-multiplier"] == "1.0"
    @test validator.parse_dn_list("10,14,20") == [10, 14, 20]
    @test_throws ArgumentError validator.parse_dn_list("10,15")

    specs = validator.parse_min_valid_specs("fixed100,frac0.50,80")
    @test specs == [
        validator.MinValidSpec(:fixed, 100, 0.0, "fixed100"),
        validator.MinValidSpec(:fraction, 0, 0.50, "frac0.50"),
        validator.MinValidSpec(:fixed, 80, 0.0, "80"),
    ]

    candidates, skipped = validator.expand_candidates([10, 14], specs)
    @test validator.CandidateSpec(10, 100) in candidates
    @test validator.CandidateSpec(10, 50) in candidates
    @test validator.CandidateSpec(10, 80) in candidates
    @test validator.CandidateSpec(14, 100) in candidates
    @test validator.CandidateSpec(14, 98) in candidates
    @test validator.CandidateSpec(14, 80) in candidates
    @test any(c -> c.dn == 14 && c.min_valid == 98 && c.min_valid_label == "frac0.50", candidates)
    @test isempty(skipped)

    candidates2, skipped2 = validator.expand_candidates(
        [10],
        [validator.MinValidSpec(:fixed, 150, 0.0, "150")],
    )
    @test isempty(candidates2)
    @test skipped2 == ["dn=10 min_valid=150 skipped because min_valid > dn^2=100"]
end

@testset "validate_vgt_parameters axial summaries" begin
    validator = Module(:ValidateVGTParametersMetricsSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    a = [0.0 pi/2 NaN; pi/4 pi/2 0.0]
    b = [0.0 0.0 0.0; -pi/4 pi/2 NaN]
    delta = validator.axial_residual_map(a, b)

    @test delta[1, 1] ≈ 0.0
    @test abs(delta[1, 2]) ≈ pi / 2
    @test isnan(delta[1, 3])
    @test abs(delta[2, 1]) ≈ pi / 2
    @test delta[2, 2] ≈ 0.0
    @test isnan(delta[2, 3])

    summary = validator.summarize_abs_deg(delta)
    @test summary.n == 4
    @test summary.median ≈ 45.0
    @test summary.max ≈ 90.0

    layer = [0.1 0.2 0.3 0.4]
    values = [10.0 20.0 30.0 40.0]
    bins = validator.binned_abs_summaries(layer, values; nbins=2)
    @test length(bins) == 2
    @test bins[1].n == 2
    @test bins[1].median ≈ 15.0
    @test bins[2].n == 2
    @test bins[2].median ≈ 35.0
end

@testset "validate_vgt_parameters recommendation logic" begin
    validator = Module(:ValidateVGTParametersRecommendSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    metrics = [
        validator.CandidateMetrics(
            candidate=validator.CandidateSpec(10, 100),
            finite_psi_fraction=0.80,
            finite_reliability_fraction=0.80,
            fit_attempts=1000,
            fit_successes=950,
            fit_failures=50,
            valid_windows=1000,
            fit_success_rate=0.95,
            median_rbar2=0.20,
            median_reliability=0.25,
            low_confidence_fraction=0.65,
            adjacent_median_deg=12.0,
            adjacent_p95_deg=65.0,
            high_reliability_median_deg=9.0,
            high_reliability_p95_deg=40.0,
        ),
        validator.CandidateMetrics(
            candidate=validator.CandidateSpec(20, 100),
            finite_psi_fraction=0.72,
            finite_reliability_fraction=0.72,
            fit_attempts=900,
            fit_successes=870,
            fit_failures=30,
            valid_windows=900,
            fit_success_rate=0.9666666667,
            median_rbar2=0.38,
            median_reliability=0.44,
            low_confidence_fraction=0.30,
            adjacent_median_deg=4.0,
            adjacent_p95_deg=18.0,
            high_reliability_median_deg=2.5,
            high_reliability_p95_deg=9.0,
        ),
    ]

    rec = validator.recommend_candidate(metrics)
    @test rec.candidate == validator.CandidateSpec(20, 100)
    @test occursin("Recommended for review: dn=20, min_valid=100, rbar2_mode=baseline", rec.text)
    @test occursin("production config is not changed automatically", rec.text)
    @test occursin("Rejected dn=10, min_valid=100, rbar2_mode=baseline", rec.text)

    no_rec = validator.recommend_candidate([
        validator.CandidateMetrics(
            candidate=validator.CandidateSpec(10, 100),
            finite_psi_fraction=0.05,
            finite_reliability_fraction=0.05,
            fit_attempts=10,
            fit_successes=5,
            fit_failures=5,
            valid_windows=10,
            fit_success_rate=0.50,
            median_rbar2=0.05,
            median_reliability=0.05,
            low_confidence_fraction=0.95,
            adjacent_median_deg=30.0,
            adjacent_p95_deg=88.0,
            high_reliability_median_deg=NaN,
            high_reliability_p95_deg=NaN,
        ),
    ])
    @test no_rec.candidate === nothing
    @test occursin("No recommendation", no_rec.text)
end

@testset "compare_vgt_runs reports baseline support partitions" begin
    comparer = Module(:CompareVGTRunsSmoke)
    Core.eval(comparer, :(using Test))
    Base.include(comparer, joinpath(@__DIR__, "..", "scripts", "compare_vgt_runs.jl"))

    baseline = (;
        label="baseline",
        dir="baseline_dir",
        psi=[0.0 NaN; pi / 4 pi / 2],
        support=[0.8 NaN; 0.6 0.7],
        reliability=[0.9 NaN; 0.5 0.6],
        rbar2=[0.4 NaN; 0.2 0.3],
    )
    candidate = (;
        label="candidate",
        dir="candidate_dir",
        psi=[0.0 pi / 3; NaN pi / 2],
        support=[0.85 0.35; NaN 0.75],
        reliability=[0.95 0.45; NaN 0.65],
        rbar2=[0.45 0.25; NaN 0.35],
    )

    masks = comparer.finite_support_masks(baseline.psi, candidate.psi)
    @test masks.common == [true false; false true]
    @test masks.gained == [false true; false false]
    @test masks.lost == [false false; true false]

    row = comparer.support_partition_row(baseline, candidate)
    @test row["common_fraction"] ≈ 0.5
    @test row["gained_fraction"] ≈ 0.25
    @test row["lost_fraction"] ≈ 0.25
    @test row["common_n"] == 2
    @test row["gained_n"] == 1
    @test row["lost_n"] == 1
    @test row["common_delta_median_deg"] ≈ 0.0
    @test row["candidate_gained_support_median"] ≈ 0.35
    @test row["baseline_lost_support_median"] ≈ 0.6
    @test row["candidate_common_reliability_median"] ≈ 0.8
    @test row["baseline_lost_rbar2_median"] ≈ 0.2
end

@testset "validate_vgt_parameters metadata matching" begin
    validator = Module(:ValidateVGTParametersMetadataSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    cfg = validator.RunConfig(
        input="data/exp_raw/example.fits",
        noise=0.16,
        fitter=:axial_vonmises,
        gradient_mode=:central_checked,
        output_dir="validation_output/example",
        keep_intermediate_chunks=false,
        noise_mode=:scalar,
        noise_tile_size=128,
        noise_sigma_multiplier=1.0,
        ny_chunks=1,
    )
    candidate = validator.CandidateSpec(20, 100)
    metadata = validator.candidate_metadata(cfg, candidate)
    cdir = validator.candidate_dir(cfg, candidate)

    @test metadata["input"] == "data/exp_raw/example.fits"
    @test metadata["noise"] == "0.16"
    @test metadata["noise_sigma_multiplier"] == "1.0"
    @test metadata["fitter"] == "axial_vonmises"
    @test metadata["gradient_mode"] == "central_checked"
    @test metadata["dn"] == "20"
    @test metadata["min_valid"] == "100"
    @test startswith(cdir, "validation_output/example")
    @test occursin("dn=20", cdir)
    @test occursin("min_valid=100", cdir)
    @test occursin("noise=0.16", cdir)

    mktempdir() do dir
        path = joinpath(dir, "run_metadata.txt")
        validator.write_metadata(path, metadata)
        @test validator.metadata_matches(path, metadata)
        changed = copy(metadata)
        changed["dn"] = "28"
        @test !validator.metadata_matches(path, changed)
    end
end

@testset "validate_vgt_parameters candidate metric extraction" begin
    validator = Module(:ValidateVGTParametersExtractionSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    diagnostics = (
        fit_attempts=10,
        fit_successes=8,
        fit_failures=2,
        valid_windows=9,
        fit_low_confidence_windows=3,
    )
    run = validator.CandidateRun(
        candidate=validator.CandidateSpec(20, 100),
        dir="unused",
        psi=[0.0 0.1; NaN 0.2],
        reliability=[0.5 0.6; NaN 0.7],
        rbar2=[0.2 0.4; NaN 0.6],
        diagnostics=diagnostics,
    )
    metrics = validator.extract_candidate_metrics(run)
    @test metrics.finite_psi_fraction ≈ 0.75
    @test metrics.finite_reliability_fraction ≈ 0.75
    @test metrics.fit_success_rate ≈ 0.8
    @test metrics.median_rbar2 ≈ 0.4
    @test metrics.low_confidence_fraction ≈ 3 / 9

    baseline = validator.CandidateRun(
        candidate=validator.CandidateSpec(14, 100),
        dir="unused",
        psi=[0.0 0.0; 0.0 0.0],
        reliability=[0.9 0.9; 0.9 0.9],
        rbar2=[0.9 0.9; 0.9 0.9],
        diagnostics=diagnostics,
    )
    pair = validator.compare_candidate_pair(baseline, run)
    @test pair.a == validator.CandidateSpec(14, 100)
    @test pair.b == validator.CandidateSpec(20, 100)
    @test pair.global.n == 3
end

@testset "validate_vgt_parameters report writing" begin
    validator = Module(:ValidateVGTParametersReportSmoke)
    Core.eval(validator, :(using Test))
    Base.include(validator, joinpath(@__DIR__, "..", "scripts", "validate_vgt_parameters.jl"))

    metrics = [validator.CandidateMetrics(
        candidate=validator.CandidateSpec(20, 100),
        finite_psi_fraction=0.72,
        finite_reliability_fraction=0.70,
        fit_attempts=100,
        fit_successes=95,
        fit_failures=5,
        valid_windows=100,
        fit_success_rate=0.95,
        median_rbar2=0.40,
        median_reliability=0.50,
        low_confidence_fraction=0.20,
        adjacent_median_deg=4.0,
        adjacent_p95_deg=18.0,
        high_reliability_median_deg=2.0,
        high_reliability_p95_deg=8.0,
    )]
    pairs = [validator.PairwiseComparison(
        a=validator.CandidateSpec(14, 100),
        b=validator.CandidateSpec(20, 100),
        var"global"=validator.AbsSummary(n=10, median=4.0, p75=6.0, p90=8.0, p95=10.0, max=12.0, mean=5.0),
        high_reliability=validator.AbsSummary(n=5, median=2.0, p75=3.0, p90=4.0, p95=5.0, max=6.0, mean=2.5),
        low_reliability=validator.AbsSummary(n=5, median=8.0, p75=9.0, p90=10.0, p95=11.0, max=12.0, mean=8.5),
    )]
    rec = validator.Recommendation(candidate=validator.CandidateSpec(20, 100), text="Recommended: dn=20, min_valid=100")

    mktempdir() do dir
        validator.write_reports(dir, metrics, pairs, rec)
        @test isfile(joinpath(dir, "summary.txt"))
        @test isfile(joinpath(dir, "summary.csv"))
        @test isfile(joinpath(dir, "pairwise_residuals.csv"))
        @test occursin("This scan selects an internally stable VGT parameter setting",
                       read(joinpath(dir, "summary.txt"), String))
        @test occursin("dn,min_valid", read(joinpath(dir, "summary.csv"), String))
        @test occursin("a_dn,a_min_valid,a_min_valid_spec,a_rbar2_mode,b_dn,b_min_valid,b_min_valid_spec,b_rbar2_mode",
                       read(joinpath(dir, "pairwise_residuals.csv"), String))
    end
end
