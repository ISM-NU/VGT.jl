using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using JLD2
using Printf
using Statistics
using VGT

Base.@kwdef struct MinValidSpec
    kind::Symbol
    value::Int
    fraction::Float64
    label::String
end

Base.:(==)(a::MinValidSpec, b::MinValidSpec) =
    a.kind == b.kind && a.value == b.value && a.fraction == b.fraction && a.label == b.label

Base.@kwdef struct Rbar2Spec
    weight::Bool
    gate::Float64
    label::String
end

struct CandidateSpec
    dn::Int
    min_valid::Int
    min_valid_label::String
    rbar2_weight::Bool
    rbar2_gate::Float64
    rbar2_label::String
end

const CANDIDATE_KER = 4
const CANDIDATE_ZERO_STOKES_POLICY = :nan_before_smoothing
const CANDIDATE_OUTPUT_SUPPORT_MODE = :smoothed_stokes_support
const CANDIDATE_GRADIENT_MODE = :central_checked

CandidateSpec(dn::Int, min_valid::Int) =
    CandidateSpec(dn, min_valid, string(min_valid), false, 0.0, "baseline")
CandidateSpec(dn::Int, min_valid::Int, min_valid_label::String) =
    CandidateSpec(dn, min_valid, min_valid_label, false, 0.0, "baseline")

Base.:(==)(a::CandidateSpec, b::CandidateSpec) =
    a.dn == b.dn && a.min_valid == b.min_valid &&
    a.rbar2_weight == b.rbar2_weight && a.rbar2_gate == b.rbar2_gate
Base.hash(c::CandidateSpec, h::UInt) = hash((c.dn, c.min_valid, c.rbar2_weight, c.rbar2_gate), h)

function default_options()
    return Dict(
        "input" => nothing,
        "noise" => nothing,
        "dn" => "10,14,20,28",
        "min-valid" => "auto",
        "rbar2-mode" => "baseline,weight,gate0.50",
        "fitter" => "axial_vonmises",
        "gradient-mode" => String(CANDIDATE_GRADIENT_MODE),
        "output-dir" => experiment_default_dir("vgt_parameter_scan"),
        "reference-dn" => nothing,
        "force" => "false",
        "keep-intermediate-chunks" => "false",
        "noise-mode" => "scalar",
        "noise-tile-size" => "128",
        "noise-sigma-multiplier" => "1.0",
        "ny-chunks" => "1",
        "n-threads" => nothing,
    )
end

function parse_args(args)
    opts = default_options()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("Unexpected positional argument: $arg"))
        key = arg[3:end]
        haskey(opts, key) || throw(ArgumentError("Unknown option --$key"))
        if key == "force"
            opts[key] = "true"
            i += 1
        else
            i == length(args) && throw(ArgumentError("Missing value for --$key"))
            opts[key] = args[i + 1]
            i += 2
        end
    end
    return opts
end

function parse_dn_list(raw::AbstractString)
    values = Int[]
    for part in split(raw, ',')
        text = strip(part)
        isempty(text) && continue
        dn = parse(Int, text)
        dn >= 2 || throw(ArgumentError("dn must be >= 2; got $dn"))
        iseven(dn) || throw(ArgumentError("dn must be even; got $dn"))
        push!(values, dn)
    end
    isempty(values) && throw(ArgumentError("--dn must contain at least one value"))
    return unique(sort(values))
end

function parse_min_valid_specs(raw::AbstractString)
    specs = MinValidSpec[]
    for part in split(raw, ',')
        text = strip(part)
        isempty(text) && continue
        if text == "auto"
            for frac in (0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
                push!(specs, MinValidSpec(;
                    kind=:fraction,
                    value=0,
                    fraction=frac,
                    label=@sprintf("auto_frac%.2f", frac),
                ))
            end
        elseif text == "fixed100"
            push!(specs, MinValidSpec(; kind=:fixed, value=100, fraction=0.0, label=text))
        elseif startswith(text, "frac")
            frac = parse(Float64, text[5:end])
            0.0 < frac <= 1.0 || throw(ArgumentError("fraction min-valid must be in (0,1]; got $text"))
            push!(specs, MinValidSpec(; kind=:fraction, value=0, fraction=frac, label=text))
        else
            value = parse(Int, text)
            value >= 1 || throw(ArgumentError("min_valid must be >= 1; got $value"))
            push!(specs, MinValidSpec(; kind=:fixed, value=value, fraction=0.0, label=text))
        end
    end
    isempty(specs) && throw(ArgumentError("--min-valid must contain at least one value"))
    return specs
end

function parse_rbar2_specs(raw::AbstractString)
    specs = Rbar2Spec[]
    for part in split(raw, ',')
        text = strip(part)
        isempty(text) && continue
        if text == "baseline"
            push!(specs, Rbar2Spec(; weight=false, gate=0.0, label=text))
        elseif text == "weight"
            push!(specs, Rbar2Spec(; weight=true, gate=0.0, label=text))
        elseif startswith(text, "gate")
            gate = parse(Float64, text[5:end])
            0.0 < gate <= 1.0 || throw(ArgumentError("rbar2 gate must be in (0,1]; got $text"))
            push!(specs, Rbar2Spec(; weight=false, gate=gate, label=text))
        else
            throw(ArgumentError("Unsupported --rbar2-mode value $text; use baseline, weight, or gate0.50"))
        end
    end
    isempty(specs) && throw(ArgumentError("--rbar2-mode must contain at least one value"))
    return specs
end

function min_valid_value(spec::MinValidSpec, dn::Int)
    if spec.kind === :fixed
        return spec.value
    elseif spec.kind === :fraction
        return max(1, round(Int, spec.fraction * dn * dn))
    end
    throw(ArgumentError("Unsupported min_valid spec kind $(spec.kind)"))
end

function expand_candidates(dn_values::Vector{Int}, specs::Vector{MinValidSpec},
                           rbar2_specs::Vector{Rbar2Spec}=Rbar2Spec[Rbar2Spec(; weight=false, gate=0.0, label="baseline")])
    candidates = CandidateSpec[]
    skipped = String[]
    seen = Set{CandidateSpec}()
    for dn in dn_values
        for spec in specs
            min_valid = min_valid_value(spec, dn)
            if min_valid > dn * dn
                push!(skipped, "dn=$dn min_valid=$min_valid skipped because min_valid > dn^2=$(dn * dn)")
                continue
            end
            for rbar2_spec in rbar2_specs
                candidate = CandidateSpec(
                    dn,
                    min_valid,
                    spec.label,
                    rbar2_spec.weight,
                    rbar2_spec.gate,
                    rbar2_spec.label,
                )
                if !(candidate in seen)
                    push!(candidates, candidate)
                    push!(seen, candidate)
                end
            end
        end
    end
    return candidates, skipped
end

Base.@kwdef struct RunConfig
    input::String
    noise::Float64
    fitter::Symbol
    gradient_mode::Symbol
    output_dir::String
    keep_intermediate_chunks::Bool
    noise_mode::Symbol
    noise_tile_size::Int
    noise_sigma_multiplier::Float64
    ny_chunks::Union{Int, Symbol}
end

function parse_bool(raw::AbstractString)
    raw in ("true", "yes", "1") && return true
    raw in ("false", "no", "0") && return false
    throw(ArgumentError("Expected boolean string, got $raw"))
end

function parse_execution_value(raw::AbstractString)
    raw == "auto" && return :auto
    return parse(Int, raw)
end

function run_config_from_options(opts)
    opts["input"] === nothing && throw(ArgumentError("--input is required"))
    opts["noise"] === nothing && throw(ArgumentError("--noise is required"))
    opts["n-threads"] === nothing ||
        @warn "`--n-threads` is ignored; start Julia with `--threads=N` or set JULIA_NUM_THREADS" maxlog=1
    return RunConfig(
        input=opts["input"],
        noise=parse(Float64, opts["noise"]),
        fitter=Symbol(opts["fitter"]),
        gradient_mode=Symbol(opts["gradient-mode"]),
        output_dir=opts["output-dir"],
        keep_intermediate_chunks=parse_bool(opts["keep-intermediate-chunks"]),
        noise_mode=Symbol(opts["noise-mode"]),
        noise_tile_size=parse(Int, opts["noise-tile-size"]),
        noise_sigma_multiplier=parse(Float64, opts["noise-sigma-multiplier"]),
        ny_chunks=parse_execution_value(opts["ny-chunks"]),
    )
end

function candidate_dir(cfg::RunConfig, candidate::CandidateSpec)
    params = Dict(
        "input" => splitext(basename(cfg.input))[1],
        "noise" => cfg.noise,
        "fitter" => String(cfg.fitter),
        "gradient_mode" => String(cfg.gradient_mode),
        "dn" => candidate.dn,
        "min_valid" => candidate.min_valid,
        "min_valid_spec" => candidate.min_valid_label,
        "rbar2_mode" => candidate.rbar2_label,
        "rbar2_gate" => candidate.rbar2_gate,
        "noise_mode" => String(cfg.noise_mode),
        "noise_sigma_multiplier" => cfg.noise_sigma_multiplier,
    )
    return parameterized_output_dir(cfg.output_dir, params)
end

function candidate_metadata(cfg::RunConfig, candidate::CandidateSpec)
    return Dict(
        "input" => cfg.input,
        "input_absolute_path" => abspath(cfg.input),
        "input_size_bytes" => string(isfile(cfg.input) ? filesize(cfg.input) : -1),
        "input_mtime" => string(isfile(cfg.input) ? stat(cfg.input).mtime : -1.0),
        "noise" => string(cfg.noise),
        "fitter" => string(cfg.fitter),
        "dn" => string(candidate.dn),
        "ker" => string(CANDIDATE_KER),
        "min_valid" => string(candidate.min_valid),
        "min_valid_spec" => string(candidate.min_valid_label),
        "rbar2_weight" => string(candidate.rbar2_weight),
        "rbar2_gate" => string(candidate.rbar2_gate),
        "rbar2_mode" => string(candidate.rbar2_label),
        "noise_mode" => string(cfg.noise_mode),
        "noise_tile_size" => string(cfg.noise_tile_size),
        "noise_sigma_multiplier" => string(cfg.noise_sigma_multiplier),
        "zero_stokes_policy" => string(CANDIDATE_ZERO_STOKES_POLICY),
        "output_support_mode" => string(CANDIDATE_OUTPUT_SUPPORT_MODE),
        "gradient_mode" => string(cfg.gradient_mode),
        "ny_chunks" => string(cfg.ny_chunks),
    )
end

function write_metadata(path::AbstractString, metadata::Dict{String, String})
    mkpath(dirname(path))
    open(path, "w") do io
        for key in sort(collect(keys(metadata)))
            println(io, key, "=", metadata[key])
        end
    end
end

function read_metadata(path::AbstractString)
    out = Dict{String, String}()
    isfile(path) || return out
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line, '='; limit=2)
        length(parts) == 2 || continue
        out[parts[1]] = parts[2]
    end
    return out
end

function metadata_matches(path::AbstractString, expected::Dict{String, String})
    actual = read_metadata(path)
    return all(get(actual, key, nothing) == value for (key, value) in expected)
end

Base.@kwdef struct CandidateRun
    candidate::CandidateSpec
    dir::String
    psi::Matrix{Float64}
    reliability::Matrix{Float64}
    rbar2::Matrix{Float64}
    diagnostics
end

function run_candidate(cfg::RunConfig, candidate::CandidateSpec; force::Bool=false)
    outdir = candidate_dir(cfg, candidate)
    metadata_path = joinpath(outdir, "run_metadata.txt")
    expected = candidate_metadata(cfg, candidate)
    should_force = force || !metadata_matches(metadata_path, expected)

    data, _ = produce_or_load(
        Dict(
            "input" => cfg.input,
            "noise" => cfg.noise,
            "fitter" => String(cfg.fitter),
            "gradient" => String(cfg.gradient_mode),
            "dn" => candidate.dn,
            "ker" => CANDIDATE_KER,
            "min_valid" => candidate.min_valid,
            "min_valid_spec" => candidate.min_valid_label,
            "rbar2_weight" => string(candidate.rbar2_weight),
            "rbar2_gate" => candidate.rbar2_gate,
            "rbar2_mode" => candidate.rbar2_label,
            "noise_mode" => String(cfg.noise_mode),
            "noise_tile_size" => cfg.noise_tile_size,
            "noise_sigma_multiplier" => cfg.noise_sigma_multiplier,
            "zero_stokes" => String(CANDIDATE_ZERO_STOKES_POLICY),
            "output_support" => String(CANDIDATE_OUTPUT_SUPPORT_MODE),
            "gradient" => String(cfg.gradient_mode),
            "ny_chunks" => string(cfg.ny_chunks),
        ),
        outdir;
        filename="candidate",
        prefix="",
        force=should_force,
        verbose=false,
    ) do _
        mkpath(outdir)
        vgt_cfg = VGT.VGTConfig(
            input_file=cfg.input,
            output_dir=outdir,
            keep_intermediate_chunks=cfg.keep_intermediate_chunks,
            noise=cfg.noise,
            dn=candidate.dn,
            ker=CANDIDATE_KER,
            ny_chunks=cfg.ny_chunks,
            min_valid=candidate.min_valid,
            fitter_mode=cfg.fitter,
            rbar2_weight=candidate.rbar2_weight,
            rbar2_gate=candidate.rbar2_gate,
            noise_mode=cfg.noise_mode,
            noise_tile_size=cfg.noise_tile_size,
            noise_sigma_multiplier=cfg.noise_sigma_multiplier,
            zero_stokes_policy=CANDIDATE_ZERO_STOKES_POLICY,
            output_support_mode=CANDIDATE_OUTPUT_SUPPORT_MODE,
            gradient_mode=cfg.gradient_mode,
        )
        run = VGT.run_pipeline(vgt_cfg)
        VGT.write_pipeline_outputs(run)
        write_metadata(metadata_path, expected)
        return Dict(
            "psi" => run.psi,
            "reliability" => Float64.(run.reliability),
            "rbar2" => run.fit_rbar2,
            "diagnostics" => run.diagnostics,
        )
    end

    return CandidateRun(
        candidate=candidate,
        dir=outdir,
        psi=data["psi"],
        reliability=data["reliability"],
        rbar2=data["rbar2"],
        diagnostics=data["diagnostics"],
    )
end

Base.@kwdef struct AbsSummary
    n::Int = 0
    median::Float64 = NaN
    p75::Float64 = NaN
    p90::Float64 = NaN
    p95::Float64 = NaN
    max::Float64 = NaN
    mean::Float64 = NaN
end

Base.@kwdef struct BinnedSummary
    lo::Float64 = NaN
    hi::Float64 = NaN
    n::Int = 0
    median::Float64 = NaN
    p95::Float64 = NaN
end

@inline axial_signed_delta(a, b) = mod(a - b + pi / 2, pi) - pi / 2

function axial_residual_map(a::AbstractMatrix, b::AbstractMatrix)
    size(a) == size(b) || throw(ArgumentError("map sizes differ: $(size(a)) vs $(size(b))"))
    delta = fill(NaN, size(a))
    @inbounds for idx in eachindex(a, b)
        if isfinite(a[idx]) && isfinite(b[idx])
            delta[idx] = axial_signed_delta(a[idx], b[idx])
        end
    end
    return delta
end

function summarize_abs_deg(values)
    clean = Float64[]
    for v in vec(values)
        isfinite(v) || continue
        push!(clean, abs(v) * 180 / pi)
    end
    isempty(clean) && return AbsSummary()
    return AbsSummary(
        n=length(clean),
        median=median(clean),
        p75=quantile(clean, 0.75),
        p90=quantile(clean, 0.90),
        p95=quantile(clean, 0.95),
        max=maximum(clean),
        mean=mean(clean),
    )
end

function summarize_abs_values(values)
    clean = filter(isfinite, collect(vec(values)))
    isempty(clean) && return AbsSummary()
    return AbsSummary(
        n=length(clean),
        median=median(clean),
        p75=quantile(clean, 0.75),
        p90=quantile(clean, 0.90),
        p95=quantile(clean, 0.95),
        max=maximum(clean),
        mean=mean(clean),
    )
end

function binned_abs_summaries(layer::AbstractArray, values::AbstractArray; nbins::Int=5)
    size(layer) == size(values) || throw(ArgumentError("layer and values sizes differ"))
    pairs = [(Float64(layer[i]), Float64(values[i])) for i in eachindex(layer, values)
             if isfinite(layer[i]) && isfinite(values[i])]
    isempty(pairs) && return BinnedSummary[]
    xs = [p[1] for p in pairs]
    ys = [p[2] for p in pairs]
    edges = quantile(xs, range(0, 1; length=nbins + 1))
    out = BinnedSummary[]
    for b in 1:nbins
        lo = edges[b]
        hi = edges[b + 1]
        selected = Float64[]
        for (x, y) in zip(xs, ys)
            in_bin = b == nbins ? (lo <= x <= hi) : (lo <= x < hi)
            in_bin && push!(selected, y)
        end
        isempty(selected) && continue
        push!(out, BinnedSummary(
            lo=lo,
            hi=hi,
            n=length(selected),
            median=median(selected),
            p95=quantile(selected, 0.95),
        ))
    end
    return out
end

Base.@kwdef struct CandidateMetrics
    candidate::CandidateSpec
    finite_psi_fraction::Float64 = NaN
    finite_reliability_fraction::Float64 = NaN
    fit_attempts::Int = 0
    fit_successes::Int = 0
    fit_failures::Int = 0
    valid_windows::Int = 0
    fit_success_rate::Float64 = NaN
    median_rbar2::Float64 = NaN
    median_reliability::Float64 = NaN
    low_confidence_fraction::Float64 = NaN
    adjacent_median_deg::Float64 = NaN
    adjacent_p95_deg::Float64 = NaN
    high_reliability_median_deg::Float64 = NaN
    high_reliability_p95_deg::Float64 = NaN
end

Base.@kwdef struct Recommendation
    candidate::Union{CandidateSpec, Nothing}
    text::String
end

function candidate_score(m::CandidateMetrics)
    coverage = isfinite(m.finite_psi_fraction) ? m.finite_psi_fraction : 0.0
    success = isfinite(m.fit_success_rate) ? m.fit_success_rate : 0.0
    stability = isfinite(m.high_reliability_p95_deg) ? max(0.0, 1.0 - m.high_reliability_p95_deg / 90.0) : 0.0
    confidence = isfinite(m.median_rbar2) ? clamp(m.median_rbar2, 0.0, 1.0) : 0.0
    low_conf_penalty = isfinite(m.low_confidence_fraction) ? m.low_confidence_fraction : 1.0
    return 0.30 * coverage + 0.25 * success + 0.30 * stability + 0.15 * confidence - 0.20 * low_conf_penalty
end

function is_recommendable(m::CandidateMetrics)
    m.finite_psi_fraction >= 0.10 || return false
    m.fit_success_rate >= 0.80 || return false
    isfinite(m.high_reliability_p95_deg) || return false
    m.high_reliability_p95_deg <= 25.0 || return false
    m.low_confidence_fraction <= 0.80 || return false
    return true
end

function rejection_reason(m::CandidateMetrics)
    reasons = String[]
    m.finite_psi_fraction >= 0.10 || push!(reasons, "low finite coverage")
    m.fit_success_rate >= 0.80 || push!(reasons, "low fit success rate")
    isfinite(m.high_reliability_p95_deg) || push!(reasons, "no high-reliability perturbation measurement")
    isfinite(m.high_reliability_p95_deg) && m.high_reliability_p95_deg <= 25.0 ||
        push!(reasons, "unstable high-reliability perturbation p95")
    m.low_confidence_fraction <= 0.80 || push!(reasons, "too many low-confidence windows")
    isempty(reasons) && push!(reasons, "lower score than recommended candidate")
    return join(unique(reasons), ", ")
end

function recommend_candidate(metrics::Vector{CandidateMetrics})
    candidates = filter(is_recommendable, metrics)
    if isempty(candidates)
        return Recommendation(
            candidate=nothing,
            text="No recommendation: no candidate satisfied coverage, fit-success, confidence, and perturbation-stability criteria.",
        )
    end
    best = candidates[argmax(candidate_score.(candidates))]
    lines = String[
        "Recommended for review: dn=$(best.candidate.dn), min_valid=$(best.candidate.min_valid), rbar2_mode=$(best.candidate.rbar2_label)",
        @sprintf("Reason: score=%.3f, finite coverage=%.3f, fit success=%.3f, high-reliability p95=%.3f deg.",
                 candidate_score(best), best.finite_psi_fraction, best.fit_success_rate,
                 best.high_reliability_p95_deg),
        "This is a reporting recommendation only; production config is not changed automatically.",
    ]
    for m in metrics
        m.candidate == best.candidate && continue
        push!(lines,
              "Rejected dn=$(m.candidate.dn), min_valid=$(m.candidate.min_valid), rbar2_mode=$(m.candidate.rbar2_label): $(rejection_reason(m)).")
    end
    return Recommendation(candidate=best.candidate, text=join(lines, "\n"))
end

Base.@kwdef struct PairwiseComparison
    a::CandidateSpec
    b::CandidateSpec
    var"global"::AbsSummary
    high_reliability::AbsSummary
    low_reliability::AbsSummary
end

function finite_fraction(x)
    return count(isfinite, x) / length(x)
end

function finite_median(x)
    vals = filter(isfinite, vec(x))
    isempty(vals) && return NaN
    return median(vals)
end

function extract_candidate_metrics(run::CandidateRun)
    d = run.diagnostics
    success_rate = d.fit_attempts > 0 ? d.fit_successes / d.fit_attempts : NaN
    low_conf_fraction = d.valid_windows > 0 ? d.fit_low_confidence_windows / d.valid_windows : NaN
    return CandidateMetrics(
        candidate=run.candidate,
        finite_psi_fraction=finite_fraction(run.psi),
        finite_reliability_fraction=finite_fraction(run.reliability),
        fit_attempts=d.fit_attempts,
        fit_successes=d.fit_successes,
        fit_failures=d.fit_failures,
        valid_windows=d.valid_windows,
        fit_success_rate=success_rate,
        median_rbar2=finite_median(run.rbar2),
        median_reliability=finite_median(run.reliability),
        low_confidence_fraction=low_conf_fraction,
    )
end

function compare_candidate_pair(a::CandidateRun, b::CandidateRun)
    delta = axial_residual_map(a.psi, b.psi)
    global_summary = summarize_abs_deg(delta)
    rel = b.reliability
    rel_vals = filter(isfinite, vec(rel))
    if isempty(rel_vals)
        return PairwiseComparison(a=a.candidate, b=b.candidate, var"global"=global_summary,
                                  high_reliability=AbsSummary(), low_reliability=AbsSummary())
    end
    high_cut = quantile(rel_vals, 0.80)
    low_cut = quantile(rel_vals, 0.20)
    high = fill(NaN, size(delta))
    low = fill(NaN, size(delta))
    @inbounds for idx in eachindex(delta, rel)
        if isfinite(delta[idx]) && isfinite(rel[idx])
            rel[idx] >= high_cut && (high[idx] = delta[idx])
            rel[idx] <= low_cut && (low[idx] = delta[idx])
        end
    end
    return PairwiseComparison(
        a=a.candidate,
        b=b.candidate,
        var"global"=global_summary,
        high_reliability=summarize_abs_deg(high),
        low_reliability=summarize_abs_deg(low),
    )
end

function _median_or_nan(values)
    clean = filter(isfinite, values)
    isempty(clean) && return NaN
    return median(clean)
end

function attach_adjacency_metrics(metrics::Vector{CandidateMetrics}, pairs::Vector{PairwiseComparison})
    updated = CandidateMetrics[]
    for m in metrics
        adjacent = [p for p in pairs if p.a == m.candidate || p.b == m.candidate]
        if isempty(adjacent)
            push!(updated, m)
            continue
        end
        med = _median_or_nan([p.global.median for p in adjacent])
        p95 = _median_or_nan([p.global.p95 for p in adjacent])
        high_med = _median_or_nan([p.high_reliability.median for p in adjacent])
        high_p95 = _median_or_nan([p.high_reliability.p95 for p in adjacent])
        push!(updated, CandidateMetrics(
            candidate=m.candidate,
            finite_psi_fraction=m.finite_psi_fraction,
            finite_reliability_fraction=m.finite_reliability_fraction,
            fit_attempts=m.fit_attempts,
            fit_successes=m.fit_successes,
            fit_failures=m.fit_failures,
            valid_windows=m.valid_windows,
            fit_success_rate=m.fit_success_rate,
            median_rbar2=m.median_rbar2,
            median_reliability=m.median_reliability,
            low_confidence_fraction=m.low_confidence_fraction,
            adjacent_median_deg=med,
            adjacent_p95_deg=p95,
            high_reliability_median_deg=high_med,
            high_reliability_p95_deg=high_p95,
        ))
    end
    return updated
end

function csv_value(x)
    if x isa Number
        return isfinite(Float64(x)) ? string(x) : ""
    end
    return string(x)
end

function write_candidate_csv(path::AbstractString, metrics::Vector{CandidateMetrics})
    open(path, "w") do io
        println(io, "dn,min_valid,min_valid_spec,rbar2_mode,rbar2_weight,rbar2_gate,finite_psi_fraction,finite_reliability_fraction,fit_attempts,fit_successes,fit_failures,valid_windows,fit_success_rate,median_rbar2,median_reliability,low_confidence_fraction,adjacent_median_deg,adjacent_p95_deg,high_reliability_median_deg,high_reliability_p95_deg")
        for m in metrics
            println(io, join(csv_value.([
                m.candidate.dn, m.candidate.min_valid, m.candidate.min_valid_label,
                m.candidate.rbar2_label, m.candidate.rbar2_weight, m.candidate.rbar2_gate,
                m.finite_psi_fraction,
                m.finite_reliability_fraction, m.fit_attempts, m.fit_successes,
                m.fit_failures, m.valid_windows, m.fit_success_rate, m.median_rbar2,
                m.median_reliability, m.low_confidence_fraction,
                m.adjacent_median_deg, m.adjacent_p95_deg,
                m.high_reliability_median_deg, m.high_reliability_p95_deg,
            ]), ","))
        end
    end
end

function write_pairwise_csv(path::AbstractString, pairs::Vector{PairwiseComparison})
    open(path, "w") do io
        println(io, "a_dn,a_min_valid,a_min_valid_spec,a_rbar2_mode,b_dn,b_min_valid,b_min_valid_spec,b_rbar2_mode,global_n,global_median_deg,global_p95_deg,high_n,high_median_deg,high_p95_deg,low_n,low_median_deg,low_p95_deg")
        for p in pairs
            println(io, join(csv_value.([
                p.a.dn, p.a.min_valid, p.a.min_valid_label,
                p.a.rbar2_label,
                p.b.dn, p.b.min_valid, p.b.min_valid_label,
                p.b.rbar2_label,
                p.global.n, p.global.median, p.global.p95,
                p.high_reliability.n, p.high_reliability.median, p.high_reliability.p95,
                p.low_reliability.n, p.low_reliability.median, p.low_reliability.p95,
            ]), ","))
        end
    end
end

function write_summary_txt(path::AbstractString, metrics::Vector{CandidateMetrics},
                           pairs::Vector{PairwiseComparison}, rec::Recommendation)
    open(path, "w") do io
        println(io, "VGT Parameter Validation Summary")
        println(io, "="^32)
        println(io)
        println(io, rec.text)
        println(io)
        println(io, "Candidate metrics:")
        for m in metrics
            @printf(io,
                    "  dn=%d min_valid=%d spec=%s rbar2_mode=%s coverage=%.4f success=%.4f median_rbar2=%.4f high_rel_p95=%.4f\n",
                    m.candidate.dn, m.candidate.min_valid, m.candidate.min_valid_label,
                    m.candidate.rbar2_label,
                    m.finite_psi_fraction,
                    m.fit_success_rate, m.median_rbar2, m.high_reliability_p95_deg)
        end
        if rec.candidate === nothing
            println(io)
            println(io, "Rejection reasons:")
            for m in metrics
                println(io, "  dn=$(m.candidate.dn) min_valid=$(m.candidate.min_valid) spec=$(m.candidate.min_valid_label) rbar2_mode=$(m.candidate.rbar2_label): $(rejection_reason(m))")
            end
        end
        println(io)
        println(io, "Pairwise residuals:")
        for p in pairs
            @printf(io,
                    "  dn=%d/min=%d/%s/%s vs dn=%d/min=%d/%s/%s global_median=%.4f global_p95=%.4f high_rel_p95=%.4f\n",
                    p.a.dn, p.a.min_valid, p.a.min_valid_label, p.a.rbar2_label,
                    p.b.dn, p.b.min_valid, p.b.min_valid_label, p.b.rbar2_label,
                    p.global.median, p.global.p95, p.high_reliability.p95)
        end
        println(io)
        println(io, "Scientific interpretation:")
        println(io, "This scan selects an internally stable VGT parameter setting for the supplied cube and configuration.")
        println(io, "It does not prove astrophysical correctness. External validation against component-matched starlight polarization or Planck Q/U remains required for foreground-science claims.")
    end
end

function write_reports(output_dir::AbstractString, metrics::Vector{CandidateMetrics},
                       pairs::Vector{PairwiseComparison}, rec::Recommendation)
    mkpath(output_dir)
    write_summary_txt(joinpath(output_dir, "summary.txt"), metrics, pairs, rec)
    write_candidate_csv(joinpath(output_dir, "summary.csv"), metrics)
    write_pairwise_csv(joinpath(output_dir, "pairwise_residuals.csv"), pairs)
end

function write_candidate_summary_jld2(output_dir::AbstractString, cfg::RunConfig,
                                      m::CandidateMetrics, rec::Recommendation)
    params = Dict(
        "input" => splitext(basename(cfg.input))[1],
        "dn" => m.candidate.dn,
        "min_valid" => m.candidate.min_valid,
        "min_valid_spec" => m.candidate.min_valid_label,
        "rbar2_mode" => m.candidate.rbar2_label,
        "rbar2_gate" => m.candidate.rbar2_gate,
        "noise" => cfg.noise,
        "fitter" => String(cfg.fitter),
    )
    path = experiment_summary_path(output_dir, params; prefix="vgt_parameter")
    mkpath(dirname(path))
    is_recommended = rec.candidate !== nothing && rec.candidate == m.candidate
    data = Dict{Symbol, Any}(
        :input => cfg.input,
        :noise => cfg.noise,
        :fitter => String(cfg.fitter),
        :gradient_mode => String(cfg.gradient_mode),
        :dn => m.candidate.dn,
        :min_valid => m.candidate.min_valid,
        :min_valid_spec => m.candidate.min_valid_label,
        :rbar2_mode => m.candidate.rbar2_label,
        :rbar2_weight => m.candidate.rbar2_weight,
        :rbar2_gate => m.candidate.rbar2_gate,
        :finite_psi_fraction => m.finite_psi_fraction,
        :reliability_finite_fraction => m.finite_reliability_fraction,
        :fit_success_rate => m.fit_success_rate,
        :median_rbar2 => m.median_rbar2,
        :median_reliability => m.median_reliability,
        :low_confidence_fraction => m.low_confidence_fraction,
        :adjacent_median_deg => m.adjacent_median_deg,
        :adjacent_p95_deg => m.adjacent_p95_deg,
        :high_reliability_median_deg => m.high_reliability_median_deg,
        :high_reliability_p95_deg => m.high_reliability_p95_deg,
        :recommended => is_recommended,
        :recommendation_text => rec.text,
    )
    tag!(data; gitpath=projectdir(), warn=false, storepatch=false)
    jldsave(path; data...)
    return path
end

function write_summary_jld2(output_dir::AbstractString, cfg::RunConfig,
                            metrics::Vector{CandidateMetrics}, rec::Recommendation)
    paths = String[]
    for m in metrics
        push!(paths, write_candidate_summary_jld2(output_dir, cfg, m, rec))
    end
    return paths
end

function main(args=ARGS)
    opts = parse_args(args)
    cfg = run_config_from_options(opts)
    dn_values = parse_dn_list(opts["dn"])
    specs = parse_min_valid_specs(opts["min-valid"])
    rbar2_specs = parse_rbar2_specs(opts["rbar2-mode"])
    candidates, skipped = expand_candidates(dn_values, specs, rbar2_specs)
    isempty(candidates) && throw(ArgumentError("No valid candidate dn/min_valid combinations"))
    force = parse_bool(opts["force"])

    mkpath(cfg.output_dir)
    println("VGT parameter validation scan")
    println("  input: ", cfg.input)
    println("  output: ", cfg.output_dir)
    println("  fitter: ", cfg.fitter)
    println("  gradient_mode: ", cfg.gradient_mode)
    println("  min_valid selection: reporting-only scan; no production config is changed automatically")
    println("  candidates: ", join(["dn=$(c.dn)/min_valid=$(c.min_valid)/rbar2=$(c.rbar2_label)" for c in candidates], ", "))
    for note in skipped
        println("  skipped: ", note)
    end

    runs = CandidateRun[]
    for candidate in candidates
        println()
        println("Running candidate dn=$(candidate.dn) min_valid=$(candidate.min_valid) rbar2_mode=$(candidate.rbar2_label)")
        push!(runs, run_candidate(cfg, candidate; force=force))
    end

    sorted_runs = sort(runs; by=r -> (r.candidate.rbar2_label, r.candidate.min_valid_label, r.candidate.dn))
    pairs = PairwiseComparison[]
    for key in unique([(r.candidate.rbar2_label, r.candidate.min_valid_label) for r in sorted_runs])
        group = sort(filter(r -> (r.candidate.rbar2_label, r.candidate.min_valid_label) == key, sorted_runs); by=r -> r.candidate.dn)
        for i in 1:(length(group) - 1)
            push!(pairs, compare_candidate_pair(group[i], group[i + 1]))
        end
    end

    metrics = [extract_candidate_metrics(r) for r in sorted_runs]
    metrics = attach_adjacency_metrics(metrics, pairs)
    rec = recommend_candidate(metrics)
    write_reports(cfg.output_dir, metrics, pairs, rec)
    summary_paths = write_summary_jld2(cfg.output_dir, cfg, metrics, rec)
    collected_path = collect_experiment_summaries(cfg.output_dir)

    println()
    println(rec.text)
    println()
    println("Wrote summary: ", joinpath(cfg.output_dir, "summary.txt"))
    println("Wrote metrics: ", joinpath(cfg.output_dir, "summary.csv"))
    println("Wrote residuals: ", joinpath(cfg.output_dir, "pairwise_residuals.csv"))
    println("Wrote summary records: ", length(summary_paths))
    collected_path !== nothing && println("Collected summaries: ", collected_path)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
