using DrWatson
quickactivate(@__DIR__, "VGT")

using Printf
using Statistics
using TOML
using VGT

function default_options()
    return Dict(
        "runs" => nothing,
        "output-dir" => experiment_default_dir("vgt_run_comparison"),
        "write-residual-maps" => "false",
    )
end

function parse_args(args)
    opts = default_options()
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--")
            key = arg[3:end]
            haskey(opts, key) || throw(ArgumentError("Unknown option --$key"))
            i == length(args) && throw(ArgumentError("Missing value for --$key"))
            opts[key] = args[i + 1]
            i += 2
        else
            push!(positional, arg)
            i += 1
        end
    end
    if opts["runs"] === nothing
        opts["runs"] = join(positional, ",")
    elseif !isempty(positional)
        opts["runs"] = string(opts["runs"], ",", join(positional, ","))
    end
    isempty(strip(opts["runs"])) && throw(ArgumentError("Provide run directories via --runs or positional arguments"))
    return opts
end

function parse_bool(raw::AbstractString)
    raw in ("true", "yes", "1") && return true
    raw in ("false", "no", "0") && return false
    throw(ArgumentError("Expected boolean string, got $raw"))
end

function split_runs(raw::AbstractString)
    runs = [strip(part) for part in split(raw, ',') if !isempty(strip(part))]
    length(runs) >= 2 || throw(ArgumentError("At least two run directories are required"))
    return runs
end

function first_match(dir::AbstractString, prefix::AbstractString; exclude_prefixes=String[])
    matches = sort(filter(readdir(dir; join=true)) do name
        base = basename(name)
        startswith(base, prefix) && endswith(base, ".fits") &&
            !any(exclude -> startswith(base, exclude), exclude_prefixes)
    end)
    isempty(matches) && throw(ArgumentError("No $prefix*.fits file found in $dir"))
    return first(matches)
end

function first_match_or_nothing(dir::AbstractString, prefix::AbstractString; exclude_prefixes=String[])
    matches = sort(filter(readdir(dir; join=true)) do name
        base = basename(name)
        startswith(base, prefix) && endswith(base, ".fits") &&
            !any(exclude -> startswith(base, exclude), exclude_prefixes)
    end)
    isempty(matches) && return nothing
    return first(matches)
end

function load_metadata(dir::AbstractString)
    path = joinpath(dir, "run_metadata.toml")
    isfile(path) || return Dict{String, Any}()
    return TOML.parsefile(path)
end

function compact_label(dir::AbstractString)
    base = basename(normpath(dir))
    m = match(r"cfg=([^_]+)", base)
    m === nothing || return string("cfg_", m.captures[1])
    return replace(base, r"[^A-Za-z0-9_.-]" => "_")[1:min(end, 48)]
end

function clean_values(x)
    vals = Float64[]
    for v in vec(x)
        isfinite(v) && push!(vals, Float64(v))
    end
    return vals
end

function finite_fraction(x)
    return count(isfinite, x) / length(x)
end

function optional_fits_image(path, dims::Tuple)
    path === nothing && return fill(NaN, dims)
    image = Float64.(VGT.read_fits_image(path))
    size(image) == dims || throw(ArgumentError("Optional FITS product $path has size $(size(image)); expected $dims"))
    return image
end

function qtile(vals, q)
    isempty(vals) && return NaN
    return quantile(vals, q)
end

@inline axial_signed_delta(a, b) = mod(a - b + pi / 2, pi) - pi / 2

function axial_residual(a, b)
    size(a) == size(b) || throw(ArgumentError("Map sizes differ: $(size(a)) vs $(size(b))"))
    out = fill(NaN, size(a))
    @inbounds for idx in eachindex(a, b)
        if isfinite(a[idx]) && isfinite(b[idx])
            out[idx] = axial_signed_delta(Float64(a[idx]), Float64(b[idx]))
        end
    end
    return out
end

function summarize_abs_deg(values)
    vals = [abs(Float64(v)) * 180 / pi for v in vec(values) if isfinite(v)]
    return Dict(
        "n" => length(vals),
        "median_deg" => qtile(vals, 0.50),
        "p75_deg" => qtile(vals, 0.75),
        "p90_deg" => qtile(vals, 0.90),
        "p95_deg" => qtile(vals, 0.95),
        "mean_deg" => isempty(vals) ? NaN : mean(vals),
    )
end

function median_on_mask(layer, mask)
    size(layer) == size(mask) || throw(ArgumentError("Layer size $(size(layer)) does not match mask size $(size(mask))"))
    vals = Float64[]
    @inbounds for idx in eachindex(layer, mask)
        if mask[idx] && isfinite(layer[idx])
            push!(vals, Float64(layer[idx]))
        end
    end
    isempty(vals) && return NaN
    return median(vals)
end

function finite_support_masks(baseline_psi, candidate_psi)
    size(baseline_psi) == size(candidate_psi) ||
        throw(ArgumentError("Map sizes differ: $(size(baseline_psi)) vs $(size(candidate_psi))"))
    baseline_finite = isfinite.(baseline_psi)
    candidate_finite = isfinite.(candidate_psi)
    return (;
        common=baseline_finite .& candidate_finite,
        gained=.!baseline_finite .& candidate_finite,
        lost=baseline_finite .& .!candidate_finite,
    )
end

function partition_count(mask)
    return count(identity, mask)
end

function partition_fraction(mask)
    return partition_count(mask) / length(mask)
end

function support_partition_row(baseline, candidate)
    masks = finite_support_masks(baseline.psi, candidate.psi)
    common_delta = axial_residual(baseline.psi, candidate.psi)
    common_stats = summarize_abs_deg(common_delta)
    return Dict{String, Any}(
        "baseline" => baseline.label,
        "candidate" => candidate.label,
        "baseline_dir" => baseline.dir,
        "candidate_dir" => candidate.dir,
        "common_fraction" => partition_fraction(masks.common),
        "gained_fraction" => partition_fraction(masks.gained),
        "lost_fraction" => partition_fraction(masks.lost),
        "common_n" => partition_count(masks.common),
        "gained_n" => partition_count(masks.gained),
        "lost_n" => partition_count(masks.lost),
        "common_delta_median_deg" => common_stats["median_deg"],
        "common_delta_p90_deg" => common_stats["p90_deg"],
        "baseline_common_support_median" => median_on_mask(baseline.support, masks.common),
        "candidate_common_support_median" => median_on_mask(candidate.support, masks.common),
        "baseline_gained_support_median" => median_on_mask(baseline.support, masks.gained),
        "candidate_gained_support_median" => median_on_mask(candidate.support, masks.gained),
        "baseline_lost_support_median" => median_on_mask(baseline.support, masks.lost),
        "candidate_lost_support_median" => median_on_mask(candidate.support, masks.lost),
        "baseline_common_reliability_median" => median_on_mask(baseline.reliability, masks.common),
        "candidate_common_reliability_median" => median_on_mask(candidate.reliability, masks.common),
        "baseline_gained_reliability_median" => median_on_mask(baseline.reliability, masks.gained),
        "candidate_gained_reliability_median" => median_on_mask(candidate.reliability, masks.gained),
        "baseline_lost_reliability_median" => median_on_mask(baseline.reliability, masks.lost),
        "candidate_lost_reliability_median" => median_on_mask(candidate.reliability, masks.lost),
        "baseline_common_rbar2_median" => median_on_mask(baseline.rbar2, masks.common),
        "candidate_common_rbar2_median" => median_on_mask(candidate.rbar2, masks.common),
        "baseline_gained_rbar2_median" => median_on_mask(baseline.rbar2, masks.gained),
        "candidate_gained_rbar2_median" => median_on_mask(candidate.rbar2, masks.gained),
        "baseline_lost_rbar2_median" => median_on_mask(baseline.rbar2, masks.lost),
        "candidate_lost_rbar2_median" => median_on_mask(candidate.rbar2, masks.lost),
    )
end

function support_partition_rows(runs)
    length(runs) >= 2 || return Dict{String, Any}[]
    baseline = first(runs)
    return [support_partition_row(baseline, candidate) for candidate in runs[2:end]]
end

function write_csv(path::AbstractString, rows, columns)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(io, join([string(get(row, col, "")) for col in columns], ","))
        end
    end
end

function read_run(dir::AbstractString)
    psi_path = first_match(dir, "psi_"; exclude_prefixes=["psi_reliability_", "psi_fit_rbar2_", "psi_support_", "psi_plot_"])
    reliability_path = first_match(dir, "psi_reliability_")
    rbar2_path = first_match(dir, "psi_fit_rbar2_")
    support_path = first_match_or_nothing(dir, "psi_support_")
    psi = Float64.(VGT.read_fits_image(psi_path))
    return (;
        dir,
        label=compact_label(dir),
        metadata=load_metadata(dir),
        psi_path,
        reliability_path,
        rbar2_path,
        support_path,
        psi,
        reliability=Float64.(VGT.read_fits_image(reliability_path)),
        rbar2=Float64.(VGT.read_fits_image(rbar2_path)),
        support=optional_fits_image(support_path, size(psi)),
    )
end

function run_summary(run)
    reliability = clean_values(run.reliability)
    rbar2 = clean_values(run.rbar2)
    diagnostics = get(run.metadata, "diagnostics", Dict{String, Any}())
    attempts = get(diagnostics, "fit_attempts", 0)
    low_conf = get(diagnostics, "fit_low_confidence_windows", 0)
    return Dict(
        "label" => run.label,
        "dir" => run.dir,
        "finite_psi_fraction" => finite_fraction(run.psi),
        "finite_reliability_fraction" => finite_fraction(run.reliability),
        "median_reliability" => qtile(reliability, 0.50),
        "p95_reliability" => qtile(reliability, 0.95),
        "median_rbar2" => qtile(rbar2, 0.50),
        "p95_rbar2" => qtile(rbar2, 0.95),
        "fit_attempts" => attempts,
        "fit_low_confidence_windows" => low_conf,
        "low_confidence_fraction" => attempts == 0 ? NaN : low_conf / attempts,
    )
end

function compare_runs(opts)
    runs = [read_run(dir) for dir in split_runs(opts["runs"])]
    output_dir = opts["output-dir"]
    mkpath(output_dir)
    summary_rows = [run_summary(run) for run in runs]
    write_csv(
        joinpath(output_dir, "summary.csv"),
        summary_rows,
        ["label", "dir", "finite_psi_fraction", "finite_reliability_fraction",
         "median_reliability", "p95_reliability", "median_rbar2", "p95_rbar2",
         "fit_attempts", "fit_low_confidence_windows", "low_confidence_fraction"],
    )

    pair_rows = Dict{String, Any}[]
    write_maps = parse_bool(opts["write-residual-maps"])
    for i in 1:length(runs)-1, j in i+1:length(runs)
        residual = axial_residual(runs[i].psi, runs[j].psi)
        stats = summarize_abs_deg(residual)
        row = Dict{String, Any}(
            "a" => runs[i].label,
            "b" => runs[j].label,
            "common_finite_fraction" => count(isfinite, residual) / length(residual),
        )
        merge!(row, stats)
        push!(pair_rows, row)
        if write_maps
            map_path = joinpath(output_dir, string("residual_", runs[i].label, "__", runs[j].label, ".fits"))
            VGT.write_fits(map_path, residual)
        end
    end
    write_csv(
        joinpath(output_dir, "pairwise_residuals.csv"),
        pair_rows,
        ["a", "b", "common_finite_fraction", "n", "median_deg", "p75_deg", "p90_deg", "p95_deg", "mean_deg"],
    )
    partition_rows = support_partition_rows(runs)
    write_csv(
        joinpath(output_dir, "support_partitions.csv"),
        partition_rows,
        ["baseline", "candidate",
         "baseline_dir", "candidate_dir",
         "common_fraction", "gained_fraction", "lost_fraction",
         "common_n", "gained_n", "lost_n",
         "common_delta_median_deg", "common_delta_p90_deg",
         "baseline_common_support_median", "candidate_common_support_median",
         "baseline_gained_support_median", "candidate_gained_support_median",
         "baseline_lost_support_median", "candidate_lost_support_median",
         "baseline_common_reliability_median", "candidate_common_reliability_median",
         "baseline_gained_reliability_median", "candidate_gained_reliability_median",
         "baseline_lost_reliability_median", "candidate_lost_reliability_median",
         "baseline_common_rbar2_median", "candidate_common_rbar2_median",
         "baseline_gained_rbar2_median", "candidate_gained_rbar2_median",
         "baseline_lost_rbar2_median", "candidate_lost_rbar2_median"],
    )
    @printf "Wrote %s, %s, and %s\n" joinpath(output_dir, "summary.csv") joinpath(output_dir, "pairwise_residuals.csv") joinpath(output_dir, "support_partitions.csv")
    return output_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        compare_runs(parse_args(ARGS))
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end
