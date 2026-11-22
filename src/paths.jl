const DEFAULT_RUN_FAMILY = "vgt"

function project_config_path(filename::AbstractString="config.toml")
    return projectdir(filename)
end

function rawdatadir(parts::AbstractString...)
    return datadir("exp_raw", parts...)
end

function processeddatadir(parts::AbstractString...)
    return datadir("exp_pro", parts...)
end

function simulationdir(parts::AbstractString...)
    return datadir("sims", parts...)
end

function experiment_default_dir(kind::AbstractString)
    return simulationdir(kind)
end

function _stringify_path_param(value)
    value isa Symbol && return String(value)
    value isa AbstractString && return String(value)
    value isa Bool && return value ? "true" : "false"
    return value
end

function run_parameters(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY)
    params = Dict(
        "label" => label,
        "input" => splitext(basename(cfg.input_file))[1],
        "noise" => cfg.noise,
        "dn" => cfg.dn,
        "ker" => cfg.ker,
        "min_valid" => cfg.min_valid,
        "fitter" => _stringify_path_param(cfg.fitter_mode),
        "rbar2_weight" => _stringify_path_param(cfg.rbar2_weight),
        "rbar2_gate" => cfg.rbar2_gate,
        "noise_mode" => _stringify_path_param(cfg.noise_mode),
        "noise_tile_size" => cfg.noise_tile_size,
        "noise_sigma_multiplier" => cfg.noise_sigma_multiplier,
        "zero_stokes" => _stringify_path_param(cfg.zero_stokes_policy),
        "output_support" => _stringify_path_param(cfg.output_support_mode),
        "gradient" => _stringify_path_param(cfg.gradient_mode),
    )
    # Only a sky frame enters the run name, so pixel-frame run directories keep their names.
    cfg.angle_frame === :pixel || (params["angle_frame"] = String(cfg.angle_frame))
    cfg.output_bin == 1 || (params["output_bin"] = cfg.output_bin)
    # Only when on, so default run names and hashes are unchanged.
    if cfg.gradient_whiten
        params["gradient_whiten"] = "true"
        params["gradient_noise_cov"] = join(string.(cfg.gradient_noise_cov), "_")
        isempty(cfg.gradient_noise_cov_file) || (params["gradient_noise_cov_file"] = cfg.gradient_noise_cov_file)
    end
    cfg.edge_mask && (params["edge_mask"] = "true")
    isempty(cfg.velocity_range_kms) || (params["vrange"] = join(string.(cfg.velocity_range_kms), "_"))
    return params
end

function run_parameter_hash(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY)
    io = IOBuffer()
    TOML.print(io, Dict("run_parameters" => run_parameters(cfg; label=label)); sorted=true)
    return bytes2hex(sha1(String(take!(io))))[1:12]
end

function run_name_parameters(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY)
    params = Dict(
        "input" => splitext(basename(cfg.input_file))[1],
        "noise" => cfg.noise,
        "dn" => cfg.dn,
        "ker" => cfg.ker,
        "min_valid" => cfg.min_valid,
        "fitter" => _stringify_path_param(cfg.fitter_mode),
        "noise_mode" => _stringify_path_param(cfg.noise_mode),
        "noise_sigma_multiplier" => cfg.noise_sigma_multiplier,
        "cfg" => run_parameter_hash(cfg; label=label),
    )
    return params
end

function run_savename(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY,
                      ext::AbstractString="")
    params = run_name_parameters(cfg; label=label)
    return isempty(ext) ? savename(params) : savename(params, ext)
end

function default_output_dir(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY)
    return simulationdir(label, run_savename(cfg; label=label))
end

function parameterized_output_dir(root::AbstractString, params; ext::AbstractString="")
    name = isempty(ext) ? savename(params) : savename(params, ext)
    return joinpath(root, name)
end

function experiment_output_dir(kind::AbstractString, params; ext::AbstractString="")
    return parameterized_output_dir(simulationdir(kind), params; ext=ext)
end

function experiment_summary_dir(output_dir::AbstractString)
    return joinpath(output_dir, "summaries")
end

function experiment_summary_path(output_dir::AbstractString, params; prefix::AbstractString="summary",
                                 ext::AbstractString="jld2")
    filename = string(prefix, "_", savename(params, ext))
    return joinpath(experiment_summary_dir(output_dir), filename)
end

function collect_experiment_summaries(output_dir::AbstractString)
    summaries = experiment_summary_dir(output_dir)
    isdir(summaries) || return nothing
    path = joinpath(output_dir, "collected_results.jld2")

    if isdefined(DrWatson, :collect_results)
        table = getproperty(DrWatson, :collect_results)(summaries)
        jldsave(path; results=table)
        return path
    end

    records = Dict{String, Any}[]
    for file in sort(readdir(summaries; join=true))
        endswith(lowercase(file), ".jld2") || continue
        loaded = JLD2.load(file)
        record = Dict{String, Any}("path" => file)
        for (key, value) in loaded
            record[string(key)] = value
        end
        push!(records, record)
    end
    jldsave(path; results=records)
    return path
end
