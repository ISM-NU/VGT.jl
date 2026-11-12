using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using JLD2
using Printf
using Profile
using Statistics

using VGT: VGT_QU_mv_2d, VGT_QU_mv_2d_profiled,
                 collect_experiment_summaries, experiment_default_dir, experiment_summary_path,
                 compute_stokes_intensity, imfilter_gaussian, load_config

function parse_args(args)
    opts = Dict(
        "input" => nothing,
        "output-dir" => experiment_default_dir("profile_pipeline"),
        "noise" => nothing,
        "dn" => nothing,
        "ker" => nothing,
        "min-valid" => nothing,
        "x1" => nothing,
        "x2" => nothing,
        "y1" => nothing,
        "y2" => nothing,
        "crop-size" => "256",
        "autocrop" => "densest",
        "ni" => "1",
        "nf" => nothing,
        "nchannels" => "32",
        "profile-delay" => "0.001",
        "warmup" => "true",
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unexpected argument: $arg")
        key = arg[3:end]
        haskey(opts, key) || error("Unknown option --$key")
        i == length(args) && error("Missing value for --$key")
        opts[key] = args[i + 1]
        i += 2
    end
    return opts
end

function load_config_defaults()
    defaults = Dict{String, Any}()
    cfg = load_config(projectdir("config.toml"))
    defaults["input"] = cfg.input_file
    defaults["noise"] = cfg.noise
    defaults["dn"] = cfg.dn
    defaults["ker"] = cfg.ker
    defaults["min-valid"] = cfg.min_valid
    return defaults
end

function resolve_opt(opts, defaults, key::String)
    value = opts[key]
    value === nothing ? defaults[key] : value
end

parse_int_opt(value) = value isa Integer ? Int(value) : parse(Int, string(value))
parse_float_opt(value) = value isa AbstractFloat ? Float64(value) : parse(Float64, string(value))

function resolve_channel_range(ni_req::Int, nf_req, nchannels::Int, nv::Int)
    nv >= 1 || throw(ArgumentError("Input cube must have at least one channel; got nv=$nv"))
    nchannels >= 1 || throw(ArgumentError("nchannels must be >= 1; got $nchannels"))

    if nf_req === nothing
        ni = clamp(ni_req, 1, nv)
        if ni + nchannels - 1 > nv
            ni = max(1, nv - nchannels + 1)
        end
        nf = min(nv, ni + nchannels - 1)
        return ni, nf
    end

    ni = clamp(ni_req, 1, nv)
    nf = clamp(Int(nf_req), ni, nv)
    return ni, nf
end

function load_subcube(path::String, x1::Int, x2::Int, y1::Int, y2::Int, ni::Int, nf::Int)
    f = FITS(path)
    cube = try
        read(f[1], x1:x2, y1:y2, ni:nf)
    finally
        close(f)
    end
    return Array{Float64}(cube)
end

function find_dense_patch(path::String, crop_size::Int, ni::Int, nf::Int, noise::Float64)
    f = FITS(path)
    cube_hdu = f[1]
    dims = size(cube_hdu)
    nx, ny = dims[1], dims[2]
    window = min(crop_size, nx, ny)
    half = max(1, div(window, 2))
    step = max(window, 128)
    thresh = noise
    sample_chans = unique(Int[
        ni,
        clamp(ni + div(max(nf - ni, 0), 2), ni, nf),
        nf,
    ])

    best_sum = -Inf
    best_cx = min(nx, half)
    best_cy = min(ny, half)
    try
        for cx in half:step:(nx - half + 1), cy in half:step:(ny - half + 1)
            x1 = max(1, cx - half + 1)
            x2 = min(nx, x1 + window - 1)
            y1 = max(1, cy - half + 1)
            y2 = min(ny, y1 + window - 1)
            local_sum = 0.0
            for chan in sample_chans
                patch = read(cube_hdu, x1:x2, y1:y2, chan:chan)
                @inbounds for v in patch
                    if isfinite(v) && v > thresh
                        local_sum += v
                    end
                end
            end
            if local_sum > best_sum
                best_sum = local_sum
                best_cx = cx
                best_cy = cy
            end
        end
    finally
        close(f)
    end

    x1 = max(1, best_cx - half + 1)
    y1 = max(1, best_cy - half + 1)
    x2 = min(nx, x1 + window - 1)
    y2 = min(ny, y1 + window - 1)
    return x1, x2, y1, y2, best_sum
end

function maybe_select_crop(opts, input_path::String, dims, ni::Int, nf::Int, noise::Float64)
    if all(!isnothing, (opts["x1"], opts["x2"], opts["y1"], opts["y2"]))
        return parse(Int, opts["x1"]), parse(Int, opts["x2"]),
               parse(Int, opts["y1"]), parse(Int, opts["y2"]), nothing
    end

    crop_size = parse(Int, opts["crop-size"])
    if opts["autocrop"] == "densest"
        return find_dense_patch(input_path, crop_size, ni, nf, noise)
    end

    nx, ny = dims[1], dims[2]
    x1 = 1
    y1 = 1
    x2 = min(nx, crop_size)
    y2 = min(ny, crop_size)
    return x1, x2, y1, y2, nothing
end

function timed_result(f)
    elapsed = @elapsed result = f()
    return result, elapsed
end

function safe_ratio(num, den)
    den == 0 ? 0.0 : num / den
end

function summarize_profile(ns_values::Vector{Int})
    isempty(ns_values) && return (total_s=0.0, mean_s=0.0, median_s=0.0, max_s=0.0)
    secs = ns_values ./ 1e9
    return (
        total_s = sum(secs),
        mean_s = mean(secs),
        median_s = median(secs),
        max_s = maximum(secs),
    )
end

function format_seconds(io, label, value)
    @printf(io, "  %-24s %.6f s\n", label, value)
end

function build_report(io, metadata, timings, vgt_profile)
    println(io, "Pipeline Profiling Report")
    println(io, "="^72)
    println(io, "Input file: $(metadata[:input])")
    println(io, "Crop: x=$(metadata[:x1]):$(metadata[:x2]), y=$(metadata[:y1]):$(metadata[:y2])")
    if haskey(metadata, :crop_signal) && metadata[:crop_signal] !== nothing
        @printf(io, "Autocrop signal sum (>noise): %.4f\n", metadata[:crop_signal])
    end
    println(io, "Channels: $(metadata[:ni]):$(metadata[:nf]) ($(metadata[:nchannels]) total)")
    println(io, "Parameters: noise=$(metadata[:noise]) dn=$(metadata[:dn]) ker=$(metadata[:ker]) min_valid=$(metadata[:min_valid])")
    println(io)
    println(io, "Stage timings")
    println(io, "-"^72)
    format_seconds(io, "read_subcube", timings[:read_s])
    format_seconds(io, "vgt_total", timings[:vgt_s])
    format_seconds(io, "smooth_q", timings[:smooth_q_s])
    format_seconds(io, "smooth_u", timings[:smooth_u_s])
    format_seconds(io, "smooth_i", timings[:smooth_i_s])
    format_seconds(io, "psi_compute", timings[:psi_s])
    format_seconds(io, "reliability_compute", timings[:reliability_s])
    format_seconds(io, "finalize_total", timings[:finalize_total_s])
    format_seconds(io, "total_measured", timings[:total_s])
    println(io)

    grad = summarize_profile(vgt_profile[:gradient_ns])
    window = summarize_profile(vgt_profile[:window_ns])
    fit_s = vgt_profile[:fit_ns] / 1e9
    println(io, "VGT core breakdown")
    println(io, "-"^72)
    format_seconds(io, "gradient_total", grad.total_s)
    format_seconds(io, "window_total", window.total_s)
    format_seconds(io, "fit_total_inside_windows", fit_s)
    @printf(io, "  %-24s %d\n", "fit_attempts", vgt_profile[:fit_attempts])
    @printf(io, "  %-24s %d\n", "fit_successes", vgt_profile[:fit_successes])
    @printf(io, "  %-24s %d\n", "fit_failures", get(vgt_profile, :fit_failures, 0))
    @printf(io, "  %-24s %d\n", "valid_windows", vgt_profile[:valid_windows])
    @printf(io, "  %-24s %.2f%%\n", "fit_success_rate", 100 * safe_ratio(vgt_profile[:fit_successes], max(vgt_profile[:fit_attempts], 1)))
    @printf(io, "  %-24s %.2f%%\n", "fit_share_of_vgt", 100 * safe_ratio(fit_s, timings[:vgt_s]))
    println(io)

    println(io, "Per-channel means")
    println(io, "-"^72)
    format_seconds(io, "gradient_mean", grad.mean_s)
    format_seconds(io, "window_mean", window.mean_s)
    format_seconds(io, "accumulation_mean", accum.mean_s)
    format_seconds(io, "window_median", window.median_s)
    format_seconds(io, "window_max", window.max_s)
end

function warmup_kernels(cube, dn::Int, noise::Float64, ker::Int, min_valid::Int)
    warm_nx = min(size(cube, 1), max(dn, 32))
    warm_ny = min(size(cube, 2), max(dn, 32))
    warm_nv = min(size(cube, 3), 2)
    warm_cube = @view cube[1:warm_nx, 1:warm_ny, 1:warm_nv]
    q_accum, u_accum = VGT_QU_mv_2d(warm_cube, dn, 1, size(warm_cube, 3), noise; min_valid=min_valid)
    intensity_accum = compute_stokes_intensity(warm_cube)
    imfilter_gaussian(q_accum, ker)
    imfilter_gaussian(u_accum, ker)
    imfilter_gaussian(intensity_accum, ker)
    return nothing
end

function write_summary_jld2(output_dir::AbstractString, metadata, timings, vgt_profile)
    params = Dict(
        "input" => splitext(basename(metadata[:input]))[1],
        "dn" => metadata[:dn],
        "noise" => metadata[:noise],
        "crop" => "$(metadata[:x1])-$(metadata[:x2])_$(metadata[:y1])-$(metadata[:y2])",
        "channels" => "$(metadata[:ni])-$(metadata[:nf])",
    )
    path = experiment_summary_path(output_dir, params; prefix="profile_pipeline")
    mkpath(dirname(path))
    data = Dict{Symbol, Any}(
        :input => metadata[:input],
        :x1 => metadata[:x1],
        :x2 => metadata[:x2],
        :y1 => metadata[:y1],
        :y2 => metadata[:y2],
        :ni => metadata[:ni],
        :nf => metadata[:nf],
        :nchannels => metadata[:nchannels],
        :noise => metadata[:noise],
        :dn => metadata[:dn],
        :ker => metadata[:ker],
        :min_valid => metadata[:min_valid],
        :crop_signal => metadata[:crop_signal] === nothing ? NaN : metadata[:crop_signal],
        :read_s => timings[:read_s],
        :vgt_s => timings[:vgt_s],
        :smooth_q_s => timings[:smooth_q_s],
        :smooth_u_s => timings[:smooth_u_s],
        :smooth_i_s => timings[:smooth_i_s],
        :psi_s => timings[:psi_s],
        :reliability_s => timings[:reliability_s],
        :finalize_total_s => timings[:finalize_total_s],
        :total_s => timings[:total_s],
        :fit_attempts => vgt_profile[:fit_attempts],
        :fit_successes => vgt_profile[:fit_successes],
        :fit_failures => get(vgt_profile, :fit_failures, 0),
        :valid_windows => vgt_profile[:valid_windows],
        :fit_s => vgt_profile[:fit_ns] / 1e9,
    )
    tag!(data; gitpath=projectdir(), warn=false, storepatch=false)
    jldsave(path; data...)
    return path
end

function main(args)
    opts = parse_args(args)
    defaults = load_config_defaults()

    input_path = abspath(string(resolve_opt(opts, defaults, "input")))
    output_dir = abspath(opts["output-dir"])
    mkpath(output_dir)

    noise = parse_float_opt(resolve_opt(opts, defaults, "noise"))
    dn = parse_int_opt(resolve_opt(opts, defaults, "dn"))
    ker = parse_int_opt(resolve_opt(opts, defaults, "ker"))
    min_valid = parse_int_opt(resolve_opt(opts, defaults, "min-valid"))

    ni_req = parse(Int, opts["ni"])
    nchannels = parse(Int, opts["nchannels"])

    f = FITS(input_path)
    dims = try
        size(f[1])
    finally
        close(f)
    end
    nv = dims[3]
    nf_req = opts["nf"] === nothing ? nothing : parse(Int, opts["nf"])
    ni, nf = resolve_channel_range(ni_req, nf_req, nchannels, nv)
    x1, x2, y1, y2, crop_signal = maybe_select_crop(opts, input_path, dims, ni, nf, noise)

    _, read_s = timed_result(() -> nothing)
    cube, read_s = timed_result(() -> load_subcube(input_path, x1, x2, y1, y2, ni, nf))

    if lowercase(opts["warmup"]) == "true"
        warmup_kernels(cube, dn, noise, ker, min_valid)
    end

    profile_out = joinpath(output_dir, "profile_vgt.txt")
    q_accum = nothing
    u_accum = nothing
    vgt_profile = Dict(
        :gradient_ns => Int[],
        :window_ns => Int[],
        :fit_ns => 0,
        :fit_attempts => 0,
        :fit_successes => 0,
        :fit_failures => 0,
        :valid_windows => 0,
    )
    Profile.clear()
    Profile.init(delay=parse(Float64, opts["profile-delay"]))
    Profile.clear()
    vgt_s = @elapsed begin
        @profile begin
            q_accum, u_accum, vgt_profile = VGT_QU_mv_2d_profiled(cube, dn, 1, size(cube, 3), noise;
                                                                  min_valid=min_valid)
        end
    end
    open(profile_out, "w") do io
        Profile.print(io; format=:flat, sortedby=:count)
    end
    intensity_accum = compute_stokes_intensity(cube)
    q_smooth, smooth_q_s = timed_result(() -> imfilter_gaussian(q_accum, ker))
    u_smooth, smooth_u_s = timed_result(() -> imfilter_gaussian(u_accum, ker))
    intensity_smooth, smooth_i_s = timed_result(() -> imfilter_gaussian(intensity_accum, ker))
    psi_s = @elapsed begin
        q_smooth[intensity_accum .== 0] .= NaN
        u_smooth[intensity_accum .== 0] .= NaN
        0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)
    end
    reliability_s = @elapsed begin
        sqrt.(q_smooth.^2 .+ u_smooth.^2) ./ intensity_smooth
    end
    finalize_s = smooth_q_s + smooth_u_s + smooth_i_s + psi_s + reliability_s

    timings = Dict(
        :read_s => read_s,
        :vgt_s => vgt_s,
        :smooth_q_s => smooth_q_s,
        :smooth_u_s => smooth_u_s,
        :smooth_i_s => smooth_i_s,
        :psi_s => psi_s,
        :reliability_s => reliability_s,
        :finalize_total_s => finalize_s,
    )
    timings[:total_s] = read_s + vgt_s + finalize_s

    report_out = joinpath(output_dir, "profile_report.txt")
    metadata = Dict(
        :input => input_path,
        :x1 => x1,
        :x2 => x2,
        :y1 => y1,
        :y2 => y2,
        :ni => ni,
        :nf => nf,
        :nchannels => nf - ni + 1,
        :noise => noise,
        :dn => dn,
        :ker => ker,
        :min_valid => min_valid,
        :crop_signal => crop_signal,
    )
    open(report_out, "w") do io
        build_report(io, metadata, timings, vgt_profile)
    end
    summary_path = write_summary_jld2(output_dir, metadata, timings, vgt_profile)
    collected_path = collect_experiment_summaries(output_dir)

    println("Profile report: $report_out")
    println("Julia Profile dump: $profile_out")
    println("Summary record: $summary_path")
    collected_path !== nothing && println("Collected summaries: $collected_path")
    println()
    open(report_out, "r") do io
        print(read(io, String))
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
