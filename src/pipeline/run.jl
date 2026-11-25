# Pipeline orchestration.
#
# This file connects configuration, chunk processing, stitching, and final
# product formation. The numerical science lives in src/science/.

Base.@kwdef struct PipelineRun
    # In-memory result of one pipeline run before FITS files are written.
    config::VGTConfig
    # psi, reliability, fit_rbar2, fit_support and header are on the output
    # grid (binned when config.output_bin > 1); Q, U, I, W are the raw stitched
    # accumulators on the input grid, before smoothing. I is the channel-summed
    # intensity (sky mask); W is the summed Stokes weight of the fitted windows
    # (reliability denominator).
    psi::Matrix{Float64}
    reliability::Matrix{Float32}
    Q::Matrix{Float64}
    U::Matrix{Float64}
    I::Matrix{Float64}
    W::Matrix{Float64}
    fit_rbar2::Matrix{Float64}
    fit_support::Matrix{Float64}
    header
    diagnostics::VGTDiagnostics
    # Input channels used (all, or the velocity_range_kms window).
    channels::UnitRange{Int} = 1:0
    performance::NamedTuple = (;
        julia_threads=Threads.nthreads(),
        ny_chunks=-1,
        elapsed_seconds=NaN,
        started_at_utc="",
        finished_at_utc="",
    )
end

function resolved_output_dir(cfg::VGTConfig; label::AbstractString=DEFAULT_RUN_FAMILY)
    # Empty/default output_dir values are expanded to a DrWatson run directory
    # keyed by the important science parameters.
    if isempty(cfg.output_dir) || cfg.output_dir == "subcube"
        return default_output_dir(cfg; label=label)
    end
    return cfg.output_dir
end

function _with_output_dir(cfg::VGTConfig, output_dir::AbstractString)
    return config_with(cfg; output_dir=String(output_dir))
end

# Header of the cube cut to `channels`: the spectral axis starts at the first
# selected channel, so products carry the window a cut cube would carry.
function window_header(header, channels::UnitRange{Int}, nv_input::Int)
    channels == 1:nv_input && return header
    h = deepcopy(header)
    haskey(h, "CRPIX3") && (h["CRPIX3"] = Float64(h["CRPIX3"]) - (first(channels) - 1))
    haskey(h, "NAXIS3") && (h["NAXIS3"] = length(channels))
    return h
end

function run_pipeline(cfg::VGTConfig)
    # Resolve execution settings, process every chunk, stitch chunk products,
    # then convert accumulated Q/U/W/I into final psi and reliability maps.
    started_at = Dates.now(Dates.UTC)
    start_ns = time_ns()
    output_dir = resolved_output_dir(cfg)
    output_dir_created = !isdir(output_dir)
    cfg = _with_output_dir(cfg, output_dir)
    mkpath(cfg.output_dir)
    resolved_cfg = cfg

    try
        isfile(cfg.input_file) || error("$(cfg.input_file) not found. Current directory: $(pwd())")
        nx, ny, nv_input, header = read_cube_metadata(cfg.input_file)
        channels = selected_channels(cfg, header, nv_input)
        nv = length(channels)
        header = window_header(header, channels, nv_input)
        resolved_cfg = resolve_execution_config(cfg, nx, ny, nv)
        # Fail before the expensive part if a sky frame is asked for a cube without one.
        cfg.angle_frame === :pixel || (_celestial_wcs(header); _map_frame(header))
        @info "Resolved execution settings" julia_threads=Threads.nthreads() chunks=resolved_cfg.ny_chunks memory_budget_gb=round(memory_budget_bytes() / 1024^3; digits=1)
        isempty(cfg.weight_file) || (read_cube_metadata(cfg.weight_file)[1:3] == (nx, ny, nv_input) ||
            error("weight_file $(cfg.weight_file) does not match the size of $(cfg.input_file)"))
        noise_model = build_noise_model(
            weight_source(resolved_cfg),
            resolved_cfg.noise,
            resolved_cfg.noise_mode;
            tile_size=resolved_cfg.noise_tile_size,
            channels=channels,
        )
        nv == nv_input || @info "Velocity window" velocity_range_kms=cfg.velocity_range_kms channels=channels
        @info "Resolved noise model" requested=resolved_cfg.noise_mode resolved=noise_model.mode scalar_noise=noise_model.scalar noise_threshold=noise_model.scalar * resolved_cfg.noise_sigma_multiplier noise_sigma_multiplier=resolved_cfg.noise_sigma_multiplier tile_size=resolved_cfg.noise_tile_size

        # Chunk processing writes intermediate Q/U/W/I and diagnostics to disk so
        # large cubes do not need a full in-memory working cube.
        for j in 1:resolved_cfg.ny_chunks
            process_chunk(j, resolved_cfg, ny, nv; noise_model=noise_model, channels=channels)
        end
        q_stitched, u_stitched, intensity_stitched, weight_stitched, diagnostics, fit_maps =
            stitch_chunk_outputs(resolved_cfg, nx, ny)
        finalized = finalize_vgt_products(q_stitched, u_stitched, intensity_stitched, resolved_cfg.ker,
                                          resolved_cfg.zero_stokes_policy,
                                          resolved_cfg.output_support_mode;
                                          edge_mask=resolved_cfg.edge_mask,
                                          fitted_weight=weight_stitched)
        products = if resolved_cfg.output_bin == 1
            (; finalized.psi, finalized.reliability, fit_maps.fit_rbar2, fit_maps.fit_support, header)
        else
            bin_products(finalized, fit_maps.fit_rbar2, fit_maps.fit_support, header, resolved_cfg.output_bin)
        end
        # Gradients and Q/U smoothing are done in the pixel grid; a sky frame is
        # applied only to the final angles, per pixel of the output grid.
        psi = resolved_cfg.angle_frame === :pixel ? products.psi :
              pixel_to_position_angle(products.psi, products.header; frame=resolved_cfg.angle_frame)
        finished_at = Dates.now(Dates.UTC)
        performance = (;
            julia_threads=Threads.nthreads(),
            ny_chunks=resolved_cfg.ny_chunks,
            elapsed_seconds=round((time_ns() - start_ns) / 1e9; digits=3),
            started_at_utc=string(Dates.format(started_at, dateformat"yyyy-mm-ddTHH:MM:SS"), "Z"),
            finished_at_utc=string(Dates.format(finished_at, dateformat"yyyy-mm-ddTHH:MM:SS"), "Z"),
        )
        return PipelineRun(
            config=resolved_cfg,
            psi=psi,
            reliability=Float32.(products.reliability),
            Q=q_stitched,
            U=u_stitched,
            I=intensity_stitched,
            W=weight_stitched,
            fit_rbar2=products.fit_rbar2,
            fit_support=products.fit_support,
            header=products.header,
            diagnostics=diagnostics,
            channels=channels,
            performance=performance,
        )
    catch err
        if is_interrupt_exception(err)
            cleanup = try
                cleanup_interrupted_run!(resolved_cfg.output_dir; remove_created_output_dir=output_dir_created)
            catch cleanup_err
                @warn "Interrupt cleanup failed. Partial chunk outputs may remain in $(resolved_cfg.output_dir)." exception=(cleanup_err, catch_backtrace())
                nothing
            end
            if cleanup === nothing
                @warn "Interrupt received. Stopping after failed cleanup."
            else
                @warn "Interrupt received. Removed partial chunk outputs." output_dir=resolved_cfg.output_dir removed_chunks=length(cleanup.removed_chunks) removed_output_dir=cleanup.removed_output_dir
            end
            rethrow(InterruptException())
        end
        rethrow()
    end
end

function cleanup_interrupted_run!(output_dir::AbstractString; remove_created_output_dir::Bool=false)
    removed_chunks = String[]
    isdir(output_dir) || return (; removed_chunks, removed_output_dir=false)

    for name in readdir(output_dir)
        if occursin(r"^chunk_2d_[0-9]+\.jld2(\.tmp\..+)?$", name)
            path = joinpath(output_dir, name)
            rm(path; force=true)
            push!(removed_chunks, path)
        end
    end

    removed_output_dir = false
    if remove_created_output_dir && isdir(output_dir) && isempty(readdir(output_dir))
        rm(output_dir)
        removed_output_dir = true
    end

    return (; removed_chunks, removed_output_dir)
end

function is_interrupt_exception(err)
    err isa InterruptException && return true
    err isa CompositeException && return any(is_interrupt_exception, err.exceptions)
    err isa TaskFailedException && return is_interrupt_exception(err.task.exception)
    return false
end
