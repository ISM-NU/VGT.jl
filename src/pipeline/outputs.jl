# Output naming and file writing.
#
# This file turns a PipelineRun into FITS/PNG products. It should not change the
# numerical content of the science products.

function write_pipeline_outputs(run::PipelineRun)
    # Write the routine science bundle: primary angle map, reliability map,
    # the core fit-coherence diagnostic, and a quick-look PNG.
    output_dir = resolved_output_dir(run.config)
    mkpath(output_dir)
    psi_path = next_available_output_path(
        output_dir,
        prefixed_input_filename(run.config.input_file, "psi_", ".fits"),
    )
    reliability_path = next_available_output_path(
        output_dir,
        prefixed_input_filename(run.config.input_file, "psi_reliability_", ".fits"),
    )
    plot_path = next_available_output_path(
        output_dir,
        prefixed_input_filename(run.config.input_file, "psi_plot_", ".png"),
    )
    fit_rbar2_path = next_available_output_path(
        output_dir,
        prefixed_input_filename(run.config.input_file, "psi_fit_rbar2_", ".fits"),
    )
    support_path = next_available_output_path(
        output_dir,
        prefixed_input_filename(run.config.input_file, "psi_support_", ".fits"),
    )
    metadata_path = next_available_output_path(output_dir, "run_metadata.toml")

    product_header(product) = vgt_product_header(run.header, run.config, product)
    write_fits(psi_path, run.psi; header=product_header(:psi))
    write_fits(reliability_path, run.reliability; header=product_header(:reliability))
    write_fits(fit_rbar2_path, run.fit_rbar2; header=product_header(:fit_rbar2))
    write_fits(support_path, Float32.(run.fit_support); header=product_header(:support))

    # An all-NaN angle map (e.g. a velocity window with no emission above the
    # noise mask) is a valid empty result, but Makie cannot plot it: the color
    # range is NaN. Skip the quick-look and record that in the metadata.
    if any(isfinite, run.psi)
        fig, ax, hm = heatmap(run.psi, colormap=:jet)
        Colorbar(fig[1, 2], hm, label="Psi (rad)")
        save(plot_path, fig)
    else
        @warn "no finite angles; plot skipped" plot_path
        plot_path = nothing
    end

    outputs = (;
        psi_path,
        reliability_path,
        plot_path,
        fit_rbar2_path,
        support_path,
        metadata_path,
    )
    write_run_metadata(metadata_path, run, outputs)

    if !run.config.keep_intermediate_chunks
        # Remove temporary per-chunk files only after final products were
        # written successfully.
        for j in 1:run.config.ny_chunks
            rm(joinpath(output_dir, "chunk_2d_$(j).jld2"); force=true)
        end
    end

    return outputs
end

function run_metadata(run::PipelineRun, outputs)
    cfg = run.config
    metadata = Dict{String, Any}(
        "generated_at_utc" => string(Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS"), "Z"),
        "project" => Dict{String, Any}(
            "name" => "VGT",
            "drwatson_projectdir" => projectdir(),
        ),
        "input" => Dict{String, Any}(
            "path" => cfg.input_file,
            "absolute_path" => abspath(cfg.input_file),
            "stem" => splitext(basename(cfg.input_file))[1],
            "size_bytes" => isfile(cfg.input_file) ? filesize(cfg.input_file) : -1,
            "weight_file" => cfg.weight_file,
            "angle_frame" => String(cfg.angle_frame),
            "channels_used" => isempty(run.channels) ? Int[] : [first(run.channels), last(run.channels)],
        ),
        "run_parameters" => run_parameters(cfg),
        "config" => _config_metadata(cfg),
        "performance" => _performance_metadata(run),
        "diagnostics" => Dict{String, Any}(
            "fit_attempts" => run.diagnostics.fit_attempts,
            "fit_successes" => run.diagnostics.fit_successes,
            "fit_failures" => run.diagnostics.fit_failures,
            "valid_windows" => run.diagnostics.valid_windows,
            "fit_low_confidence_windows" => run.diagnostics.fit_low_confidence_windows,
        ),
        "outputs" => Dict{String, Any}(
            "psi_path" => outputs.psi_path,
            "reliability_path" => outputs.reliability_path,
            "fit_rbar2_path" => outputs.fit_rbar2_path,
            "support_path" => outputs.support_path,
            "metadata_path" => outputs.metadata_path,
        ),
    )
    # TOML has no null: the key is absent when the plot was skipped.
    outputs.plot_path === nothing || (metadata["outputs"]["plot_path"] = outputs.plot_path)
    tag!(metadata; gitpath=projectdir(), warn=false, storepatch=false)
    return metadata
end

const _VGT_PRODUCT_CARDS = Dict(
    :psi => ("rad", "VGT magnetic-field orientation proxy (axial)"),
    :reliability => ("", "smoothed Stokes coherence P = sqrt(Q^2+U^2)/I"),
    :fit_rbar2 => ("", "fitter doubled-angle coherence rbar^2"),
    :support => ("", "fraction of channels with a successful VGT fit"),
)

"""
    vgt_product_header(input_header, cfg, product) -> FITSHeader

Header for a 2D VGT product: the input cube's celestial WCS (see
`celestial_header`) plus the product's own unit, angle convention, and the VGT
parameters that set its effective resolution. The grid is the input grid, or
`output_bin` x `output_bin` blocks of it (`input_header` is then the binned
header from `bin_products`); `dn` and `ker` are always in input pixels.
"""
function vgt_product_header(input_header, cfg::VGTConfig, product::Symbol)
    haskey(_VGT_PRODUCT_CARDS, product) || throw(ArgumentError("Unknown VGT product $(product)"))
    header = input_header === nothing ? FITSHeader(String[], [], String[]) : celestial_header(input_header)
    bunit, description = _VGT_PRODUCT_CARDS[product]
    cards = Tuple{String, Any, String}[
        ("BUNIT", bunit, product === :psi ? "radians" : "dimensionless"),
        ("VGTPROD", String(product), description),
        ("VGTDN", cfg.dn, "VGT gradient window side, input pixels"),
        ("VGTKER", Float64(cfg.ker), "Q/U Gaussian smoothing sigma, input pixels"),
        ("VGTKFWHM", 2 * sqrt(2 * log(2)) * cfg.ker, "Q/U Gaussian smoothing FWHM, input pixels"),
        ("VGTMINV", cfg.min_valid, "min valid gradients per window"),
        ("VGTFIT", String(cfg.fitter_mode), "window angle estimator"),
        ("VGTGRAD", String(cfg.gradient_mode), "gradient mode"),
        ("VGTNMODE", String(cfg.noise_mode), "noise model"),
        ("VGTNOISE", Float64(cfg.noise), "noise threshold input (K)"),
        ("VGTBIN", cfg.output_bin, "output pixel = VGTBIN x VGTBIN input pixels"),
    ]
    if !isempty(cfg.velocity_range_kms)
        push!(cards, ("VGTVLO", cfg.velocity_range_kms[1], "velocity window low end (km/s)"))
        push!(cards, ("VGTVHI", cfg.velocity_range_kms[2], "velocity window high end (km/s)"))
    end
    if product === :psi
        if cfg.angle_frame === :pixel
            push!(cards, ("PSIFRAME", "PIXEL", "CCW from +axis 1 in pixel grid; not east of N"))
        else
            push!(cards, ("PSIFRAME", uppercase(String(cfg.angle_frame)), "IAU position angle, east of north in this frame"))
        end
        push!(cards, ("PSIANGLE", "BFIELD", "gradient angle + 90 deg (B-field, not pol.)"))
    end
    for (key, value, comment) in cards
        header[key] = value
        set_comment!(header, key, comment)
    end
    return header
end

function write_run_metadata(path::AbstractString, run::PipelineRun, outputs)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, run_metadata(run, outputs); sorted=true)
    end
    return path
end

function _config_metadata(cfg::VGTConfig)
    return Dict{String, Any}(
        "input_file" => cfg.input_file,
        "output_dir" => cfg.output_dir,
        "keep_intermediate_chunks" => cfg.keep_intermediate_chunks,
        "noise" => cfg.noise,
        "dn" => cfg.dn,
        "ker" => cfg.ker,
        "ny_chunks" => _metadata_value(cfg.ny_chunks),
        "min_valid" => cfg.min_valid,
        "fitter_mode" => _metadata_value(cfg.fitter_mode),
        "rbar2_weight" => cfg.rbar2_weight,
        "rbar2_gate" => cfg.rbar2_gate,
        "noise_mode" => _metadata_value(cfg.noise_mode),
        "noise_tile_size" => cfg.noise_tile_size,
        "noise_sigma_multiplier" => cfg.noise_sigma_multiplier,
        "zero_stokes_policy" => _metadata_value(cfg.zero_stokes_policy),
        "output_support_mode" => _metadata_value(cfg.output_support_mode),
        "gradient_mode" => _metadata_value(cfg.gradient_mode),
        "angle_frame" => _metadata_value(cfg.angle_frame),
        "output_bin" => cfg.output_bin,
        "velocity_range_kms" => cfg.velocity_range_kms,
        "gradient_whiten" => cfg.gradient_whiten,
        "gradient_noise_cov" => cfg.gradient_noise_cov,
        "gradient_noise_cov_file" => cfg.gradient_noise_cov_file,
    )
end

function _performance_metadata(run::PipelineRun)
    return Dict{String, Any}(
        "julia_threads" => run.performance.julia_threads,
        "ny_chunks" => run.performance.ny_chunks,
        "elapsed_seconds" => run.performance.elapsed_seconds,
        "started_at_utc" => run.performance.started_at_utc,
        "finished_at_utc" => run.performance.finished_at_utc,
    )
end

_metadata_value(value::Symbol) = String(value)
_metadata_value(value) = value

function prefixed_input_filename(input_path::AbstractString, prefix::AbstractString, ext::AbstractString)
    stem = splitext(basename(input_path))[1]
    return string(prefix, stem, ext)
end

function next_available_output_path(dir::AbstractString, filename::AbstractString)
    # Avoid overwriting an existing product from a previous run.
    candidate = joinpath(dir, filename)
    isfile(candidate) || return candidate

    stem, ext = splitext(filename)
    suffix = 1
    while true
        candidate = joinpath(dir, "$(stem)_$(suffix)$(ext)")
        isfile(candidate) || return candidate
        suffix += 1
    end
end
