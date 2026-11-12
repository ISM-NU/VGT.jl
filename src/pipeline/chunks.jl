# Stitch per-chunk intermediate products into full-size maps.
#
# Chunking is an execution strategy for memory control. The scientific output
# should match a full-cube run once padded chunks are trimmed and stitched.

function stitch_chunk_outputs(cfg::VGTConfig, nx::Int, ny::Int)
    q_stitched = zeros(Float64, nx, ny)
    u_stitched = zeros(Float64, nx, ny)
    intensity_stitched = zeros(Float64, nx, ny)
    fit_rbar2_stitched = fill(NaN, nx, ny)
    fit_attempt_fraction_stitched = zeros(Float64, nx, ny)
    fit_support_stitched = zeros(Float64, nx, ny)
    fit_attempts = 0
    fit_successes = 0
    fit_failures = 0
    valid_windows = 0
    fit_low_confidence_windows = 0

    for j in 1:cfg.ny_chunks
        filename = joinpath(cfg.output_dir, "chunk_2d_$(j).jld2")
        if !isfile(filename)
            continue
        end

        y_bounds = chunk_y_bounds(j, ny, cfg.ny_chunks)
        y_bounds === nothing && continue
        y_start, _ = y_bounds

        jldopen(filename, "r") do file
            # Each chunk stores only its unpadded sky rows, so it can be copied
            # directly into the matching full-map y range.
            q_chunk_accum = file["q"]
            y_end = y_start + size(q_chunk_accum, 2) - 1
            q_stitched[:, y_start:y_end] = q_chunk_accum
            u_stitched[:, y_start:y_end] = file["u"]
            intensity_stitched[:, y_start:y_end] = file["i"]
            if haskey(file, "fit_rbar2") && file["fit_rbar2"] !== nothing
                fit_rbar2_stitched[:, y_start:y_end] = file["fit_rbar2"]
            end
            if haskey(file, "fit_attempt_fraction") && file["fit_attempt_fraction"] !== nothing
                fit_attempt_fraction_stitched[:, y_start:y_end] = file["fit_attempt_fraction"]
            end
            if haskey(file, "fit_support") && file["fit_support"] !== nothing
                fit_support_stitched[:, y_start:y_end] = file["fit_support"]
            end
            fit_attempts += haskey(file, "fit_attempts") ? file["fit_attempts"] : 0
            fit_successes += haskey(file, "fit_successes") ? file["fit_successes"] : 0
            fit_failures += haskey(file, "fit_failures") ? file["fit_failures"] : 0
            valid_windows += haskey(file, "valid_windows") ? file["valid_windows"] : 0
            fit_low_confidence_windows += haskey(file, "fit_low_confidence_windows") ? file["fit_low_confidence_windows"] : 0
        end
    end

    diagnostics = VGTDiagnostics(
        # Summed diagnostics describe the whole run, not individual chunks.
        fit_attempts=fit_attempts,
        fit_successes=fit_successes,
        fit_failures=fit_failures,
        valid_windows=valid_windows,
        fit_low_confidence_windows=fit_low_confidence_windows,
    )
    fit_maps = (;
        fit_rbar2=fit_rbar2_stitched,
        fit_attempt_fraction=fit_attempt_fraction_stitched,
        fit_support=fit_support_stitched,
    )
    return q_stitched, u_stitched, intensity_stitched, diagnostics, fit_maps
end
