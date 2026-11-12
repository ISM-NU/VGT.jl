using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "VGT")

using VGT: load_config, read_cube_metadata, getcoor, is_interrupt_exception, run_pipeline, write_pipeline_outputs

function main()
    # main.jl is intentionally thin: it loads configuration, prints run
    # metadata, calls the package pipeline, and writes products.
    # Optional first argument: config path (HPC jobs pass one per cube).
    cfg = load_config(isempty(ARGS) ? "config.toml" : ARGS[1])

    try
        nx, ny, nv, header = read_cube_metadata(cfg.input_file)
        println("File: $(abspath(cfg.input_file))")
        println("Size: $(round(filesize(cfg.input_file) / 1024^3, digits=2)) GB")
        println("Data size: ($nx, $ny, $nv)")

        try
            getcoor(header, nx, ny, nv)
        catch err
            is_interrupt_exception(err) && rethrow()
            @warn "Coordinate extraction failed. Skipping metadata."
        end

        println("Starting Low-Memory Parallel VGT")
        run = run_pipeline(cfg)

        println("Writing outputs...")
        outputs = write_pipeline_outputs(run)
        println("Success. Results saved to $(outputs.psi_path)")
        println("Reliability map saved to $(outputs.reliability_path)  (P=sqrt(Q²+U²)/I, range 0-1)")
        println("Fit diagnostics: attempts=$(run.diagnostics.fit_attempts) successes=$(run.diagnostics.fit_successes) failures=$(run.diagnostics.fit_failures)")
        return 0
    catch err
        if err isa InterruptException
            return 130
        end
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
