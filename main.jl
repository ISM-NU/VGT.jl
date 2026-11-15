using DrWatson
Base.exit_on_sigint(false)
quickactivate(@__DIR__, "VGT")

using VGT: load_config, read_cube_metadata, getcoor, is_interrupt_exception, run_pipeline, write_pipeline_outputs

function run_one(cfg)
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
end

function main(args=ARGS)
    # main.jl is intentionally thin: it loads configuration, prints run
    # metadata, calls the package pipeline, and writes products.
    # Arguments: config paths (default config.toml). Several configs run one
    # after another in this Julia process, so Julia starts and compiles once;
    # for short runs (velocity windows, layer cubes) that start-up is most of
    # the wall time. All configs are loaded first so a bad one fails before any
    # run; a run that fails is reported and the next one still runs.
    paths = isempty(args) ? ["config.toml"] : collect(args)
    cfgs = [load_config(p) for p in paths]
    failed = String[]
    for (path, cfg) in zip(paths, cfgs)
        length(paths) > 1 && println("== $(path)")
        try
            run_one(cfg)
        catch err
            (err isa InterruptException || is_interrupt_exception(err)) && return 130
            length(paths) == 1 && rethrow()
            @error "Run failed" config=path exception=(err, catch_backtrace())
            push!(failed, path)
        end
    end
    isempty(failed) || println("Failed configs: $(join(failed, ", "))")
    return isempty(failed) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
