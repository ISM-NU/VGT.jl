using Test

@testset "astroVGT" begin
    include("configuration.jl")
    include("fitters.jl")
    include("fitter_results.jl")
    include("noise_thresholds.jl")
    include("script_smoke.jl")
    include("interrupt_handling.jl")
    include("vgt_edge_cases.jl")
    include("write_outputs.jl")
    include("sky_angles.jl")
    include("output_bin.jl")
    include("velocity_range.jl")
    include("gradient_whitening.jl")
    include("main_entry.jl")
    include("regrid.jl")
    include("pipeline_fit_diagnostics.jl")
    include("chunk_seams.jl")
end
