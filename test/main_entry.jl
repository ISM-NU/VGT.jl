using FITSIO
using Test

# main.jl with several configs: all run in one process, a bad config fails before any run,
# and a run that fails does not stop the next one.
@testset "main.jl runs several configs" begin
    entry = Module(:MainEntry)
    Base.include(entry, joinpath(@__DIR__, "..", "main.jl"))
    mktempdir() do dir
        cube = zeros(Float32, 16, 14, 4)
        @inbounds for k in 1:4, j in 1:14, i in 1:16
            cube[i, j, k] = 3f0 + 0.2f0 * i + 0.1f0 * j * k + 0.01f0 * i * j
        end
        path = joinpath(dir, "cube.fits")
        FITS(f -> write(f, cube), path, "w")
        conf(name, input) = (p = joinpath(dir, "$name.toml");
                             write(p, """
                                 [input]
                                 input_file = "$input"
                                 output_dir = "$(joinpath(dir, name))"
                                 [noise]
                                 noise = 1e-6
                                 [vgt]
                                 dn = 4
                                 ker = 1
                                 min_valid = 4
                                 [execution]
                                 ny_chunks = 1
                                 """); p)
        a, b = conf("a", path), conf("b", path)
        @test redirect_stdout(devnull) do; entry.main([a, b]); end == 0
        @test all(isfile(joinpath(dir, n, "run_metadata.toml")) for n in ("a", "b"))

        missing_input = conf("missing", joinpath(dir, "nope.fits"))
        c = conf("c", path)
        @test redirect_stdout(devnull) do; entry.main([missing_input, c]); end == 1   # failed run, next still runs
        @test isfile(joinpath(dir, "c", "run_metadata.toml"))

        bad = joinpath(dir, "bad.toml")
        write(bad, "[noise]\nnoise = 1.0\n")                                       # no input_file
        d = conf("d", path)
        @test_throws Exception entry.main([d, bad])
        @test !isdir(joinpath(dir, "d"))                                         # nothing ran
    end
end
