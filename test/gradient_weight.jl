using VGT
using Random
using Test

# gradient_weight = :amplitude2 weights each gradient's doubled angle by |g|^2 (the
# structure tensor); docs/validation-faint-gas.md. The default :unit path must not change.

@testset "gradient_weight fitter" begin
    rng = MersenneTwister(3)
    gx = randn(rng, 20, 20) .+ 0.7
    gy = randn(rng, 20, 20) .- 0.2
    gx[3, 4] = NaN
    unit_default = VGT.fit_window_result(:axial_vonmises, gx, gy)
    unit_explicit = VGT.fit_window_result(:axial_vonmises, gx, gy; gradient_weight=:unit)
    @test unit_default == unit_explicit

    amp = VGT.fit_window_result(:axial_vonmises, gx, gy; gradient_weight=:amplitude2)
    ok = isfinite.(gx) .& isfinite.(gy)
    c = sum(gx[ok] .^ 2 .- gy[ok] .^ 2)
    s = sum(2 .* gx[ok] .* gy[ok])
    @test amp.fit_ok
    @test amp.angle ≈ 0.5 * atan(s, c) atol = 1e-12          # principal axis of the structure tensor
    @test amp.rbar2 ≈ hypot(s, c) / sum(gx[ok] .^ 2 .+ gy[ok] .^ 2)
    @test amp.effective_n == count(ok)

    # One strong gradient at 0 rad against many weak ones at pi/2: |g|^2 weighting follows
    # the strong one, unit weighting the many.
    wx = vcat(10.0, zeros(9))
    wy = vcat(0.0, fill(1.0, 9))
    @test abs(VGT.fit_window_result(:axial_vonmises, wx, wy; gradient_weight=:amplitude2).angle) < 1e-12
    @test abs(abs(VGT.fit_window_result(:axial_vonmises, wx, wy).angle) - pi / 2) < 1e-12

    # Equal magnitudes: the two weightings agree.
    phi = 0.3 .+ 0.05 .* randn(rng, 50)
    a1 = VGT.fit_window_result(:axial_vonmises, cos.(phi), sin.(phi))
    a2 = VGT.fit_window_result(:axial_vonmises, cos.(phi), sin.(phi); gradient_weight=:amplitude2)
    @test a1.angle ≈ a2.angle atol = 1e-12
    @test a1.rbar2 ≈ a2.rbar2 atol = 1e-12

    @test_throws ArgumentError VGT.fit_window_result(:axial_vonmises, gx, gy; gradient_weight=:amplitude)
    @test !VGT.fit_window_result(:axial_vonmises, [0.0, 0.0], [0.0, 0.0]; gradient_weight=:amplitude2).fit_ok
end

@testset "gradient_weight kernel" begin
    # Noise-free ramp: one gradient direction everywhere, so both weightings give the same map.
    nx, ny, nv = 40, 40, 3
    cube = [10.0 + 0.3 * i + 0.1 * j + k for i in 1:nx, j in 1:ny, k in 1:nv]
    q1, u1 = VGT.VGT_QU_mv_2d(cube, 10, 1, nv, 1.0; min_valid=50)
    q2, u2 = VGT.VGT_QU_mv_2d(cube, 10, 1, nv, 1.0; min_valid=50, gradient_weight=:amplitude2)
    @test q1 ≈ q2 atol = 1e-9
    @test u1 ≈ u2 atol = 1e-9
    @test_throws ArgumentError VGT.VGT_QU_mv_2d(cube, 10, 1, nv, 1.0; gradient_weight=:bogus)
end

@testset "gradient_weight config" begin
    mktempdir() do dir
        path = joinpath(dir, "config.toml")
        base = """
        [input]
        input_file = "data/exp_raw/example.fits"
        output_dir = ""

        [noise]
        noise = 0.17

        [vgt]
        dn = 20
        ker = 4
        min_valid = 100

        [execution]
        ny_chunks = 1
        """
        write(path, base)
        cfg = VGT.load_config(path)
        @test cfg.gradient_weight === :unit
        @test !haskey(VGT.run_parameters(cfg), "gradient_weight")   # default run names unchanged
        write(path, replace(base, "min_valid = 100" => "min_valid = 100\ngradient_weight = \"amplitude2\""))
        acfg = VGT.load_config(path)
        @test acfg.gradient_weight === :amplitude2
        @test VGT.run_parameters(acfg)["gradient_weight"] == "amplitude2"
        @test VGT.config_with(acfg; dn=10).gradient_weight === :amplitude2
        write(path, replace(base, "min_valid = 100" => "min_valid = 100\ngradient_weight = \"amplitude\""))
        @test_throws ArgumentError VGT.load_config(path)
    end
end
