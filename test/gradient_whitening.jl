using VGT
using Test
using FITSIO
using Random
using Statistics
using LinearAlgebra: det

# Deterministic smooth cube (no RNG) for the bitwise-default check.
function gradient_whitening_reference_cube(nx=40, ny=40, nv=3)
    cube = zeros(Float64, nx, ny, nv)
    for k in 1:nv, j in 1:ny, i in 1:nx
        cube[i, j, k] = 20 + 0.3 * (cos(0.5) * i + sin(0.5) * j) +
                        0.4 * sin(1.7i + 2.3j + 3k) + 0.3 * cos(0.9i * j / 7 + k)
    end
    return cube
end

# White noise convolved with an elliptical Gaussian (long axis at `long_axis`).
function anisotropic_noise(nx, ny, nv, rng; σ_long=1.6, σ_short=1.0, long_axis=0.8)
    r = ceil(Int, 4 * σ_long)
    c, s = cos(long_axis), sin(long_axis)
    kern = [exp(-0.5 * (((x * c + y * s) / σ_long)^2 + ((-x * s + y * c) / σ_short)^2))
            for x in -r:r, y in -r:r]
    kern ./= sqrt(sum(abs2, kern))
    out = Array{Float64}(undef, nx, ny, nv)
    for k in 1:nv
        white = randn(rng, nx + 2r, ny + 2r)
        for j in 1:ny, i in 1:nx
            acc = 0.0
            for jj in -r:r, ii in -r:r
                acc += kern[ii + r + 1, jj + r + 1] * white[i + r + ii, j + r + jj]
            end
            out[i, j, k] = acc
        end
    end
    return out
end

@testset "gradient_whiten is off by default and leaves output unchanged" begin
    cube = gradient_whitening_reference_cube()
    # Values from the kernel before the option existed (origin/main).
    # rtol allows for FMA differences between CPUs.
    reference = Dict(
        :central_checked => (21.817856594602002, 81.01800180202324),
        :sobel_checked => (52.51534797677656, 65.750181254463),
        :nan_edge_backward => (-4.793401294596194, 83.61593584590709),
    )
    for (mode, (q_ref, u_ref)) in reference
        q, u = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_mode=mode)
        @test q[20, 20] ≈ q_ref rtol = 1e-10
        @test u[20, 20] ≈ u_ref rtol = 1e-10
        q0, u0 = VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_mode=mode,
                              gradient_whiten=false, gradient_noise_cov=[1.0, 2.0, 0.3])
        @test q0 == q && u0 == u
    end
end

@testset "whitening matrices" begin
    for cov in ([1.0, 1.0, 0.0], [0.3, 0.2, 0.05], [5.0, 0.4, -1.2])
        m = VGT._gradient_whitening(cov)
        W = [m.w11 m.w12; m.w12 m.w22]
        S = [m.s11 m.s12; m.s12 m.s22]
        C = [cov[1] cov[3]; cov[3] cov[2]]
        @test W * S ≈ [1.0 0.0; 0.0 1.0] atol = 1e-12
        @test det(S) ≈ 1 atol = 1e-12
        # W maps the covariance to a multiple of the identity (unit-determinant scale).
        @test W * C * W ≈ sqrt(det(C)) * [1.0 0.0; 0.0 1.0] atol = 1e-12
        for θ in (-1.2, 0.0, 0.4, 1.5)
            gx, gy = VGT._whiten_gradient(m, cos(θ), sin(θ))
            back = VGT._unwhiten_angle(m, atan(gy, gx))
            @test VGT.headless_angle_difference(back, θ) ≈ 0 atol = 1e-12
        end
    end
end

@testset "whitening with isotropic covariance equals off" begin
    rng = MersenneTwister(3)
    cube = 50.0 .+ [0.3 * i + 0.2 * j for i in 1:48, j in 1:48, k in 1:2] .+ randn(rng, 48, 48, 2)
    for mode in (:central_checked, :sobel_checked)
        q0, u0 = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_mode=mode)
        for cov in ([1.0, 1.0, 0.0], [0.37, 0.37, 0.0])  # scale does not matter
            q1, u1 = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_mode=mode,
                                  gradient_whiten=true, gradient_noise_cov=cov)
            @test q1 ≈ q0 rtol = 1e-12 atol = 1e-9
            @test u1 ≈ u0 rtol = 1e-12 atol = 1e-9
        end
    end
end

@testset "whitening removes the preferred direction of anisotropic noise" begin
    rng = MersenneTwister(11)
    noise = anisotropic_noise(96, 96, 12, rng)
    m = measure_gradient_noise(noise)
    @test m.anisotropy > 0.1
    A_check, dir_check = gradient_anisotropy(m.cov)
    @test A_check == m.anisotropy && dir_check == m.direction_deg
    w = VGT._gradient_whitening(m.cov)
    nx, ny, nv = size(noise)
    z_off = 0.0im; z_on = 0.0im; n = 0
    zw_off = 0.0im; zw_on = 0.0im; p_off = 0.0; p_on = 0.0
    for k in 1:nv, j in 2:(ny - 1), i in 2:(nx - 1)
        gx = (noise[i - 1, j, k] - noise[i + 1, j, k]) / 2
        gy = (noise[i, j - 1, k] - noise[i, j + 1, k]) / 2
        wx, wy = VGT._whiten_gradient(w, gx, gy)
        n += 1
        z_off += cis(2 * atan(gy, gx))
        z_on += cis(2 * atan(wy, wx))
        zw_off += complex(gx, gy)^2; p_off += gx^2 + gy^2
        zw_on += complex(wx, wy)^2; p_on += wx^2 + wy^2
    end
    # Unit-weight mean resultant of the doubled angle (what VGT averages). The
    # noise gradients are correlated over a few pixels, so allow a generous null scale.
    null_sd = 1 / sqrt(n) * 4
    @test abs(z_off) / n > 10 * null_sd
    @test abs(z_on) / n < 3 * null_sd
    @test abs(zw_on) / p_on < 0.05 * (abs(zw_off) / p_off)
end

@testset "whitening reduces the bias toward anisotropic noise in the kernel" begin
    rng = MersenneTwister(5)
    nx = 96
    noise = anisotropic_noise(nx, nx, 6, rng; σ_long=1.8, σ_short=1.2, long_axis=0.9)
    m = measure_gradient_noise(noise)
    θ_true = deg2rad(20.0)
    G = 0.5 * m.sigma_x
    cube = [1000.0 + G * (cos(θ_true) * (i - 48) + sin(θ_true) * (j - 48)) for i in 1:nx, j in 1:nx, k in 1:6] .+ noise
    function bias(; whiten)
        q, u = VGT_QU_mv_2d(cube, 20, 1, 6, 1.0; min_valid=100,
                            gradient_whiten=whiten, gradient_noise_cov=whiten ? m.cov : Float64[])
        keep = (q .!= 0) .| (u .!= 0)
        z = mean(cis.(2 .* (0.5 .* atan.(u[keep], q[keep]) .- θ_true)))
        return abs(rad2deg(0.5 * angle(z)))
    end
    b_off = bias(whiten=false)
    b_on = bias(whiten=true)
    @test b_off > 8
    @test b_on < b_off / 2
end

@testset "gradient_whiten config validation" begin
    base = (input_file="c.fits", noise=1.0, dn=4, ker=1, min_valid=4)
    cfg = VGTConfig(; base...)
    @test !cfg.gradient_whiten && isempty(cfg.gradient_noise_cov)
    on = VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov=[0.3, 0.2, 0.05])
    @test on.gradient_whiten && on.gradient_noise_cov == [0.3, 0.2, 0.05]
    @test VGT.config_with(on; dn=6).gradient_noise_cov == [0.3, 0.2, 0.05]
    err = try VGTConfig(; base..., gradient_whiten=true); nothing catch e e end
    @test err isa ArgumentError && occursin("gradient_noise_cov", err.msg)
    for bad in ([1.0, 1.0], [1.0, 1.0, 1.0], [1.0, -1.0, 0.0], [NaN, 1.0, 0.0], [1.0, 1.0, 2.0, 0.0])
        @test_throws ArgumentError VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov=bad)
    end
    @test_throws ArgumentError VGTConfig(; base..., gradient_whiten=true,
                                         gradient_noise_cov=[1.0, 1.0, 0.0],
                                         gradient_mode=:nan_edge_backward)
    cube = gradient_whitening_reference_cube()
    @test_throws ArgumentError VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_whiten=true)
    @test_throws ArgumentError VGT_QU_mv_2d(cube, 6, 1, 3, 1.0; min_valid=10, gradient_whiten=true,
                                            gradient_noise_cov=[1.0, 1.0, 0.0],
                                            gradient_mode=:nan_edge_backward)

    # Default run names and hashes are unchanged; the option changes them.
    @test !haskey(VGT.run_parameters(cfg), "gradient_whiten")
    @test VGT.run_parameter_hash(cfg) != VGT.run_parameter_hash(on)
    meta = VGT._config_metadata(on)
    @test meta["gradient_whiten"] === true && meta["gradient_noise_cov"] == [0.3, 0.2, 0.05]

    mktempdir() do dir
        path = joinpath(dir, "config.toml")
        write(path, """
        [input]
        input_file = "x.fits"
        output_dir = ""
        [noise]
        noise = 0.17
        gradient_whiten = true
        gradient_noise_cov = [0.3, 0.2, 0.05]
        [vgt]
        dn = 20
        ker = 4
        min_valid = 100
        [execution]
        ny_chunks = 1
        """)
        loaded = load_config(path)
        @test loaded.gradient_whiten && loaded.gradient_noise_cov == [0.3, 0.2, 0.05]
        write(path, replace(read(path, String), "gradient_noise_cov = [0.3, 0.2, 0.05]\n" => ""))
        @test_throws ArgumentError load_config(path)
    end
end

@testset "spatially varying covariance map" begin
    rng = MersenneTwister(21)
    nx = 64
    # Constant map equals the constant covariance option.
    cube = 50.0 .+ [0.3 * i + 0.2 * j for i in 1:nx, j in 1:nx, k in 1:2] .+ randn(rng, nx, nx, 2)
    cov = [0.30, 0.20, 0.05]
    constmap = cat(fill(cov[1], nx, nx), fill(cov[2], nx, nx), fill(cov[3], nx, nx); dims=3)
    for mode in (:central_checked, :sobel_checked)
        q1, u1 = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_mode=mode,
                              gradient_whiten=true, gradient_noise_cov=cov)
        q2, u2 = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_mode=mode,
                              gradient_whiten=true, gradient_noise_cov_map=constmap)
        @test q2 ≈ q1 rtol = 1e-12 atol = 1e-9
        @test u2 ≈ u1 rtol = 1e-12 atol = 1e-9
    end
    # Bad maps fail clearly.
    bad = copy(constmap); bad[5, 5, 2] = -1.0
    @test_throws ArgumentError VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20,
                                            gradient_whiten=true, gradient_noise_cov_map=bad)
    @test_throws ArgumentError VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20,
                                            gradient_whiten=true, gradient_noise_cov_map=constmap[1:10, :, :])
    @test_throws ArgumentError VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_whiten=true,
                                            gradient_noise_cov=cov, gradient_noise_cov_map=constmap)

    # Pixels without a covariance (NaN, e.g. outside the survey coverage) drop out
    # instead of failing the run; windows that do not reach them are unchanged.
    nanmap = copy(constmap); nanmap[1:15, :, :] .= NaN
    qn, un = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_whiten=true, gradient_noise_cov_map=nanmap)
    qc, uc = VGT_QU_mv_2d(cube, 12, 1, 2, 1.0; min_valid=20, gradient_whiten=true, gradient_noise_cov_map=constmap)
    @test all(isfinite, qn) && all(isfinite, un)
    @test all(iszero, qn[1:15, :]) && all(iszero, un[1:15, :])
    @test qn[23:end, :] ≈ qc[23:end, :] rtol = 1e-12 atol = 1e-9
    @test any(!iszero, qn[23:end, :])

    # The helper recovers a constant covariance and ignores NaN pixels.
    noise = anisotropic_noise(80, 80, 16, MersenneTwister(2))
    m = measure_gradient_noise(noise)
    mp = gradient_noise_cov_map(noise; fwhm_pixels=30)
    @test size(mp) == (80, 80, 3)
    @test all(isfinite, mp)
    @test mean(mp[:, :, 1]) ≈ m.cov[1] rtol = 0.05
    @test mean(mp[:, :, 3]) ≈ m.cov[3] rtol = 0.2
    holes = copy(noise); holes[30:40, 30:40, :] .= NaN
    mh = gradient_noise_cov_map(holes; fwhm_pixels=30)
    @test all(isfinite, mh)
    @test mh[10, 10, 1] ≈ mp[10, 10, 1] rtol = 0.1

    # Noise whose anisotropy direction differs between the left and right halves:
    # the per-pixel map leaves no preferred direction in either half, a single
    # global covariance does not.
    left = anisotropic_noise(96, 96, 16, MersenneTwister(31); long_axis=0.1)
    right = anisotropic_noise(96, 96, 16, MersenneTwister(32); long_axis=1.4)
    mix(a, b) = cat(a[1:48, :, :], b[49:96, :, :]; dims=1)
    ref = mix(left, right)
    test_noise = mix(anisotropic_noise(96, 96, 16, MersenneTwister(33); long_axis=0.1),
                     anisotropic_noise(96, 96, 16, MersenneTwister(34); long_axis=1.4))
    cmap = gradient_noise_cov_map(ref; fwhm_pixels=24)
    gmap = VGT._whitening_maps(cmap)
    gglob = VGT._gradient_whitening(measure_gradient_noise(ref).cov)
    resultant(half, whiten) = begin
        z = 0.0im; n = 0
        for k in 1:16, j in 2:95, i in (half == 1 ? (8:40) : (57:89))
            gx = (test_noise[i - 1, j, k] - test_noise[i + 1, j, k]) / 2
            gy = (test_noise[i, j - 1, k] - test_noise[i, j + 1, k]) / 2
            wx, wy = whiten(i, j, gx, gy)
            z += cis(2 * atan(wy, wx)); n += 1
        end
        abs(z) / n
    end
    for half in (1, 2)
        r_map = resultant(half, (i, j, gx, gy) -> VGT._whiten_gradient_at(gmap, i, j, gx, gy))
        r_glob = resultant(half, (i, j, gx, gy) -> VGT._whiten_gradient(gglob, gx, gy))
        r_off = resultant(half, (i, j, gx, gy) -> (gx, gy))
        @test r_map < 0.4 * r_off
        @test r_map < 0.5 * r_glob
    end

    # Pipeline: map file read per chunk equals the in-memory run.
    mktempdir() do dir
        nxp, nyp, nvp = 24, 40, 3
        cubep = Float32.(30 .+ [0.2 * i + 0.1 * j for i in 1:nxp, j in 1:nyp, k in 1:nvp] .+ randn(MersenneTwister(4), nxp, nyp, nvp))
        covp = cat([0.3 + 0.002 * i for i in 1:nxp, j in 1:nyp], fill(0.2, nxp, nyp),
                   [0.02 + 0.001 * j for i in 1:nxp, j in 1:nyp]; dims=3)
        covfile = joinpath(dir, "cov.fits")
        write_gradient_noise_cov_map(covfile, covp)
        @test size(read_gradient_noise_cov_map(covfile)) == (nxp, nyp, 3)
        input = joinpath(dir, "cube.fits")
        FITS(input, "w") do f
            write(f, cubep)
        end
        cfgp = VGTConfig(input_file=input, output_dir=joinpath(dir, "out"), noise=1e-6, dn=4, ker=1,
                         ny_chunks=2, min_valid=4, gradient_whiten=true, gradient_noise_cov_file=covfile)
        mkpath(cfgp.output_dir)
        for j in 1:2
            VGT.process_chunk(j, cfgp, nyp, nvp)
        end
        qs, us, = VGT.stitch_chunk_outputs(cfgp, nxp, nyp)
        qd, ud = VGT_QU_mv_2d(cubep, 4, 1, nvp, 1e-6; min_valid=4, gradient_whiten=true,
                              gradient_noise_cov_map=covp)
        @test qs ≈ qd rtol = 1e-10 atol = 1e-8
        @test us ≈ ud rtol = 1e-10 atol = 1e-8
        @test any(!iszero, qs[:, 26:nyp])
        # resolve_execution_config checks the file against the input grid.
        @test VGT.resolve_execution_config(cfgp, nxp, nyp, nvp).gradient_noise_cov_file == covfile
        @test_throws ArgumentError VGT.resolve_execution_config(cfgp, nxp, nyp + 4, nvp)
        missing_cfg = VGT.config_with(cfgp; gradient_noise_cov_file=joinpath(dir, "nope.fits"))
        @test_throws ArgumentError VGT.resolve_execution_config(missing_cfg, nxp, nyp, nvp)
    end

    base = (input_file="c.fits", noise=1.0, dn=4, ker=1, min_valid=4)
    @test_throws ArgumentError VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov=[1.0, 1.0, 0.0],
                                         gradient_noise_cov_file="m.fits")
    @test VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov_file="m.fits").gradient_noise_cov_file == "m.fits"
    @test VGT.run_parameter_hash(VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov_file="m.fits")) !=
          VGT.run_parameter_hash(VGTConfig(; base..., gradient_whiten=true, gradient_noise_cov_file="n.fits"))
end
