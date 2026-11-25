using VGT
using FITSIO
using Random
using Test

function make_linear_cube(nx::Int, ny::Int, nv::Int)
    cube = zeros(Float64, nx, ny, nv)
    @inbounds for k in 1:nv, j in 1:ny, i in 1:nx
        cube[i, j, k] = 10.0 * k + 2.0 * i + 3.0 * j
    end
    return cube
end

function reference_finalize_vgt_products(q_accum, u_accum, intensity_accum, fitted_weight,
                                         ker::Int, zero_stokes_policy::Symbol,
                                         output_support_mode::Symbol)
    invalid_intensity = .!isfinite.(intensity_accum) .| (intensity_accum .<= 0)
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))

    q_before_smoothing = copy(q_accum)
    u_before_smoothing = copy(u_accum)
    if zero_stokes_policy === :nan_before_smoothing
        q_before_smoothing[q_accum .== 0] .= NaN
        u_before_smoothing[u_accum .== 0] .= NaN
    end
    q_before_smoothing[invalid_intensity] .= NaN
    u_before_smoothing[invalid_intensity] .= NaN
    q_smooth = VGT.imfilter_gaussian(q_before_smoothing, ker)
    u_smooth = VGT.imfilter_gaussian(u_before_smoothing, ker)
    q_smooth[invalid_intensity] .= NaN
    u_smooth[invalid_intensity] .= NaN
    if output_support_mode === :direct_fit_support
        direct_invalid = (q_accum .== 0.0) .& (u_accum .== 0.0)
        q_smooth[direct_invalid] .= NaN
        u_smooth[direct_invalid] .= NaN
    end
    psi = 0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)

    # Reliability: Q, U and the fitted weight W smoothed over one joint mask.
    joint_invalid = copy(invalid_intensity)
    zero_stokes_policy === :nan_before_smoothing && (joint_invalid .|= fitted_weight .== 0)
    q_r, u_r, w_r = copy(q_accum), copy(u_accum), copy(fitted_weight)
    for a in (q_r, u_r, w_r)
        a[joint_invalid] .= NaN
    end
    q_r = VGT.imfilter_gaussian(q_r, ker)
    u_r = VGT.imfilter_gaussian(u_r, ker)
    weight_smooth = VGT.imfilter_gaussian(w_r, ker)
    q_r[isnan.(psi)] .= NaN
    u_r[isnan.(psi)] .= NaN
    reliability = VGT.compute_reliability(q_r, u_r, weight_smooth)
    return (; q_smooth, u_smooth, weight_smooth, psi, reliability)
end

function test_same_finalized_products(actual, expected)
    @test isequal(actual.q_smooth, expected.q_smooth)
    @test isequal(actual.u_smooth, expected.u_smooth)
    @test isequal(actual.weight_smooth, expected.weight_smooth)
    @test isequal(actual.psi, expected.psi)
    @test isequal(actual.reliability, expected.reliability)
end

@testset "dn and window validation" begin
    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=21,
        ker=4,
        ny_chunks=1,
        min_valid=100,
    )

    @test_throws ArgumentError VGT.VGTConfig(
        input_file="dummy.fits",
        output_dir="out",
        keep_intermediate_chunks=false,
        noise=1.0,
        dn=20,
        ker=4,
        ny_chunks=1,
        min_valid=401,
    )

    cube = make_linear_cube(32, 32, 3)
    @test_throws ArgumentError VGT.VGT_QU_mv_2d(cube, 21, 1, 3, 0.1)
    @test_throws ArgumentError VGT.VGT_QU_mv_2d(cube, 40, 1, 3, 0.1)
end

@testset "array Noise preserves scalar behavior when constant" begin
    cube = make_linear_cube(32, 32, 3)
    dn = 8
    noise = 0.1
    min_valid = 16

    q_scalar, u_scalar = VGT.VGT_QU_mv_2d(cube, dn, 1, 3, noise; min_valid=min_valid)
    q_array, u_array = VGT.VGT_QU_mv_2d(cube, dn, 1, 3, fill(noise, size(cube)); min_valid=min_valid)

    @test q_array == q_scalar
    @test u_array == u_scalar
    @test_throws ArgumentError VGT.VGT_QU_mv_2d(cube, dn, 1, 3, fill(noise, 4, 4, 3); min_valid=min_valid)
end

@testset "central checked gradient rejects incomplete stencils" begin
    cube = make_linear_cube(12, 12, 1)
    dn = 4
    noise = 0.1

    q_central, u_central = VGT.VGT_QU_mv_2d(
        cube, dn, 1, 1, noise;
        min_valid=4,
        gradient_mode=:central_checked,
    )
    q_backward, u_backward = VGT.VGT_QU_mv_2d(
        cube, dn, 1, 1, noise;
        min_valid=4,
        gradient_mode=:nan_edge_backward,
    )
    @test count(!=(0.0), q_central) <= count(!=(0.0), q_backward)
    @test count(!=(0.0), u_central) <= count(!=(0.0), u_backward)
end

@testset "min_valid controls window support" begin
    sparse = zeros(Float64, 8, 8, 1)
    sparse[2:3, 2:3, 1] .= 10.0
    q_min, u_min = VGT.VGT_QU_mv_2d(
        sparse, 4, 1, 1, 1.0;
        min_valid=4,
        fitter_mode=:axial_vonmises,
        gradient_mode=:nan_edge_backward,
    )
    q_above, u_above = VGT.VGT_QU_mv_2d(
        sparse, 4, 1, 1, 1.0;
        min_valid=5,
        fitter_mode=:axial_vonmises,
        gradient_mode=:nan_edge_backward,
    )
    @test count(!=(0.0), q_min) + count(!=(0.0), u_min) >
          count(!=(0.0), q_above) + count(!=(0.0), u_above)
end

@testset "output support mode controls smoothed halo" begin
    q = zeros(Float64, 9, 9)
    u = zeros(Float64, 9, 9)
    i = ones(Float64, 9, 9)
    q[5, 5] = 1.0
    u[5, 5] = 0.5

    halo = VGT.finalize_vgt_products(
        q, u, i, 1, :nan_before_smoothing, :smoothed_stokes_support
    )
    no_halo = VGT.finalize_vgt_products(
        q, u, i, 1, :nan_before_smoothing, :direct_fit_support
    )

    @test count(isfinite, halo.psi) > count(isfinite, no_halo.psi)
    @test count(isfinite, halo.psi) == 25
    @test count(isfinite, no_halo.psi) == 1

    zero_smoothed = VGT.finalize_vgt_products(
        q, u, i, 1, :keep_zero_before_smoothing, :smoothed_stokes_support
    )
    @test count(isfinite, zero_smoothed.psi) == 81
end

@testset "finalization masks zero-intensity sky like legacy path" begin
    q = ones(Float64, 9, 9)
    u = ones(Float64, 9, 9)
    i = ones(Float64, 9, 9)

    i[:, 1:3] .= 0.0
    i[1:2, :] .= NaN

    finalized = VGT.finalize_vgt_products(q, u, i, 1; fitted_weight=fill(2.0, 9, 9))

    @test all(isnan, finalized.q_smooth[:, 1:3])
    @test all(isnan, finalized.u_smooth[:, 1:3])
    @test all(isnan, finalized.q_smooth[1:2, :])
    @test all(isnan, finalized.u_smooth[1:2, :])
    @test all(isnan, finalized.psi[:, 1:3])
    @test all(isnan, finalized.reliability[:, 1:3])
    @test all(isnan, finalized.psi[1:2, :])
    @test all(isnan, finalized.reliability[1:2, :])
    @test isfinite(finalized.q_smooth[9, 9])
    @test isfinite(finalized.u_smooth[9, 9])
    @test isfinite(finalized.psi[9, 9])
    @test isfinite(finalized.reliability[9, 9])
end

@testset "finalization loop masks preserve vectorized behavior" begin
    q = reshape(Float64[
        0, 1, 2, 3, 4,
        5, 0, 7, 8, 9,
        10, 11, 0, 13, 14,
        15, 16, 17, 0, 19,
        20, 21, 22, 23, 0,
    ], 5, 5)
    u = reshape(Float64[
        0, -1, -2, -3, -4,
        -5, 0, -7, -8, -9,
        -10, -11, 0, -13, -14,
        -15, -16, -17, 0, -19,
        -20, -21, -22, -23, 0,
    ], 5, 5)
    i = ones(Float64, 5, 5)
    i[1, 1] = 0.0
    i[2, 3] = -1.0
    i[4, 2] = NaN
    i[5, 5] = 3.0
    w = 1.5 .* hypot.(q, u)

    # u2 is exactly zero at a fitted pixel (W > 0), so Q and U have different
    # masks and the reliability uses its own joint-mask smoothing.
    u2 = copy(u)
    u2[3, 2] = 0.0
    for (uu, label) in ((u, "shared masks"), (u2, "Q/U masks differ")),
        zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing),
        output_support_mode in (:smoothed_stokes_support, :direct_fit_support)

        actual = VGT.finalize_vgt_products(
            q, uu, i, 1, zero_stokes_policy, output_support_mode; fitted_weight=w
        )
        expected = reference_finalize_vgt_products(
            q, uu, i, w, 1, zero_stokes_policy, output_support_mode
        )
        test_same_finalized_products(actual, expected)
    end
end

@testset "finalization does not mutate accumulated products" begin
    q = reshape(Float64[0, 1, 2, 3, 4, 5, 0, 7, 8], 3, 3)
    u = reshape(Float64[0, -1, -2, -3, -4, -5, 0, -7, -8], 3, 3)
    i = ones(Float64, 3, 3)
    i[1, 2] = 0.0
    i[3, 1] = NaN
    q_before = copy(q)
    u_before = copy(u)
    i_before = copy(i)

    VGT.finalize_vgt_products(q, u, i, 1,
                                    :nan_before_smoothing, :direct_fit_support)

    @test isequal(q, q_before)
    @test isequal(u, u_before)
    @test isequal(i, i_before)
end

@testset "finalization excludes invalid intensity before smoothing" begin
    q = ones(Float64, 9, 9)
    u = ones(Float64, 9, 9)
    i = ones(Float64, 9, 9)
    q[5, 5] = 100.0
    u[5, 5] = 100.0
    i[5, 5] = 0.0

    w = fill(2.0, 9, 9)
    w[5, 5] = 200.0
    strict = VGT.finalize_vgt_products(
        q, u, i, 1, :keep_zero_before_smoothing, :smoothed_stokes_support; fitted_weight=w
    )

    @test isnan(strict.psi[5, 5])
    # The masked pixel's large Q/U must not leak into its neighbours.
    @test strict.q_smooth[5, 6] ≈ 1.0
    @test strict.u_smooth[5, 6] ≈ 1.0
    @test strict.weight_smooth[5, 6] ≈ 2.0
end

@testset "weight_cube decouples mask and weights from the gradient cube" begin
    nx, ny, nv = 24, 24, 3
    grad = Float64[sin(0.4i) + cos(0.3j) + 0.1k for i in 1:nx, j in 1:ny, k in 1:nv]
    q0, u0 = VGT_QU_mv_2d(grad, 8, 1, nv, -10.0; min_valid=20)
    q1, u1 = VGT_QU_mv_2d(grad, 8, 1, nv, -10.0; min_valid=20, weight_cube=grad)
    @test isequal(q0, q1) && isequal(u0, u1)
    # Negative gradient-cube values are kept when the weight cube passes the mask,
    # and the angles are those of the gradient cube with the weight cube's weights.
    w = fill(2.0, nx, ny, nv)
    q2, u2 = VGT_QU_mv_2d(grad .- 5, 8, 1, nv, 1.0; min_valid=20, weight_cube=w)
    qr, ur = VGT_QU_mv_2d(grad .+ 100, 8, 1, nv, 1.0; min_valid=20)
    m = (qr .!= 0) .& (q2 .!= 0)
    @test count(m) > 0
    @test atan.(u2[m], q2[m]) ≈ atan.(ur[m], qr[m])
    @test_throws ArgumentError VGT_QU_mv_2d(grad, 8, 1, nv, 1.0; weight_cube=zeros(2, 2, 2))
end

# Durable contract (decision 2026-11-25 in docs/decisions.md): reliability is
# sqrt(Q^2+U^2)/W with W the summed Stokes weight of the fitted windows, so it
# never exceeds 1 beyond rounding, whatever the fit coverage, noise or chunking.
@testset "psi_reliability stays in [0, 1]" begin
    nx, ny, nv = 48, 44, 8
    rng = MersenneTwister(7)
    cube = zeros(Float32, nx, ny, nv)
    for k in 1:nv, j in 1:ny, i in 1:nx
        # Filaments and a blob that moves with channel over a background below
        # the noise threshold, so each pixel is fitted in only some channels;
        # the background straddles zero, so negative centre pixels occur too.
        fil = 3f0 * exp(-((0.6i + 0.8j - 30 - k) / 4)^2) + 2.5f0 * exp(-((0.9i - 0.4j - 5 + k) / 3)^2)
        blob = 4f0 * exp(-((i - 10 - 3k)^2 + (j - 35 + 2k)^2) / (2 * 8^2))
        cube[i, j, k] = fil + blob + 0.3f0 * randn(rng, Float32)
    end
    @test any(<(0), cube)
    mktempdir() do dir
        path = joinpath(dir, "c.fits")
        FITS(f -> write(f, cube), path, "w")
        base = (; input_file=path, keep_intermediate_chunks=false, noise=0.6, noise_mode=:scalar,
                dn=8, ker=2, min_valid=16)
        for (name, extra) in (("one", (; ny_chunks=1)), ("split", (; ny_chunks=3)),
                              ("bin", (; ny_chunks=1, output_bin=3)),
                              ("direct", (; ny_chunks=1, output_support_mode=:direct_fit_support)))
            run = run_pipeline(VGTConfig(; base..., extra..., output_dir=joinpath(dir, name)))
            rel = filter(isfinite, run.reliability)
            @test length(rel) > 0.3 * count(isfinite, run.psi)
            @test all(>=(0), rel)
            @test maximum(rel) <= 1 + 1e-6          # Float32 output
            @test all(run.W .>= hypot.(run.Q, run.U) .* (1 - 1e-12))
            # Partial coverage: some sky with intensity has no fitted weight.
            name == "one" && @test any((run.W .== 0) .& (run.I .> 0))
        end
    end
    # Exact-zero U at fitted pixels (an axis-aligned fitted angle) splits the Q
    # and U masks; the reliability must still use one joint mask.
    q = fill(1.0, 12, 12)
    u = zeros(12, 12)
    u[7:12, :] .= 1.0
    q[1:3, 1:3] .= 0.0
    u[1:3, 1:3] .= 0.0
    w = hypot.(q, u)
    fin = VGT.finalize_vgt_products(q, u, ones(12, 12), 2; fitted_weight=w)
    rel = filter(isfinite, fin.reliability)
    @test !isempty(rel)
    @test maximum(rel) <= 1 + 1e-12
end

@testset "negative centre intensity gives zero Stokes weight" begin
    nx, ny, nv = 24, 24, 2
    cube = Float64[0.5 * i + 0.2 * j + k for i in 1:nx, j in 1:ny, k in 1:nv]
    weight = fill(2.0, nx, ny, nv)
    weight[12, 12, 1] = -3.0                       # negative window centre in channel 1
    # Threshold below every weight, so the window still fits; only its weight changes.
    q, u, prof = VGT_QU_mv_2d_profiled(cube, 8, 1, nv, -10.0; min_valid=20, weight_cube=weight)
    w = prof[:fitted_weight]
    @test prof[:fit_support][12, 12] == 1.0        # both channels fitted
    @test w[12, 12] == 2.0                         # channel 1 contributes zero, not -3
    @test w[13, 12] == 4.0
    # Same orientation as a neighbour: the clamped window does not flip the vote.
    @test atan(u[12, 12], q[12, 12]) ≈ atan(u[13, 12], q[13, 12]) atol = 1e-6
    @test hypot(q[12, 12], u[12, 12]) ≈ 2.0 rtol = 1e-6
end
