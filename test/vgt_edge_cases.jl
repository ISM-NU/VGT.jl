using VGT
using Test

function make_linear_cube(nx::Int, ny::Int, nv::Int)
    cube = zeros(Float64, nx, ny, nv)
    @inbounds for k in 1:nv, j in 1:ny, i in 1:nx
        cube[i, j, k] = 10.0 * k + 2.0 * i + 3.0 * j
    end
    return cube
end

function reference_finalize_vgt_products(q_accum, u_accum, intensity_accum,
                                         ker::Int, zero_stokes_policy::Symbol,
                                         output_support_mode::Symbol)
    invalid_intensity = .!isfinite.(intensity_accum) .| (intensity_accum .<= 0)
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))

    q_before_smoothing = copy(q_accum)
    u_before_smoothing = copy(u_accum)
    intensity_before_smoothing = copy(intensity_accum)
    if zero_stokes_policy === :nan_before_smoothing
        q_before_smoothing[q_accum .== 0] .= NaN
        u_before_smoothing[u_accum .== 0] .= NaN
    end
    q_before_smoothing[invalid_intensity] .= NaN
    u_before_smoothing[invalid_intensity] .= NaN
    intensity_before_smoothing[invalid_intensity] .= NaN
    q_smooth = VGT.imfilter_gaussian(q_before_smoothing, ker)
    u_smooth = VGT.imfilter_gaussian(u_before_smoothing, ker)
    intensity_smooth = VGT.imfilter_gaussian(intensity_before_smoothing, ker)
    q_smooth[invalid_intensity] .= NaN
    u_smooth[invalid_intensity] .= NaN
    if output_support_mode === :direct_fit_support
        direct_invalid = (q_accum .== 0.0) .& (u_accum .== 0.0)
        q_smooth[direct_invalid] .= NaN
        u_smooth[direct_invalid] .= NaN
    end
    psi = 0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)
    reliability = VGT.compute_reliability(q_smooth, u_smooth, intensity_smooth)
    return (; q_smooth, u_smooth, intensity_smooth, psi, reliability)
end

function test_same_finalized_products(actual, expected)
    @test isequal(actual.q_smooth, expected.q_smooth)
    @test isequal(actual.u_smooth, expected.u_smooth)
    @test isequal(actual.intensity_smooth, expected.intensity_smooth)
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

    finalized = VGT.finalize_vgt_products(q, u, i, 1)

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

    for zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing),
        output_support_mode in (:smoothed_stokes_support, :direct_fit_support)

        actual = VGT.finalize_vgt_products(
            q, u, i, 1, zero_stokes_policy, output_support_mode
        )
        expected = reference_finalize_vgt_products(
            q, u, i, 1, zero_stokes_policy, output_support_mode
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

    strict = VGT.finalize_vgt_products(
        q, u, i, 1, :keep_zero_before_smoothing, :smoothed_stokes_support
    )

    @test isnan(strict.psi[5, 5])
    # The masked pixel's large Q/U must not leak into its neighbours.
    @test strict.q_smooth[5, 6] ≈ 1.0
    @test strict.u_smooth[5, 6] ≈ 1.0
    @test strict.intensity_smooth[5, 6] ≈ 1.0
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
