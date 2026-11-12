using VGT
using Test

headless_abs_diff(a, b) = abs(mod(a - b + pi / 2, pi) - pi / 2)

function synthetic_gradients(angles::AbstractVector{<:Real})
    phi = Float64.(angles)
    return reshape(cos.(phi), :, 1), reshape(sin.(phi), :, 1)
end

@testset "gradient samples stay paired" begin
    # A non-finite Gx drops the whole (Gx, Gy) pair, never just one component.
    grad_x = [1.0 NaN 1.0 1.0]
    grad_y = [0.0 1.0 0.2 0.1]
    @test VGT.fit_window_result(:axial_vonmises, grad_x, grad_y).effective_n == 3
    @test_throws DimensionMismatch VGT.fit_window_result(:axial_vonmises, grad_x, grad_y[:, 1:3])
    @test_throws ArgumentError VGT.fit_window_result(:lsqfit_reference, grad_x, grad_y)
end

@testset "axial von Mises fitter" begin
    true_phi = 0.42
    signal = [true_phi + 0.02 * sin(0.19 * i) for i in 1:200]
    grad_x, grad_y = synthetic_gradients(signal)
    fitted, _ = VGT.fit_axial_vonmises(grad_x, grad_y)
    @test headless_abs_diff(fitted, true_phi) < 0.04

    boundary_phi = pi / 2 - 0.03
    boundary_signal = [isodd(i) ? boundary_phi - 0.01 * sin(0.31 * i) :
                                  -boundary_phi + 0.01 * cos(0.17 * i) for i in 1:200]
    boundary_grad_x, boundary_grad_y = synthetic_gradients(boundary_signal)
    fitted_b, _ = VGT.fit_axial_vonmises(boundary_grad_x, boundary_grad_y)
    @test headless_abs_diff(fitted_b, boundary_phi) < 0.05

    result = VGT.fit_window_result(:axial_vonmises, grad_x, grad_y)
    @test result.fit_ok
    @test headless_abs_diff(result.angle, true_phi) < 0.04
    @test 0.0 <= result.rbar2 <= 1.0

    perfect_grad_x, perfect_grad_y = synthetic_gradients(fill(true_phi, 64))
    perfect_result = VGT.fit_window_result(:axial_vonmises, perfect_grad_x, perfect_grad_y)
    @test perfect_result.fit_ok
    @test perfect_result.rbar2 ≈ 1.0

    invalid_grad_x = [1.0 Inf NaN]
    invalid_grad_y = [0.0 1.0 1.0]
    invalid_result = VGT.fit_window_result(:axial_vonmises, invalid_grad_x, invalid_grad_y)
    @test invalid_result.effective_n == 1
    @test headless_abs_diff(invalid_result.angle, 0.0) < 1e-12

    low_conf_grad_x, low_conf_grad_y = synthetic_gradients(vcat(fill(0.0, 74), fill(pi / 2, 26)))
    low_conf_result = VGT.fit_window_result(:axial_vonmises, low_conf_grad_x, low_conf_grad_y)
    @test low_conf_result.fit_ok
    @test low_conf_result.rbar2 < 0.5
    @test low_conf_result.low_confidence

    diffuse_grad_x, diffuse_grad_y = synthetic_gradients(vcat(fill(0.0, 50), fill(pi / 2, 50)))
    diffuse_result = VGT.fit_window_result(:axial_vonmises, diffuse_grad_x, diffuse_grad_y)
    @test !diffuse_result.fit_ok
    @test diffuse_result.failure_reason == :degenerate_resultant
    @test diffuse_result.low_confidence
    @test diffuse_result.rbar2 < 1e-12
end
