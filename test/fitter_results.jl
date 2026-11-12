using VGT
using Test

headless_abs_diff(a, b) = abs(mod(a - b + pi / 2, pi) - pi / 2)

function synthetic_gradients_for_results(angles::AbstractVector{<:Real})
    phi = Float64.(angles)
    return reshape(cos.(phi), :, 1), reshape(sin.(phi), :, 1)
end


@testset "structured fitter result failure path" begin
    grad_x = fill(NaN, 8, 8)
    grad_y = fill(NaN, 8, 8)
    result = VGT.fit_window_result(:axial_vonmises, grad_x, grad_y)
    @test result isa VGT.WindowFitResult
    @test !result.fit_ok
    @test isnan(result.angle)
    @test isnan(result.rbar2)
    @test result.effective_n == 0
    @test result.failure_reason == :empty_sample
end

@testset "fitter confidence flag follows rbar2" begin
    concentrated_grad_x, concentrated_grad_y = synthetic_gradients_for_results(fill(0.25, 128))
    concentrated = VGT.fit_window_result(:axial_vonmises, concentrated_grad_x, concentrated_grad_y)
    @test concentrated.fit_ok
    @test !concentrated.low_confidence

    bimodal_angles = vcat(fill(-0.5, 64), fill(0.5, 64))
    bimodal_grad_x, bimodal_grad_y = synthetic_gradients_for_results(bimodal_angles)
    bimodal = VGT.fit_window_result(:axial_vonmises, bimodal_grad_x, bimodal_grad_y)
    @test bimodal.fit_ok
    @test bimodal.rbar2 < concentrated.rbar2
end
