using VGT
using Test
using Random
using Statistics

# Draw `n` independent axial angles whose doubled-angle mean resultant is `rt`
# in expectation: a fraction rt of the angles at 0, the rest uniform.
function mixture_angles(rng, n, rt)
    return [rand(rng) < rt ? 0.0 : π * (rand(rng) - 0.5) for _ in 1:n]
end
resultant(a) = abs(mean(cis.(2 .* a)))

@testset "coherence debiasing" begin
    rng = MersenneTwister(11)
    n, trials = 6, 40_000       # few independent angles per beam: the regime of a 10' beam on VGT maps
    f2 = 1 / n
    for rt in (0.0, 0.3, 0.6, 0.9)
        r = [resultant(mixture_angles(rng, n, rt)) for _ in 1:trials]
        # raw R^2 follows E[R^2] = rt^2 (1 - f^2) + f^2
        @test isapprox(mean(abs2, r), rt^2 * (1 - f2) + f2; atol=0.01)
        # debiased R^2 is unbiased for rt^2 (before the clip at 0, which only matters near rt = 0)
        d2 = (abs2.(r) .- f2) ./ (1 - f2)
        @test isapprox(mean(d2), rt^2; atol=0.01)
        rt >= 0.6 && @test abs(mean(debias_coherence(r, f2)) - rt) < abs(mean(r) - rt)
        # cross product of two independent realizations is unbiased for rt^2
        za = [mean(cis.(2 .* mixture_angles(rng, n, rt))) for _ in 1:trials]
        zb = [mean(cis.(2 .* mixture_angles(rng, n, rt))) for _ in 1:trials]
        @test isapprox(mean(real.(za .* conj.(zb))), rt^2; atol=0.01)
    end

    # floor from a pure-noise angle map: f^2 ~ 1 / (independent angles per beam)
    psi = π .* (rand(MersenneTwister(3), 200, 200) .- 0.5)
    fwhm = 8.0
    nind = 2π * (fwhm / (2 * sqrt(2 * log(2))))^2 * 2   # (sum w)^2 / sum w^2 of a 2D Gaussian = 4 pi s^2
    @test isapprox(coherence_floor2(psi, fwhm), 1 / nind; rtol=0.15)
    @test coherence_floor2(psi, fwhm) < 0.02

    # a uniform angle map has coherence 1 everywhere it is defined, and edges are NaN only where the beam is cut
    z = coherence_phasor(fill(0.3, 50, 50), 6.0)
    @test all(x -> isnan(x) || isapprox(abs(x), 1.0; atol=1e-12), z)
    @test isapprox(angle(z[25, 25]) / 2, 0.3; atol=1e-12)
    @test count(isfinite, z) > 0.8 * length(z)

    # API contract
    @test debias_coherence(0.2, 0.25) == 0.0
    @test isnan(debias_coherence(NaN, 0.25))
    @test debias_coherence(1.0, 0.25) ≈ 1.0
    @test_throws ArgumentError debias_coherence([0.5], 1.0)
    @test_throws ArgumentError cross_coherence(zeros(ComplexF64, 2), zeros(ComplexF64, 3))
    @test cross_coherence([1.0 + 0im], [-1.0 + 0im]) == [0.0]
end
