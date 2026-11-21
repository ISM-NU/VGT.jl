using VGT
using Test
using Random
using Statistics

# Monte Carlo noise uncertainty of the angle map (src/science/uncertainty.jl,
# docs/validation-resolution-uncertainty.md).

lagcor(a, dx, dy) = cor(vec(a[1:end - dx, 1:end - dy, :]), vec(a[1 + dx:end, 1 + dy:end, :]))

@testset "correlated_noise" begin
    rng = Xoshiro(11)
    z = correlated_noise(rng, 128, 128, 6; fwhm_major_px=3.0)
    @test std(z) ≈ 1 rtol = 0.05
    s = 3 / sqrt(8log(2))
    @test lagcor(z, 1, 0) ≈ exp(-1 / (4s^2)) atol = 0.03   # Gaussian-beam ACF at lag 1
    @test lagcor(z, 0, 1) ≈ exp(-1 / (4s^2)) atol = 0.03
    e = correlated_noise(rng, 128, 128, 6; fwhm_major_px=4.0, fwhm_minor_px=2.0, angle=0.0)
    @test lagcor(e, 1, 0) > lagcor(e, 0, 1) + 0.1         # longer correlation along the major axis (axis 1)
    sp = correlated_noise(rng, 64, 64, 40; fwhm_major_px=0.0, spectral_lag1=0.23)
    @test cor(vec(sp[:, :, 1:end - 1]), vec(sp[:, :, 2:end])) ≈ 0.23 atol = 0.03
    m = fill(2.0, 64, 64)
    m[1:32, :] .= 0.5
    w = correlated_noise(rng, 64, 64, 20; fwhm_major_px=3.0, sigma=m)
    @test std(w[33:end, :, :]) / std(w[1:32, :, :]) ≈ 4 rtol = 0.1
end

@testset "circular_spread" begin
    a = fill(0.3, 4, 4)
    r = circular_spread([a, a, a])
    @test all(r.sigma .≈ 0)
    @test all(r.mean_angle .≈ 0.3)
    # Wrap: +-(pi/2 - 0.01) are 0.02 rad apart as orientations.
    r = circular_spread([fill(pi / 2 - 0.01, 2, 2), fill(-pi / 2 + 0.01, 2, 2)])
    @test all(isapprox.(r.sigma, 0.01; rtol=0.01))
    r = circular_spread([fill(NaN, 2, 2), fill(0.1, 2, 2)])
    @test all(isnan.(r.sigma)) && all(r.n .== 1)
end

@testset "psi_noise_mc" begin
    rng = Xoshiro(5)
    n, nv = 64, 3
    # Stripes along axis 2 (structures along axis 2, gradients along axis 1) with small-scale texture.
    cube = [10 + 3sin(0.7i + 0.2v) + 0.5sin(0.31j + 1.3i) for i in 1:n, j in 1:n, v in 1:nv]
    q, u = VGT_QU_mv_2d(cube, 8, 1, nv, 1.0; min_valid=16)
    ref = finalize_vgt_products(q, u, dropdims(sum(cube; dims=3); dims=3), 1).psi
    lo = psi_noise_mc(cube, 8, 1, 1.0; realizations=6, noise_sigma=0.05, fwhm_major_px=3.0, min_valid=16, rng=Xoshiro(1))
    hi = psi_noise_mc(cube, 8, 1, 1.0; realizations=6, noise_sigma=0.2, fwhm_major_px=3.0, min_valid=16, rng=Xoshiro(1))
    ok = isfinite.(ref)
    @test isequal(lo.psi, ref)                       # the returned angle is the plain VGT map
    @test median(lo.sigma_psi[ok]) < deg2rad(2)
    # Linear regime: four times the noise, about four times the spread.
    @test median(hi.sigma_psi[ok]) / median(lo.sigma_psi[ok]) ≈ 4 rtol = 0.3
    # scale halves the added noise and divides the spread back.
    half = psi_noise_mc(cube, 8, 1, 1.0; realizations=6, noise_sigma=0.2, fwhm_major_px=3.0, min_valid=16,
                        rng=Xoshiro(1), scale=0.5)
    @test median(half.sigma_psi[ok]) / median(hi.sigma_psi[ok]) ≈ 1 rtol = 0.3
    @test all(lo.n[ok] .== 6)
end
