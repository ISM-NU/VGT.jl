using VGT
using Test
using Random

# Grid-locked angle response of the gradient stencils (docs/validation-angle-calibration.md).

# Gradient angle (rad, mod pi) that a stencil returns on the plane wave cos(k (x cos t + y sin t)).
# Every interior pixel gives the same axial angle, so one pixel suffices.
function plane_wave_angle(mode, k, t; n=9)
    img = [2 + cos(k * (i * cos(t) + j * sin(t)) + 0.3) for i in 1:n, j in 1:n]
    gx, gy = VGT._stencil_gradient(img, 5, 5, mode)
    return atan(gy, gx)
end

axial_deg(a) = rad2deg(mod(a + pi / 2, pi) - pi / 2)

@testset "stencil angle response on plane waves" begin
    t = deg2rad(22.5)  # where the sin(4 theta) grid term is largest
    for k in (0.5, 1.0)
        exact(a) = atan(sin(k * sin(t)) * ((1 - 2a) + 2a * cos(k * cos(t))),
                        sin(k * cos(t)) * ((1 - 2a) + 2a * cos(k * sin(t))))
        for (mode, a) in ((:central_checked, 0.0), (:sobel_checked, 0.25), (:isotropic_checked, 1 / 6))
            @test axial_deg(plane_wave_angle(mode, k, t) - exact(a)) ≈ 0 atol = 1e-9
        end
    end
    # Leading terms: central +k^2/24 sin(4t), Sobel -k^2/48 sin(4t), isotropic O(k^4).
    k = 0.5
    @test axial_deg(plane_wave_angle(:central_checked, k, t) - t) ≈ rad2deg(k^2 / 24) rtol = 0.05
    @test axial_deg(plane_wave_angle(:sobel_checked, k, t) - t) ≈ -rad2deg(k^2 / 48) rtol = 0.05
    @test abs(axial_deg(plane_wave_angle(:isotropic_checked, k, t) - t)) < 0.01
    @test abs(axial_deg(plane_wave_angle(:isotropic_checked, 1.0, t) - t)) < 0.1
    # On the axes and diagonals no stencil has a bias.
    for mode in (:central_checked, :sobel_checked, :isotropic_checked), tt in (0.0, pi / 4, pi / 2)
        @test abs(axial_deg(plane_wave_angle(mode, 1.0, tt) - tt)) < 1e-9
    end
end

# The full kernel (stencil, mask, window fit, Q/U smoothing) has the symmetry of the square grid, so it cannot
# produce a uniform angle offset; any grid bias is odd and 90-degree periodic (a sin(4 n theta) series).
@testset "VGT kernel is equivariant under grid mirrors" begin
    rng = MersenneTwister(11)
    nx, nv, dn, ker = 48, 3, 8, 2
    cube = 10 .+ randn(rng, nx, nx, nv)
    psi_of(c, mode; w=nothing) = begin
        # the noise threshold applies to the weight cube when one is given (all ones here)
        q, u = VGT.VGT_QU_mv_2d(c, dn, 1, nv, w === nothing ? 7.0 : 0.5; min_valid=20, gradient_mode=mode, weight_cube=w)
        intensity = dropdims(sum(w === nothing ? c : w; dims=3); dims=3)
        VGT.finalize_vgt_products(q, u, intensity, ker).psi
    end
    for mode in (:central_checked, :sobel_checked, :isotropic_checked)
        # Transpose (mirror about the diagonal): exact, a field angle t becomes pi/2 - t.
        p = psi_of(cube, mode)
        pt = psi_of(permutedims(cube, (2, 1, 3)), mode)
        ok = isfinite.(p) .& isfinite.(permutedims(pt))
        @test count(ok) > 500
        @test maximum(abs.(axial_deg.((pi / 2 .- p) .- permutedims(pt))[ok])) < 1e-8
        # Flip of axis 1 with uniform Stokes weights: t becomes -t, and the even window moves by one pixel
        # (pixel p of the flipped map matches pixel nx - p). The default weight is the window's own pixel,
        # half a pixel off the window centre, so the flip is exact only with uniform weights.
        w = ones(size(cube))
        pf = psi_of(reverse(cube; dims=1), mode; w=reverse(w; dims=1))
        p1 = psi_of(cube, mode; w=w)
        a = p1[1:nx - 1, :]
        b = pf[nx - 1:-1:1, :]
        ok = isfinite.(a) .& isfinite.(b)
        @test count(ok) > 500
        @test maximum(abs.(axial_deg.(-a .- b)[ok])) < 1e-8
    end
end
