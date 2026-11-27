using VGT
using Test

# A NaN voxel must remove only its own channel at that pixel. Before the fix the
# fit centred on it got weight max(NaN, 0) = NaN, which made Q/U (and the
# channel-summed sky intensity) NaN for the whole pixel.

function nan_centre_cube(nx=40, ny=36, nv=6)
    cube = zeros(Float64, nx, ny, nv)
    for k in 1:nv, j in 1:ny, i in 1:nx
        cube[i, j, k] = 20 + 0.3 * (cos(0.5 + 0.2k) * i + sin(0.5 + 0.2k) * j) +
                        0.4 * sin(1.7i + 2.3j + 3k) + 0.3 * cos(0.9i * j / 7 + k)
    end
    return cube
end

function run_products(cube; dn=6, noise=1.0, min_valid=10)
    nv = size(cube, 3)
    q, u, prof = VGT_QU_mv_2d_profiled(cube, dn, 1, nv, noise; min_valid=min_valid)
    fin = VGT.finalize_vgt_products(copy(q), copy(u), VGT.compute_stokes_intensity(cube), 2;
                                    fitted_weight=prof[:fitted_weight])
    return q, u, prof[:fitted_weight], fin
end

@testset "a NaN voxel drops only its own channel" begin
    cube = nan_centre_cube()
    holes = [(20, 18, 3), (21, 18, 3), (12, 25, 5), (30, 10, 1)]
    with_nan = copy(cube)
    with_low = copy(cube)
    for (i, j, k) in holes
        with_nan[i, j, k] = NaN
        # -1 is below the noise threshold: its gradients are masked like a NaN's,
        # and a fit centred on it has Stokes weight max(-1, 0) = 0.
        with_low[i, j, k] = -1.0
    end
    qn, un, wn, fn = run_products(with_nan)
    ql, ul, wl, fl = run_products(with_low)
    # Q/U/W: the NaN channel is skipped, the -1 channel adds exactly zero weight.
    @test qn == ql && un == ul && wn == wl
    @test all(isfinite, qn) && all(isfinite, un)
    for (i, j, _) in holes
        @test isfinite(fn.psi[i, j])
    end
    # The sky intensity of the NaN pixels is the sum of their finite channels.
    I = VGT.compute_stokes_intensity(with_nan)
    for (i, j, k) in holes
        @test I[i, j] ≈ sum(cube[i, j, kk] for kk in axes(cube, 3) if (i, j, kk) ∉ holes)
    end
    @test isequal(fn.psi, fl.psi)
    @test isequal(fn.reliability, fl.reliability)
end

@testset "a pixel blank in every channel stays masked" begin
    cube = nan_centre_cube()
    cube[15, 15, :] .= NaN
    q, u, w, fin = run_products(cube)
    @test isnan(fin.psi[15, 15])
    @test VGT.compute_stokes_intensity(cube)[15, 15] == 0
end

@testset "intensity sum is unchanged where no voxel is blank" begin
    cube = Float32.(nan_centre_cube())
    @test VGT.compute_stokes_intensity(cube) == sum(cube; dims=3)[:, :, 1]
end
