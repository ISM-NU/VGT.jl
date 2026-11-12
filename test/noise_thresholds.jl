using VGT
using Test

@testset "robust negative noise estimator" begin
    sigma = 0.42
    sample = fill(-sigma * VGT.NEGATIVE_MAD_NORMAL, 10, 10)

    @test VGT.robust_negative_noise(sample) ≈ sigma
    @test isnan(VGT.robust_negative_noise(fill(-1.0, 9, 10)))
end

@testset "dynamic Noise builders" begin
    nx, ny, nv = 20, 20, 4
    cube = zeros(Float64, nx, ny, nv)
    @inbounds for k in 1:nv, j in 1:ny, i in 1:nx
        sigma = 0.1 * k
        cube[i, j, k] = isodd(i + j + k) ? -VGT.NEGATIVE_MAD_NORMAL * sigma : sigma
    end

    channel = VGT.estimate_channel_noise(cube; fallback=0.1, min_count=1)
    @test channel[4] > channel[1]

    spatial = VGT.estimate_spatial_noise(cube; fallback=0.1, tile_size=10, min_count=1)
    @test size(spatial) == (nx, ny)
    @test all(isfinite, spatial)

    scalar = VGT.build_noise_cube(cube, 0.1, :scalar; tile_size=10)
    @test scalar == 0.1
    @test VGT.build_noise_cube(cube, 0.1, :scalar; tile_size=10, multiplier=3.0) ≈ 0.3

    channel_cube = VGT.build_noise_cube(cube, 0.1, :channel; tile_size=10)
    @test size(channel_cube) == size(cube)
    @test channel_cube[1, 1, 4] > channel_cube[1, 1, 1]
    channel_cube_3sigma = VGT.build_noise_cube(
        cube, 0.1, :channel; tile_size=10, multiplier=3.0
    )
    @test channel_cube_3sigma ≈ 3.0 .* channel_cube

    separable = VGT.build_noise_cube(cube, 0.1, :separable; tile_size=10)
    @test size(separable) == size(cube)
    @test all(isfinite, separable)
    separable_3sigma = VGT.build_noise_cube(
        cube, 0.1, :separable; tile_size=10, multiplier=3.0
    )
    @test separable_3sigma ≈ 3.0 .* separable

    model = VGT.build_noise_model(cube, 0.1, :separable; tile_size=10)
    @test model.mode == :separable
    chunk_noise = VGT.noise_for_chunk(model, nx, 5:12, nv)
    @test chunk_noise == separable[:, 5:12, :]
    @test VGT.noise_for_chunk(model, nx, 5:12, nv; multiplier=2.5) ≈
          2.5 .* chunk_noise
    @test_throws ArgumentError VGT.noise_for_chunk(model, nx, 5:12, nv; multiplier=0.0)

    auto_model = VGT.build_noise_model(cube, 0.1, :auto; tile_size=10)
    @test auto_model.mode in (:scalar, :channel, :spatial, :separable)
end
