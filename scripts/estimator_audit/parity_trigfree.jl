# Does the trig-free doubled-angle mean (est_weighted p=0) reproduce production Q/U?
ENV["AUDIT_SEEDS"] = "0:-1"; ENV["AUDIT_OUT"] = tempname()
include(joinpath(@__DIR__, "run_audit.jl"))
for (label, b) in (("crossing_filaments", MW.build_crossing_filaments), ("sparse_co_like", MW.build_sparse_co_like))
    cube, truth, noise = b(; nx=96, ny=96)
    q0, u0, _ = VGT_QU_mv_2d_profiled(cube, 12, 1, size(cube, 3), noise; min_valid=50, gradient_mode=:central_checked)
    patch_kernel!()
    @eval Main estimate(mode::Symbol, gx, gy) = mode === :tf ? est_weighted(gx, gy, 0) : vm_production(gx, gy)
    q1, u1, _ = Base.invokelatest(VGT_QU_mv_2d_profiled, cube, 12, 1, size(cube, 3), noise; min_valid=50, fitter_mode=:tf, gradient_mode=:central_checked)
    a0 = 0.5 .* atan.(u0, q0); a1 = 0.5 .* atan.(u1, q1)
    d = filter(isfinite, rad2deg.(abs.(headless_angle_difference.(a0, a1))))
    @printf("%-20s max |dQ| %.2e  max |dpsi| %.2e deg  zero-Stokes pixels equal: %s\n", label,
            maximum(abs.(q0 .- q1)), maximum(d), ((q0 .== 0) == (q1 .== 0)))
end
