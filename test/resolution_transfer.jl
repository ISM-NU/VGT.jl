using Test

# Effective-resolution estimator of scripts/validate_resolution_uncertainty.jl (part tf): the half-transfer FWHM must
# recover a Gaussian response and must ignore the sidelobes of a box (sinc) response, which bias a Gaussian fit.

const RES = Module(:ResolutionTransferTest)
Core.eval(RES, :(using Test))
Core.eval(RES, :(include(p) = Base.include($RES, p)))   # the script includes its generator
Base.include(RES, joinpath(@__DIR__, "..", "scripts", "validate_resolution_uncertainty.jl"))

@testset "half-transfer FWHM" begin
    ks = 2pi ./ (512 ./ [2, 4, 6, 8, 12, 16, 24, 32, 40])     # the part tf wavelength grid
    for F in (8.0, 16.0, 24.0)
        s = F / sqrt(8log(2))
        lh, fw = RES.half_transfer(ks, exp.(-ks .^ 2 .* s^2 ./ 2))
        @test fw ≈ F rtol = 0.06                            # log-linear interpolation on a coarse grid
        @test lh ≈ F * pi / (2log(2)) rtol = 0.06
    end
    # Box of width dn: T = |sinc(k dn / 2)|, first zero at lambda = dn, sidelobes up to 0.22 below it.
    for dn in (20, 32, 48)
        T = abs.(sin.(ks .* dn ./ 2) ./ (ks .* dn ./ 2))
        lh, fw = RES.half_transfer(ks, T)
        @test lh ≈ pi * dn / 1.8955 rtol = 0.06             # sinc(x) = 0.5 at x = 1.8955
        @test fw > RES.gaussfit_fwhm(ks, T) || dn < 32      # the Gaussian fit is pulled low by the sidelobes
    end
    # Half the plateau: a Gaussian scaled by 0.94 crosses 0.47 where the unscaled one crosses 0.5.
    s16 = 16 / sqrt(8log(2))
    @test RES.half_transfer(ks, 0.94 .* exp.(-ks .^ 2 .* s16^2 ./ 2); level=0.47)[2] ≈ 16 rtol = 0.06
    @test all(isnan, RES.half_transfer([0.1, 0.2], [0.9, 0.8]))
    @test_throws ArgumentError RES.half_transfer([0.2, 0.1], [0.9, 0.1])
end

@testset "fit_tf writes the half-transfer FWHM where part sweep reads it" begin
    mktempdir() do dir
        runs, out = joinpath(dir, "tf_runs.csv"), joinpath(dir, "tf_fit.csv")
        ks = 2pi ./ (512 ./ [2, 4, 6, 8, 12, 16, 24, 32, 40])
        s = 16 / sqrt(8log(2))
        open(runs, "w") do io
            println(io, "dn,ker,rb,snr,seed,lambda,k,A,B,T")
            for k in ks
                println(io, "20,4,0,Inf,1,$(2pi / k),$k,0,0,$(exp(-k^2 * s^2 / 2))")
            end
        end
        RES.fit_tf(runs, out)
        ls = readlines(out)
        @test split(ls[1], ",")[5] == "fwhm_px"
        r = split(ls[2], ",")
        @test parse(Float64, r[5]) ≈ parse(Float64, r[8]) * 2log(2) / pi rtol = 1e-3
        @test RES.read_tf_fit(out)[(20, 4, 0)] ≈ 16 rtol = 0.06
        old = joinpath(dir, "old.csv")
        write(old, "dn,ker,rb,snr,fwhm_fit_px,fwhm_formula_px,lambda_half_px,T_lowk,npoints\n20,4,0,Inf,16.5,16.3,37.9,0.94,5\n")
        @test_throws ArgumentError RES.read_tf_fit(old)
    end
end
