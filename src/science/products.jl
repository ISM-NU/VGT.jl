# Final science products from accumulated VGT Q/U/I maps.
#
# The VGT kernel accumulates gradient orientations. This file smooths those
# Stokes-like maps, rotates gradients into magnetic-field orientation proxies,
# and computes reliability/coherence diagnostics.

function imfilter_gaussian(d::AbstractMatrix, p::Number)
    # Gaussian smoothing pools nearby orientation evidence before converting
    # Q/U into the final angle map.
    return imfilter(d, Kernel.gaussian(p), NA())
end

function compute_stokes_intensity(cube)
    # Collapse the velocity axis to make the intensity normalization map.
    return sum(cube; dims=3)[:, :, 1]
end

function compute_reliability(q_smooth::AbstractMatrix, u_smooth::AbstractMatrix, intensity_smooth::AbstractMatrix)
    # Reliability is the fractional polarized/coherent amplitude of the
    # smoothed orientation signal: sqrt(Q^2+U^2)/I.
    reliability = Array{Float64}(undef, size(q_smooth))
    @inbounds for i in eachindex(reliability)
        intensity_value = intensity_smooth[i]
        reliability[i] = (isfinite(intensity_value) && intensity_value > 0.0) ?
            sqrt(q_smooth[i]^2 + u_smooth[i]^2) / intensity_value :
            NaN
    end
    return reliability
end

@inline function _invalid_finalize_intensity(value)
    # The mask rejects non-finite and non-positive intensity, because
    # reliability is normalized by intensity.
    return !isfinite(value) || value <= 0
end

function finalize_vgt_products(q_accum::AbstractMatrix, u_accum::AbstractMatrix, intensity_accum::AbstractMatrix,
                               ker::Int,
                               zero_stokes_policy::Symbol=:nan_before_smoothing,
                               output_support_mode::Symbol=:smoothed_stokes_support)
    # Convert raw accumulated Q/U/I maps into the final psi and reliability
    # products. The policy arguments should not be changed casually, because they alter output support and masks.
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))
    q_before_smoothing = copy(q_accum)
    u_before_smoothing = copy(u_accum)
    intensity_before_smoothing = copy(intensity_accum)
    if zero_stokes_policy === :nan_before_smoothing
        @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing, q_accum, u_accum)
            q_accum[idx] == 0 && (q_before_smoothing[idx] = NaN)
            u_accum[idx] == 0 && (u_before_smoothing[idx] = NaN)
        end
    end
    @inbounds for idx in eachindex(q_before_smoothing, u_before_smoothing,
                                   intensity_before_smoothing, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_before_smoothing[idx] = NaN
            u_before_smoothing[idx] = NaN
            intensity_before_smoothing[idx] = NaN
        end
    end
    q_smooth = imfilter_gaussian(q_before_smoothing, ker)
    u_smooth = imfilter_gaussian(u_before_smoothing, ker)
    intensity_smooth = imfilter_gaussian(intensity_before_smoothing, ker)
    @inbounds for idx in eachindex(q_smooth, u_smooth, intensity_accum)
        if _invalid_finalize_intensity(intensity_accum[idx])
            q_smooth[idx] = NaN
            u_smooth[idx] = NaN
        end
    end
    if output_support_mode === :direct_fit_support
        @inbounds for idx in eachindex(q_smooth, u_smooth, q_accum, u_accum)
            if q_accum[idx] == 0.0 && u_accum[idx] == 0.0
                q_smooth[idx] = NaN
                u_smooth[idx] = NaN
            end
        end
    end
    # atan(U,Q)/2 recovers the axial gradient orientation. VGT uses gradients as
    # a magnetic-field proxy after a 90-degree rotation, hence +pi/2.
    psi = 0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)
    reliability = compute_reliability(q_smooth, u_smooth, intensity_smooth)
    return (; q_smooth, u_smooth, intensity_smooth, psi, reliability)
end

function headless_angle_difference(a, b)
    # Difference between two axial angles, wrapped so 0 and pi represent the
    # same orientation.
    return mod.(a .- b .+ pi / 2, pi) .- pi / 2
end

"""
    block_mean(a, b; valid=isfinite.(a)) -> Matrix{Float64}

Mean of `a` over non-overlapping `b` x `b` blocks, using only entries where
`valid` is true; a block with no valid entry is NaN. Partial blocks at the
right and top edges are dropped, so output pixel (i, j) covers input pixels
(i-1)b+1:ib by (j-1)b+1:jb and the binned WCS stays exact.
"""
function block_mean(a::AbstractMatrix, b::Int; valid::AbstractMatrix{Bool}=isfinite.(a))
    b >= 1 || throw(ArgumentError("block size must be positive; got $b"))
    nx, ny = size(a) .÷ b
    out = Matrix{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        total = 0.0
        count = 0
        for jj in (j - 1) * b + 1:j * b, ii in (i - 1) * b + 1:i * b
            valid[ii, jj] || continue
            total += a[ii, jj]
            count += 1
        end
        out[i, j] = count == 0 ? NaN : total / count
    end
    return out
end

"""
    bin_header(header, b) -> header

Copy of a cube or map header for a grid of `b` x `b` blocks of axes 1 and 2:
`NAXIS1/2`, `CRPIX1/2`, `CDELT1/2` and `CDi_j` (j = 1, 2) are rescaled so each
binned pixel centre has the world coordinates of its block centre. `PC` and
the reference value are unchanged. Returns `nothing` for `nothing`.
"""
bin_header(::Nothing, ::Int) = nothing
function bin_header(header, b::Int)
    b == 1 && return header
    out = deepcopy(header)
    for axis in 1:2
        n = "NAXIS$axis"
        haskey(out, n) && (out[n] = Int(out[n]) ÷ b)
        # Input pixel x lies at binned coordinate (x - 0.5) / b + 0.5.
        c = "CRPIX$axis"
        haskey(out, c) && (out[c] = (Float64(out[c]) - 0.5) / b + 0.5)
        d = "CDELT$axis"
        haskey(out, d) && (out[d] = Float64(out[d]) * b)
        for i in 1:2
            cd = "CD$(i)_$(axis)"
            haskey(out, cd) && (out[cd] = Float64(out[cd]) * b)
        end
    end
    return out
end

"""
    bin_products(finalized, fit_rbar2, fit_support, header, b)

Products on the `b` x `b` binned grid: Q, U and I are averaged over the pixels
where all three are valid, then psi and reliability are recomputed from the
averages (never by averaging angles). The fit diagnostics are plain block
means. Returns `(; psi, reliability, fit_rbar2, fit_support, header)`.
"""
function bin_products(finalized, fit_rbar2::AbstractMatrix, fit_support::AbstractMatrix, header, b::Int)
    valid = isfinite.(finalized.q_smooth) .& isfinite.(finalized.u_smooth) .&
            isfinite.(finalized.intensity_smooth)
    q = block_mean(finalized.q_smooth, b; valid)
    u = block_mean(finalized.u_smooth, b; valid)
    intensity = block_mean(finalized.intensity_smooth, b; valid)
    psi = 0.5 .* atan.(u, q) .+ (pi / 2)
    return (; psi, reliability=compute_reliability(q, u, intensity),
            fit_rbar2=block_mean(fit_rbar2, b), fit_support=block_mean(fit_support, b),
            header=bin_header(header, b))
end

# ---------------------------------------------------------------------------
# Noise-debiased coherence (opt-in helpers; no production output uses them).
#
# A coherence R = |<exp(2i psi)>| over a beam is a resultant length of unit
# phasors, so, like a polarization fraction, it is biased by noise (Montier et
# al. 2015a, A&A 574, A135, for the polarization analogue). Write each VGT
# angle as the noise-free angle plus an independent error with attenuation
# g = E[exp(2i delta)], and let R_t be the coherence of the noise-free angles in
# the same beam (what the dust in that beam sees). Then, with f^2 = E[R^2] of
# pure noise (1 / N_eff independent angles per beam):
#   E[R^2]                    = g^2 R_t^2 + (1 - g^2) f^2
#   E[Re(z_a conj z_b)]       = g_a g_b R_t^2   (two independent noise realizations)
#   E[<cos 2(psi_a - psi_b)>] = g_a g_b         (per-pixel angle agreement)
# so split_coherence's ratio of the last two estimates R_t^2 without a noise
# model. Plain floor subtraction (debias_coherence) assumes g = 0 for the noise
# term and g = 1 for the signal; it is exact only for i.i.d. phasors and
# otherwise a lower bound on R_t (it also lowers a noise-free R).
# Synthetic validation: scripts/validate_coherence_debias.jl.
# ---------------------------------------------------------------------------

function _beam_mean(re::AbstractMatrix, im::AbstractMatrix, ok::AbstractMatrix{Bool}, fwhm_px::Real)
    fwhm_px > 0 || throw(ArgumentError("fwhm_px must be positive; got $fwhm_px"))
    k = Kernel.gaussian(fwhm_px / (2 * sqrt(2 * log(2))))
    w = Float64.(ok)
    rs, is_, ws = (imfilter(a, k, Fill(0.0)) for a in (re .* w, im .* w, w))
    full = imfilter(ones(size(re)), k, Fill(0.0))
    z = Matrix{ComplexF64}(undef, size(re))
    @inbounds for i in eachindex(z)
        z[i] = (ok[i] && ws[i] >= 0.5 * full[i]) ? complex(rs[i], is_[i]) / ws[i] : complex(NaN, NaN)
    end
    return z
end

"""
    coherence_phasor(psi, fwhm_px) -> Matrix{ComplexF64}

Mean doubled-angle phasor `<exp(2i psi)>` of an angle map (radians) over a
Gaussian beam of FWHM `fwhm_px` pixels, using only finite angles (each beam is
normalised by the kernel weight of its finite pixels). Pixels whose beam holds
less than half of the full kernel weight, or whose own angle is NaN, are NaN.
`abs` of the result is the beam coherence R.
"""
function coherence_phasor(psi::AbstractMatrix, fwhm_px::Real)
    ok = isfinite.(psi)
    p0 = ifelse.(ok, psi, 0.0)
    return _beam_mean(cos.(2 .* p0), sin.(2 .* p0), ok, fwhm_px)
end

"""
    coherence_floor2(psi_noise, fwhm_px; mask=trues(size(psi_noise))) -> Float64

Noise floor `f^2 = E[R^2]` of the beam coherence on a VGT angle map of pure
noise (same `dn`, `ker`, stencil and beam as the map it is compared with; e.g.
emission-free channels plus a constant offset so the intensity mask never
acts), averaged over finite pixels inside `mask`. `1 / f^2` is the effective
number of independent angles per beam; VGT angles are correlated over the
window, so it is far below the pixel count.
"""
function coherence_floor2(psi_noise::AbstractMatrix, fwhm_px::Real; mask::AbstractMatrix{Bool}=trues(size(psi_noise)))
    size(mask) == size(psi_noise) || throw(ArgumentError("mask size must match psi_noise"))
    r2 = abs2.(coherence_phasor(psi_noise, fwhm_px))
    v = [r2[i] for i in eachindex(r2) if mask[i] && isfinite(r2[i])]
    isempty(v) && throw(ArgumentError("no finite pure-noise coherence inside mask"))
    return sum(v) / length(v)
end

"""
    debias_coherence(R, floor2) -> same shape as R

Floor-subtracted coherence `sqrt(max(R^2 - f^2, 0) / (1 - f^2))`. Unbiased for
the expected-phasor coherence of i.i.d. angles; for VGT maps it is a lower
bound on the noise-free beam coherence (see the comment block above), kept as
a cheap diagnostic. `floor2` is a scalar or an array like `R`, `0 <= f^2 < 1`.
NaN stays NaN.
"""
function debias_coherence(R::AbstractArray, floor2)
    all(f -> isnan(f) || 0 <= f < 1, floor2) || throw(ArgumentError("floor2 must lie in [0, 1)"))
    return @. ifelse(isnan(R), NaN, sqrt(max(R^2 - floor2, 0.0) / (1 - floor2)))
end
debias_coherence(R::Real, floor2::Real) = debias_coherence([R], floor2)[1]

"""
    cross_coherence(z_a, z_b) -> real array

`sqrt(max(Re(z_a conj(z_b)), 0))` for beam phasors of two maps of the same sky
with independent noise. Free of the noise floor, but attenuated: it estimates
`sqrt(g_a g_b) R_t`, so it is low where the angles are noisy.
"""
function cross_coherence(z_a::AbstractArray{<:Complex}, z_b::AbstractArray{<:Complex})
    size(z_a) == size(z_b) || throw(ArgumentError("z_a and z_b must have the same size"))
    return @. sqrt(max(real(z_a * conj(z_b)), 0.0))
end

"""
    split_coherence(psi_a, psi_b, fwhm_px; gmin=0.05) -> (; R, g2)

Noise-free beam coherence from two VGT angle maps of the same sky with
independent noise (two noise realizations, or two halves of the channels):
`R = sqrt(Re(z_a conj z_b) / g2)` with `z` the beam phasors and
`g2 = <cos 2(psi_a - psi_b)>` over the same beam, the product of the two maps'
angle attenuations. Needs no noise model. NaN where `g2 < gmin` (angles too
noisy to correct); can exceed 1 by noise, not clipped. Correlated noise between
the two maps (e.g. spectral correlation between neighbouring channels) biases
both the numerator and `g2` up and partly cancels.
"""
function split_coherence(psi_a::AbstractMatrix, psi_b::AbstractMatrix, fwhm_px::Real; gmin::Real=0.05)
    size(psi_a) == size(psi_b) || throw(ArgumentError("psi_a and psi_b must have the same size"))
    ok = isfinite.(psi_a) .& isfinite.(psi_b)
    a0, b0 = ifelse.(ok, psi_a, 0.0), ifelse.(ok, psi_b, 0.0)
    za = _beam_mean(cos.(2 .* a0), sin.(2 .* a0), ok, fwhm_px)
    zb = _beam_mean(cos.(2 .* b0), sin.(2 .* b0), ok, fwhm_px)
    g2 = real.(_beam_mean(cos.(2 .* (a0 .- b0)), zeros(size(a0)), ok, fwhm_px))
    num = real.(za .* conj.(zb))
    # a branch, not ifelse: ifelse evaluates sqrt also where g2 < 0
    R = map((n, g) -> (isfinite(g) && g >= gmin) ? sqrt(max(n, 0.0) / g) : NaN, num, g2)
    return (; R, g2)
end
