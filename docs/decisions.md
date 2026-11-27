# Decisions

Status: active reference. Dated decisions that change the meaning of a
production output. Newest first. Each entry says what changed, why, and what
it does to results made before it.

## 2026-11-27: a NaN voxel removes only its own channel

Decided by Denis Andreyev. Code: `src/science/vgt.jl` (kernel),
`src/science/products.jl` (`compute_stokes_intensity`). Test:
`test/nan_centre.jl`.

Old behaviour. A window fit centred on a non-finite voxel got the Stokes weight
`max(NaN, 0) = NaN`, so Q/U/W became NaN for the whole pixel, and the
channel-summed sky intensity was NaN there as well. One NaN channel removed the
pixel from psi and reliability, even though the other channels had valid fits,
so a per-channel NaN flag acted as a pixel mask.

New behaviour. The kernel skips the fit of a channel whose window centre is
non-finite, and the sky intensity sums the finite voxels where the plain sum is
non-finite. Pixels blank in every channel sum to 0 and stay masked as before;
pixels without a non-finite voxel keep the plain sum bit for bit. A NaN voxel
now gives the same Q/U/W, psi and reliability as a voxel set below the noise
threshold (weight 0); the test checks this.

Effect on results. Output is bitwise identical wherever each pixel is finite in
all channels of the run or non-finite in all of them (blank sky as NaN or as 0).
Only inputs with NaN in some channels of a pixel change: those pixels now get a
direction from their other channels.

## 2026-11-25: `psi_reliability` is normalised by the fitted Stokes weight

Decided by Denis Andreyev. Code: `src/science/vgt.jl` (kernel),
`src/science/products.jl` (`finalize_vgt_products`, `bin_products`).

Old behaviour. The map was `P = sqrt(Qs^2 + Us^2) / Is`. Qs and Us were Q/U
smoothed with the NaN-aware Gaussian after setting pixels without a fit to NaN,
so they were averages over fitted pixels only. Is was the intensity summed over
every channel, fitted or not, smoothed over every positive-intensity pixel.
Numerator and denominator had different masks and different sums, so P was not
a fraction. An earlier audit of synthetic scenes found P > 1 in 25-70% of
pixels, with maxima of 4-27; the scenes measured for this change gave 0.4%
(max 1.03), 15% (max 3.6) and 17% (max 2.1).
The Stokes weight was the raw window-centre intensity, so a negative centre
pixel (noise near zero emission) entered Q/U with a negative weight, which is a
vote for the orientation rotated by 90 degrees.

New definition. The kernel also accumulates `W`, the sum of the same Stokes
weights that go into Q and U (fitted windows only, including the `rbar2_weight`
factor when that option is on). `W` is smoothed with the same filter and the
same NaN mask as Q and U, and `psi_reliability = sqrt(Qs^2 + Us^2) / Ws`.
With non-negative weights this is in [0, 1] up to rounding (triangle
inequality), also after `output_bin` block averaging. The Stokes weight is now
`max(centre intensity, 0)`. The channel-summed intensity is still used, only to
mask sky. The clamp can drop a pixel from the psi support if every fitted
window there had a negative centre (its Q/U become exactly zero); with the
default `smoothed_stokes_support` the halo usually refills it, and on the
negative-background scene below the psi NaN mask did not change.
`direct_fit_support` was not checked.

Effect on results.

- psi is bitwise identical when no fitted window has a negative centre pixel
  (checked on a 96x88x16 synthetic scene, one-piece, chunked and binned runs).
  Where centre pixels were negative, psi changes slightly: on the same scene
  shifted so that 26% of pixels were negative, 44% of psi pixels moved, median
  0.001 deg, maximum 0.6 deg.
- Reliability changes scale everywhere, not only where it was above 1. In the
  positive-background scene the median went from 0.43 to 0.77 and the share of
  pixels above 1 from 0.4% to 0, because channels without a fit no longer
  dilute the denominator. In the negative-background scene the share above 1
  went from 15% (max 3.6) to 0.
- Old reliability maps, and any threshold or high/low split tuned on them, are
  not comparable with new ones. Rerun. New maps carry the header card
  `RELNORM = 'FITWEIGHT'`; old maps do not have it.
- Chunk files now store `w`. Chunk files from older runs cannot be stitched by
  the new code (the pipeline always rewrites its chunks, so this matters only
  when stitching old intermediates by hand).

Tests: `test/vgt_edge_cases.jl`, testsets "psi_reliability stays in [0, 1]"
and "negative centre intensity gives zero Stokes weight".

Two follow-up fixes in the same change.

- Zero-Stokes masking. Under `:nan_before_smoothing`, Q and U were each set to
  NaN where that map alone was exactly zero. A fitted angle of exactly 0 or
  90 deg gives U == 0 or Q == 0, so a perfectly axis-aligned field had an
  all-NaN psi (a pure x-ramp: 289 fitted pixels, 0 finite psi). Now a pixel is
  masked only when Q and U are both zero, which is what a missing fit looks
  like. On noisy data a single exact zero essentially never occurs, so psi is
  unchanged there; the existing parity tests pass unchanged.
- `psi_support` counts only fits that contributed weight to Q/U. A fit whose
  weight was clamped to zero (negative centre pixel) no longer counts.

Tests: "axis-aligned field keeps psi (U == 0 is not a missing fit)" and the
updated "negative centre intensity gives zero Stokes weight".
