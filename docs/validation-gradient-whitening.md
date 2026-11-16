# Gradient whitening against anisotropic correlated noise: synthetic validation

Status: validation evidence. Synthetic data only. Script:
`scripts/validate_gradient_whitening.jl` (run one case at a time with
`VGT_GW_CASES=ramp|ramp_mis|isotropic|varying|structured`; each takes 1 to 5 minutes;
`VGT_GW_SPAN` and `VGT_GW_REF_NV` set the rotation and the covariance channels of `varying`).

## Problem

Interferometric noise is spatially correlated and slightly anisotropic, so noise
gradients have a preferred direction (anisotropy
`A = |sum g^2 e^(2i theta)| / sum g^2`, a few per cent for an elliptical noise grain; the synthetic noise below has 0.0475).
Where the gradients are not signal dominated, each window's axial mean is pulled
toward that direction. The intensity mask does not help: the intensity can be far
above the noise while the gradients are not.

## Option: `gradient_whiten`

Config (`[noise]`): `gradient_whiten = false` (default, output bitwise
unchanged) and `gradient_noise_cov = [sxx, syy, sxy]`, the covariance of one
gradient vector in emission-free channels, required and positive definite when on.

With `C_n = C / sqrt(det C)` every gradient is mapped by `W = C_n^(-1/2)` after
the existing validity checks, so noise gradients become isotropic and have no
preferred angle. The determinant is 1, so magnitudes stay comparable. Each
window's fitted angle is mapped back with `S = C_n^(1/2)`: the unit vector
at the fitted angle is multiplied by `S` and its angle is reported. That is exact
when the window has one signal direction and approximate otherwise.

`gradient_noise_cov_file` replaces the constant with a per-pixel map when the
anisotropy changes across the map (see "Spatially varying noise anisotropy");
set one of the two. `C` must come from the same stencil as `gradient_mode` (central or Sobel; the
backward-difference mode rejects the option). `measure_gradient_noise(cube;
gradient_mode)` and `scripts/analyze_inputs.jl cube.fits --noise-vrange=LO,HI`
print `gradient_noise_cov` ready to paste, plus A and its direction.

## Setup

- Noise: white Gaussian noise convolved with an elliptical Gaussian, FWHM 3 px
  (geometric mean), axis ratio 1.05, long axis 45 deg from the truth. Measured on
  an independent cube: sigma_gx = 0.486, sigma_gy = 0.475, A = 0.0475, preferred
  gradient direction -29.3 deg (the gradients prefer the short axis of the
  kernel), so `D = -49.3` deg is the full pull from the true gradient axis (20 deg)
  toward the noise direction.
- Signal on a baseline of 1000 (the intensity mask never acts). Ramp: plane with
  the true gradient axis, slope `G` scaled to `G/sigma_g`, `sigma_g` = rms of the two
  measured component deviations; channel envelope `1 + 0.4 sin(2 pi k / nv)`.
  Structured: the `build_many_filaments` and `build_curved_fibers` generators of
  `validate_mw_like_synthetic.jl` without their noise, scaled so the median
  central-difference gradient magnitude of the signal on the evaluated pixels is
  `S/N * sigma_g`.
- VGT as in production: 128x128x8 cubes, `dn = 20`, `min_valid = 100`,
  `central_checked`, then Q/U smoothing with `ker = 4` and direct-fit support.
  4 noise seeds. `±` is the standard error over seeds.
- Bias: doubled-angle mean offset from the truth; `bias/D` is the fraction of the full
  pull toward the noise direction. Scatter: rms headless difference from the mean
  direction. Coverage: finite fraction of interior pixels (ramp) or of evaluated
  pixels (structured). Structured error: headless difference from the planted
  tangent on pixels above the 55th percentile of the signal with truth
  coherence of at least 0.75; "toward-noise bias" is the mean signed error with
  the sign set by the side of the truth on which the noise axis lies.

## Ramp, anisotropic noise

| G/σ_g | whitening | bias (deg) | bias/D | scatter (deg) | coverage |
|---|---|---|---|---|---|
| 0.10 | off | -47.83 ± 4.15 | +0.97 | 30.32 | 1.000 |
| 0.10 | on | -20.93 ± 28.18 | +0.42 | 47.21 | 1.000 |
| 0.25 | off | -30.64 ± 2.01 | +0.62 | 28.37 | 1.000 |
| 0.25 | on | -1.59 ± 2.29 | +0.03 | 37.93 | 1.000 |
| 0.50 | off | -10.01 ± 0.56 | +0.20 | 9.66 | 1.000 |
| 0.50 | on | -0.71 ± 0.67 | +0.01 | 9.89 | 1.000 |
| 1.00 | off | -2.73 ± 0.22 | +0.06 | 3.39 | 1.000 |
| 1.00 | on | -0.22 ± 0.17 | +0.00 | 3.30 | 1.000 |
| 2.00 | off | -0.42 ± 0.07 | +0.01 | 0.79 | 1.000 |
| 2.00 | on | -0.17 ± 0.07 | +0.00 | 0.77 | 1.000 |
| 4.00 | off | +0.07 ± 0.02 | -0.00 | 0.25 | 1.000 |
| 4.00 | on | -0.04 ± 0.02 | +0.00 | 0.25 | 1.000 |

Whitening removes the pull: bias/D falls from 0.62 to 0.03 at G/sigma_g = 0.25 and
from 0.20 to 0.01 at 0.5. At G/sigma_g = 0.1 the gradients carry no signal and the
angle is random after whitening (bias ±28 deg over seeds, scatter 47 deg), where
without it the angle sits near the noise direction with lower scatter (30 deg).
The lower scatter below G/sigma_g = 0.5 is the pull itself, not accuracy.
At G/sigma_g >= 1 scatter is equal and at 4 the results agree to 0.1 deg. Coverage
is unchanged (whitening does not change validity).

## Covariance mis-estimated by 20% in A

| G/σ_g | whitening | bias (deg) | bias/D | scatter (deg) | coverage |
|---|---|---|---|---|---|
| 0.25 | on, A -20% | -10.52 ± 1.48 | +0.21 | 37.34 | 1.000 |
| 0.25 | on, A +20% | +6.81 ± 3.07 | -0.14 | 37.46 | 1.000 |
| 0.50 | on, A -20% | -2.60 ± 0.66 | +0.05 | 9.94 | 1.000 |
| 0.50 | on, A +20% | +1.14 ± 0.66 | -0.02 | 9.78 | 1.000 |
| 1.00 | on, A -20% | -0.72 ± 0.18 | +0.01 | 3.33 | 1.000 |
| 1.00 | on, A +20% | +0.27 ± 0.16 | -0.01 | 3.27 | 1.000 |
| 2.00 | on, A -20% | -0.22 ± 0.07 | +0.00 | 0.77 | 1.000 |
| 2.00 | on, A +20% | -0.12 ± 0.07 | +0.00 | 0.76 | 1.000 |

A 20% error in A leaves a residual bias of about 20% of the original pull where
the pull is large (0.25: -10.5 and +6.8 deg against -30.6 deg without whitening)
and under 1 deg from G/sigma_g = 1. Whitening with an estimate off by 20% is still
much better than not whitening.

## Isotropic noise: result must equal the default

| G/σ_g | whitening | max abs difference in psi vs off (deg) | bias (deg) | scatter (deg) |
|---|---|---|---|---|
| 1.00 | off | 0.00e+00 | -0.42 | 3.17 |
| 1.00 | on, C = identity | 2.54e-14 | -0.42 | 3.17 |
| 1.00 | on, measured C | 4.48e-01 | -0.19 | 3.17 |
| 4.00 | off | 0.00e+00 | -0.04 | 0.24 |
| 4.00 | on, C = identity | 1.27e-14 | -0.04 | 0.24 |
| 4.00 | on, measured C | 1.18e-02 | -0.05 | 0.24 |

With `C = identity` whitening reproduces the default to rounding (1e-14 deg). With
the C measured on isotropic correlated noise (A = 0.004, sampling noise only) the
maps differ by up to 0.45 deg at G/sigma_g = 1 and 0.01 deg at 4.

## Spatially varying noise anisotropy

Gradient noise whose anisotropy changes across the map is handled by
`gradient_noise_cov_file`: a (nx, ny, 3) FITS cube of (sxx, syy, sxy) used per
pixel. Gradients are mapped with `C^(-1/2)` of their own pixel, and the fitted
angle at the window centre is mapped back with `C^(1/2)` of that pixel.
`gradient_noise_cov_map` (and `analyze_inputs.jl --cov-map`) builds it from
emission-free channels: per-pixel `gx^2`, `gy^2`, `gx*gy` averaged over the
channels and smoothed with a NaN-aware Gaussian (default FWHM 3 dn).

Setup: four independent noise fields whose kernel long axis is spaced evenly
over the stated rotation, blended across the map with sqrt of linear hat weights
(covariance interpolated linearly, variance kept). Ramp signal as above. The
covariance (global and map) is measured on a separate cube of the stated number
of emission-free channels. Main run: 30 deg rotation, 48 channels, 4 seeds.

Local noise measured on the reference cube (left / middle / right third):

  third 1: A=0.0545 direction=-33.1 deg (D=-53.1 deg from the truth)
  third 2: A=0.0478 direction=-23.9 deg (D=-43.9 deg from the truth)
  third 3: A=0.0381 direction=-14.7 deg (D=-34.7 deg from the truth)
global: A=0.0451 direction=-24.9 deg

Residual preferred direction of whitened pure-noise gradients on an independent
cube (unit-weight resultant of the doubled angle; 0 = none):

| correction | left third | middle third | right third |
|---|---|---|---|
| off | 0.0217 | 0.0234 | 0.0178 |
| global C | 0.0045 | 0.0011 | 0.0053 |
| map FWHM 3 dn | 0.0044 | 0.0012 | 0.0068 |
| map FWHM 1.5 dn | 0.0052 | 0.0019 | 0.0076 |

Ramp bias per third and overall:

| G/σ_g | correction | bias left (deg) | bias middle | bias right | mean abs bias | scatter (deg) | coverage |
|---|---|---|---|---|---|---|---|
| 0.10 | off | -52.52 | -47.35 | -44.63 | 48.17 | 29.88 | 1.000 |
| 0.10 | global C | -12.67 | +30.04 | -7.21 | 16.64 | 47.56 | 1.000 |
| 0.10 | map FWHM 3 dn | -0.21 | +31.08 | -15.29 | 15.53 | 47.56 | 1.000 |
| 0.10 | map FWHM 1.5 dn | -1.76 | +31.61 | -16.93 | 16.77 | 47.92 | 1.000 |
| 0.25 | off | -33.15 | -27.86 | -27.79 | 29.60 | 27.55 | 1.000 |
| 0.25 | global C | -13.66 | +0.48 | -6.31 | 6.82 | 39.14 | 1.000 |
| 0.25 | map FWHM 3 dn | -11.06 | +3.16 | -13.95 | 9.39 | 39.11 | 1.000 |
| 0.25 | map FWHM 1.5 dn | -10.92 | +4.15 | -16.20 | 10.43 | 39.31 | 1.000 |
| 0.50 | off | -10.13 | -9.02 | -9.13 | 9.43 | 9.38 | 1.000 |
| 0.50 | global C | -0.20 | +0.16 | +0.21 | 0.19 | 9.83 | 1.000 |
| 0.50 | map FWHM 3 dn | +0.99 | +0.72 | -0.83 | 0.85 | 9.83 | 1.000 |
| 0.50 | map FWHM 1.5 dn | +1.20 | +1.10 | -1.29 | 1.20 | 10.00 | 1.000 |
| 1.00 | off | -2.48 | -2.38 | -1.91 | 2.26 | 3.09 | 1.000 |
| 1.00 | global C | -0.03 | +0.12 | +0.66 | 0.27 | 3.08 | 1.000 |
| 1.00 | map FWHM 3 dn | +0.26 | +0.25 | +0.39 | 0.30 | 3.07 | 1.000 |
| 1.00 | map FWHM 1.5 dn | +0.32 | +0.33 | +0.29 | 0.31 | 3.09 | 1.000 |
| 2.00 | off | -0.14 | -0.28 | -0.10 | 0.18 | 0.79 | 1.000 |
| 2.00 | global C | +0.10 | -0.04 | +0.14 | 0.09 | 0.78 | 1.000 |
| 2.00 | map FWHM 3 dn | +0.12 | -0.03 | +0.11 | 0.09 | 0.78 | 1.000 |
| 2.00 | map FWHM 1.5 dn | +0.12 | -0.02 | +0.10 | 0.08 | 0.78 | 1.000 |

Sensitivity to the size of the variation and to the number of channels used for the
covariance (same columns as above, `varying` case):

| axis rotation (deg) | covariance channels | residual resultant, mean of thirds: off / global C / map 3 dn | mean abs bias at G/σ_g = 0.5: off / global C / map 3 dn / map 1.5 dn | at G/σ_g = 1 |
|---|---|---|---|---|
| 30 | 48 | 0.0210 / 0.0036 / 0.0041 | 9.43 / 0.19 / 0.85 / 1.20 | 2.26 / 0.27 / 0.30 / 0.31 |
| 30 | 192 | 0.0210 / 0.0028 / 0.0032 | 9.43 / 0.47 / 0.50 / 0.65 | 2.26 / 0.23 / 0.29 / 0.32 |
| 60 | 48 | 0.0189 / 0.0062 / 0.0038 | 8.41 / 0.53 / 0.71 / 1.19 | 1.98 / 0.30 / 0.28 / 0.29 |
| 90 | 96 | 0.0166 / 0.0088 / 0.0031 | 6.83 / 1.14 / 0.36 / 0.29 | 1.55 / 0.46 / 0.26 / 0.22 |

Reading:

- The per-pixel map did **not** beat one global covariance for a 30 deg rotation.
  Both remove the pull (mean absolute bias 9.4 deg off; 0.19 deg with global C and
  0.85 / 1.20 deg with the map at 3 dn / 1.5 dn at G/sigma_g = 0.5, where one seed set
  has a standard error of roughly 0.5 deg), and the residual direction of the
  whitened noise is similar (0.0045 / 0.0011 / 0.0053 against 0.0044 / 0.0012 / 0.0068).
  At this size of variation (the measured local direction changes by about 20 deg
  between the left and right thirds, A 0.04 to 0.05) one global covariance
  leaves a residual pull of under 1 deg, about what the estimate of the map
  itself costs (finite channels, smoothing).
- The map wins when the variation is large: with a 90 deg rotation (local
  direction -47 to 0 deg between the thirds) the mean absolute bias at G/sigma_g = 0.5
  is 1.14 deg with global C and 0.36 (FWHM 3 dn) or 0.29 (1.5 dn) with the map, at 1 it is
  0.46 against 0.26 and 0.22, and the residual direction of pure-noise gradients
  is 0.0100, 0.0070, 0.0094 with global C against 0.0010, 0.0033, 0.0051 with the map.
  At 60 deg the residual resultant also favours the map (0.0045 / 0.0011 / 0.0059 against
  0.0090 / 0.0031 / 0.0066 for global C) but the angle bias at G/sigma_g >= 0.5 is
  about 1 deg or less for both.
- The two smoothing scales (3 dn and 1.5 dn) give the same results within
  noise on these 128 pixel maps; at 3 dn the smoothing width is half the map.
- At G/sigma_g = 0.1 and 0.25 the bias is dominated by the random angle of noise
  windows and no correction is distinguishable (differences of 10 deg between
  corrections are within seed scatter).
- On these synthetic maps the global covariance leaves no clearly visible local
  direction; the synthetic rotation is small compared with the smoothing on 128
  pixel maps. On maps where the measured local direction or A changes a lot, the
  per-pixel map is the safer choice; it costs nothing at high S/N (all
  corrections are within 0.1 deg of each other at G/sigma_g = 2).
- Windows without signal: the back-mapping with `S = C_n^(1/2)` turns uniformly
  random window angles into angles that favour the noise axis, with a resultant of
  about A/2. Whitened pure noise therefore keeps a small preferred direction along
  the original one, weaker than the unwhitened pull (compare the pure-noise rows
  of the varying-anisotropy tables).
- The covariance must come from channels independent of those it is tested on:
  a map built from the same channels removes their particular noise pattern, and
  tests on those channels look better than on new ones.

## Structured emission

Filament morphologies at a gradient S/N from 0.25 to 4 (`A` +/- 20% rows show the
mis-estimation sensitivity).

| scenario | S/N | whitening | median err (deg) | p90 err (deg) | toward-noise bias (deg) | coverage |
|---|---|---|---|---|---|---|
| many_filaments | 0.25 | off | 17.06 | 49.51 | +15.71 | 1.000 |
| many_filaments | 0.25 | on | 17.70 | 63.16 | +3.29 | 1.000 |
| many_filaments | 0.25 | on, A +20% | 19.53 | 66.49 | -0.79 | 1.000 |
| many_filaments | 0.25 | on, A -20% | 16.26 | 60.06 | +6.71 | 1.000 |
| many_filaments | 0.50 | off | 7.23 | 26.08 | +7.52 | 1.000 |
| many_filaments | 0.50 | on | 6.22 | 23.73 | +0.07 | 1.000 |
| many_filaments | 0.50 | on, A +20% | 6.41 | 24.78 | -1.79 | 1.000 |
| many_filaments | 0.50 | on, A -20% | 6.12 | 23.67 | +1.80 | 1.000 |
| many_filaments | 1.00 | off | 3.64 | 10.88 | +2.30 | 1.000 |
| many_filaments | 1.00 | on | 3.25 | 10.28 | -0.55 | 1.000 |
| many_filaments | 1.00 | on, A +20% | 3.31 | 10.54 | -1.19 | 1.000 |
| many_filaments | 1.00 | on, A -20% | 3.25 | 10.14 | +0.06 | 1.000 |
| many_filaments | 2.00 | off | 2.84 | 8.48 | +0.45 | 1.000 |
| many_filaments | 2.00 | on | 2.66 | 8.62 | -0.69 | 1.000 |
| many_filaments | 2.00 | on, A +20% | 2.65 | 8.78 | -0.99 | 1.000 |
| many_filaments | 2.00 | on, A -20% | 2.67 | 8.50 | -0.44 | 1.000 |
| many_filaments | 4.00 | off | 2.69 | 8.72 | -0.30 | 1.000 |
| many_filaments | 4.00 | on | 2.50 | 8.85 | -1.06 | 1.000 |
| many_filaments | 4.00 | on, A +20% | 2.50 | 8.92 | -1.29 | 1.000 |
| many_filaments | 4.00 | on, A -20% | 2.52 | 8.82 | -0.86 | 1.000 |

| scenario | S/N | whitening | median err (deg) | p90 err (deg) | toward-noise bias (deg) | coverage |
|---|---|---|---|---|---|---|
| curved_fibers | 0.25 | off | 6.41 | 32.22 | +0.56 | 1.000 |
| curved_fibers | 0.25 | on | 6.34 | 27.78 | -4.13 | 1.000 |
| curved_fibers | 0.25 | on, A +20% | 6.29 | 27.49 | -4.91 | 1.000 |
| curved_fibers | 0.25 | on, A -20% | 6.35 | 28.33 | -3.34 | 1.000 |
| curved_fibers | 0.50 | off | 5.46 | 25.18 | -0.59 | 1.000 |
| curved_fibers | 0.50 | on | 5.53 | 23.14 | -3.99 | 1.000 |
| curved_fibers | 0.50 | on, A +20% | 5.59 | 23.14 | -4.53 | 1.000 |
| curved_fibers | 0.50 | on, A -20% | 5.49 | 23.23 | -3.42 | 1.000 |
| curved_fibers | 1.00 | off | 5.09 | 22.60 | -1.03 | 1.000 |
| curved_fibers | 1.00 | on | 5.19 | 21.45 | -3.78 | 1.000 |
| curved_fibers | 1.00 | on, A +20% | 5.32 | 21.49 | -4.24 | 1.000 |
| curved_fibers | 1.00 | on, A -20% | 5.10 | 21.46 | -3.28 | 1.000 |
| curved_fibers | 2.00 | off | 4.92 | 21.36 | -1.31 | 1.000 |
| curved_fibers | 2.00 | on | 5.07 | 20.58 | -3.61 | 1.000 |
| curved_fibers | 2.00 | on, A +20% | 5.17 | 20.69 | -4.01 | 1.000 |
| curved_fibers | 2.00 | on, A -20% | 4.99 | 20.59 | -3.19 | 1.000 |
| curved_fibers | 4.00 | off | 4.87 | 20.58 | -1.50 | 1.000 |
| curved_fibers | 4.00 | on | 5.00 | 20.19 | -3.47 | 1.000 |
| curved_fibers | 4.00 | on, A +20% | 5.05 | 20.23 | -3.82 | 1.000 |
| curved_fibers | 4.00 | on, A -20% | 4.97 | 20.15 | -3.11 | 1.000 |

- many_filaments (field near 33 deg with 10 deg scatter): the pull toward the noise
  axis is +15.7 deg at S/N 0.25 and +7.5 at 0.5 without whitening and 3.3 and
  0.1 with it; median error falls from 7.2 to 6.2 deg at 0.5 and from 3.6 to 3.3 at 1.
  At 0.25 the p90 error is worse (63 vs 50 deg), as in the ramp: the angles are
  noise dominated and no longer pulled to one direction.
- curved_fibers (field direction varies over more than 100 deg): the pull has no net
  sign, and whitening does not help. It adds a small systematic: median error
  higher by 0.1 to 0.15 deg at every S/N, toward-noise bias about 2 deg more
  negative at S/N 4 (-1.5 to -3.5 deg), consistent with the back-mapping being exact
  only for a single direction per window. p90 is lower by 0.4 to 4 deg.
- Mis-estimating A by 20% changes these numbers by less than 1.5 deg in the
  toward-noise bias and less than 0.2 deg in median error.

## Tested and rejected: gradient S/N mask

An earlier option dropped gradients with `|g| < k sqrt(2) sigma_g` (pure-noise
survival `exp(-k^2)`). On the same ramp setup (no Q/U smoothing, 6 seeds; the option was removed after
this test) it did not reduce the
pull toward the noise direction:

| G/sigma_g | no mask bias / scatter (deg) | k = 1 | k = 2 |
|---|---|---|---|
| 0.5 | -9.6 / 12.8 | -10.1 / 11.5 | no fit survives min_valid |
| 1 | -2.6 / 4.3 | -3.1 / 3.8 | no fit survives min_valid |
| 2 | -0.4 / 1.1 | -0.5 / 1.1 | -0.7 / 1.3 |

The reason, measured on the noise alone: a radial cut on an anisotropic Gaussian keeps
proportionally more gradients along its major axis, so the unit-weight axial
anisotropy of the survivors rises (0.024 unmasked, 0.046 at k = 1, 0.135 at
k = 2). Removing weak gradients raises the noise anisotropy per gradient. The
option was removed; whitening addresses the covariance directly.

## Limits

One noise model (A = 0.0475, one kernel orientation, Gaussian), a plane ramp and
two filament morphologies, central stencil only, 4 seeds. In curved fields whitening
adds a small systematic and no benefit. The covariance has to be measured in
emission-free channels of the same cube. The per-pixel map was tested with smoothly
rotating anisotropy only, on 128 pixel maps; noise that differs between channels
(non-scalar noise models) is not tested. Testing the pipeline path found that the
chunk padding was one row short for the central stencil; that is fixed in b279d00
(`test/chunk_seams.jl`), and the whitened chunked run now equals the one-piece run.

## Literature (ADS abstracts only; full texts not read)

Established for VGT in the Lazarian-group papers found: averaging gradient
orientations over sub-blocks and fitting their distribution (Yuen & Lazarian
2017, ApJ 837, L24; Lazarian et al. 2018, ApJ 865, 46), spectral (PCA) filtering to
improve the tracing (Hu, Yuen & Lazarian 2018, MNRAS 480, 1333), and intensity
gradients as a companion tracer (Hu et al. 2019, ApJ 886, 17). In the abstracts I
found no gradient S/N mask and no treatment of spatially correlated, anisotropic
noise. Whitening with a measured noise covariance, and this test, are new to this
repository, not taken from those papers.
