# Angle calibration of the gradient stencils: synthetic validation

Status: validation evidence (2026-11-17), independently reviewed; corrections applied. Answer: no VGT angle bias above
1 deg is shown; the default stencil adds a grid-locked term of up to +0.85 deg,
which `isotropic_checked` removes; the extra -0.9 deg term of k^-11/3 fields is
an artefact of the generator's low-k modes (gone when they are cut). Synthetic data only. Scripts:
`scripts/validate_angle_bias.jl` (runs), `scripts/summarize_angle_bias.py`
(harmonic fits); unit tests in `test/angle_calibration.jl`.

## Question

Does VGT.jl's own computation rotate the angles it returns, by how much, and
does a better stencil remove it? Foreground templates fitted to dust
polarization need mean directions accurate to about a degree, so any
systematic of that size would matter.

## What can and cannot bias the angle

- The kernel (stencil, intensity mask, window fit, Q/U smoothing) is
  equivariant under the square grid's mirrors: a transposed cube returns
  90 deg minus the angle to 1e-8 deg, and a cube flipped along axis 1 returns
  minus the angle. The transpose test is exact with the default weights, so
  VGT.jl has no uniform angle offset. The flip test is exact only with uniform
  Stokes weights (the default weight is the window's own pixel, half a pixel
  off the even window's centre, which shifts position, not angle), so the
  90-deg period rests on it with that caveat. Any bias is therefore locked to
  the grid, odd and 90-deg periodic in the field angle t (counter-clockwise
  from +axis 1): b(t) = A4 sin 4t + A8 sin 8t + ... Symmetry breakers on the
  input side (an elliptical beam, a reprojection, striped noise) are outside
  VGT.jl and can add terms this note does not cover.
- The axial von Mises fitter (doubled-angle mean of unit vectors) is rotation
  equivariant; isotropic noise dilutes the window coherence but does not rotate it.
- The stencil is the grid-locked term. On a plane wave of wavenumber k
  (rad/pixel) with gradient angle theta, central differences return
  theta + (k^2/24) sin 4theta + O(k^4) (towards the diagonals), Sobel
  theta - (k^2/48) sin 4theta, and central differences with [1, 4, 1]/6
  smoothing across them (`isotropic_checked`) cancel the k^2 term. At
  theta = 22.5 deg: central +0.61 / +1.60 / +2.58 deg at k = 0.5 / 0.8 / 1.0,
  Sobel -0.30 / -0.78 / -1.24, isotropic +0.005 / +0.035 / +0.089. The general
  problem of rotation-invariant discrete gradients is Farid & Simoncelli
  (2004, IEEE TIP 13, 496).

## Setup

`scripts/validate_angle_bias.jl`; every run is VGT as in production
(`VGT_QU_mv_2d_profiled`, `dn = 20`, `min_valid = dn^2 / 4`, `ker = 4`, 3-sigma
mask) and the bias is half the angle of the sum of exp(2i (psi - truth)) over
the evaluated pixels. Declared before running: a bias counts as material if
|A4| > 1 deg; then `isotropic_checked` would be the correction if it brings
|A4| below 0.3 deg at unchanged alignment.

- Gaussian random fields (exact truth, no interpolation): 512^2 maps, 8
  independent channels, power k^-beta with anisotropy exp(-zeta cos^2(phi_k - t))
  (structures along t), 3 px FWHM beam, level 60 sigma (the mask never acts),
  beam-correlated noise; S/N r = rms of the signal minus its 20 px boxcar mean
  over the noise rms. Main grid: beta 8/3 (thin channel of Kolmogorov
  turbulence), zeta 3, r = Inf, 4, 2, 1; side runs at r = 2: beta 11/3, zeta 1
  and 8, dn 10 and 40. t = 0, 5, ..., 85 deg; the same random phases at every t.
  Harmonics fitted per seed, mean and standard error over seeds.
- Declared before running (2026-11-17, step added after the k^-11/3 result below, before the run): an independent
  review predicts that the about -0.9 deg term of the k^-11/3 fields is a property of the generator, not of any
  stencil or of the window. On a 512^2 lattice the few lowest modes have quantised directions, and the generator's
  own gradient-power-weighted orientation gives A4 -1.16 and A8 -0.24 deg for k^-11/3 and -0.05 for k^-8/3 (a
  cut at |k| < 2 pi/32 drops it to -0.05; 1024^2 fields still give -0.87). Test: the k^-11/3 case (r = 2, zeta 3,
  dn 20) with the modes |k| < 2 pi/32 zeroed (`--lowk-cut 32`), the same 17 seeds and settings, all stencils.
  Expected: isotropic A4 near 0 and central A4 near +0.3 to +0.85 deg (the paired stencil term, +0.30). If the
  term stays near -0.9 deg, the window or the estimator is implicated and the prediction fails.
- CATS MHD (Burkhart et al. 2020), run b1p1, line of sight 3, 218 channels,
  thermal width 1.0, as the earlier calibration runs: the beam-smoothed cube is
  rotated by t (4x FFT upsampling, Keys cubic sampling; identity at 0 deg 8e-9,
  +45/-45 deg round trip 1.2e-3 of the peak), then beam-correlated noise (a new
  realisation per angle) at peak S/N Inf, 46 and 8. Evaluated on the disk of
  radius 96 px about the rotation centre, so the same gas is measured at every
  angle; c0 is then this snapshot's own VGT-to-field offset (it rotates with the
  sky) and only A4, A8 are grid bias. Fits by least squares over the angles,
  errors from the residuals.
- `scripts/summarize_angle_bias.py` gives the fits and the paired differences
  between stencils on the same cubes.

## Results (2026-11-17)

Random fields, deg (mean ± s.e. over seeds; A4 is the amplitude of
sin 4t, positive = towards the diagonals):

| case | central A4 | Sobel A4 | isotropic A4 | central - isotropic A4 (paired) | AM | fit rbar2 |
|---|---|---|---|---|---|---|
| no noise (5 seeds) | +0.834 ± 0.095 | -0.403 ± 0.096 | +0.006 ± 0.095 | +0.829 ± 0.002 | 0.993 | 0.355 |
| r = 4 (5 seeds) | +0.837 ± 0.120 | -0.397 ± 0.119 | +0.011 ± 0.119 | +0.827 ± 0.001 | 0.992 | 0.307 |
| r = 2 (5 seeds) | +0.791 ± 0.126 | -0.356 ± 0.126 | +0.024 ± 0.126 | +0.767 ± 0.003 | 0.990 | 0.225 |
| r = 1 (5 seeds) | +0.589 ± 0.101 | -0.164 ± 0.093 | +0.091 ± 0.096 | +0.498 ± 0.011 | 0.970 | 0.130 |
| r = 2, k^-11/3 (17 seeds) | -0.756 ± 0.235 | -1.197 ± 0.239 | -1.052 ± 0.238 | +0.295 ± 0.006 | 0.956 | 0.293 |
| r = 2, weak anisotropy (zeta 1) (5 seeds) | +0.704 ± 0.247 | -0.297 ± 0.246 | +0.046 ± 0.249 | +0.658 ± 0.014 | 0.920 | 0.139 |
| r = 2, strong anisotropy (zeta 8) (17 seeds) | +0.825 ± 0.030 | -0.441 ± 0.030 | -0.022 ± 0.030 | +0.847 ± 0.002 | 0.997 | 0.321 |
| r = 2, dn 10 (5 seeds) | +0.786 ± 0.106 | -0.319 ± 0.105 | +0.048 ± 0.106 | +0.739 ± 0.005 | 0.978 | 0.283 |
| r = 2, dn 40 (6 seeds) | +0.828 ± 0.107 | -0.424 ± 0.106 | -0.010 ± 0.107 | +0.838 ± 0.004 | 0.997 | 0.209 |

Seeds 1-4 and 8 (dn 40 also 7; k^-11/3 and zeta 8 also 11-22); the other seeds
were lost to a full disk or stopped. c0 is zero within 1.5 s.e. everywhere. The
cos 4t term, which the symmetry forbids, is -0.04 to -0.13 (up to 3.5 s.e.) in
the 5-seed rows and the same for all three stencils (paired 0.00 to 0.02); with
17 seeds it is zero (-0.02 ± 0.03, +0.04 ± 0.12), so the 5-seed errors are
lower bounds. A8 is at most 0.16 deg.

With k^-11/3 every stencil has an extra grid-locked term of about -0.9 deg
(central -0.76 ± 0.24, isotropic -1.05 ± 0.24, 17 seeds) on top of the paired
stencil difference (+0.30), towards the axes. Cause: the random-field generator.
On a 512^2 lattice the few lowest modes have quantised directions, and a steep
spectrum puts most of its gradient power in them. The check declared above was
run with modes |k| < 2 pi/32 zeroed (`--lowk-cut 32`, same 17 seeds, jobs on
the cluster; seeds 1-4, 8, 11-22):

| k^-11/3, r = 2 | central A4 | Sobel A4 | isotropic A4 | central - isotropic (paired) |
|---|---|---|---|---|
| all modes (17 seeds) | -0.756 ± 0.235 | -1.197 ± 0.239 | -1.052 ± 0.238 | +0.295 ± 0.006 |
| modes below 2 pi/32 zeroed (17 seeds) | +0.469 ± 0.069 | -0.327 ± 0.070 | -0.061 ± 0.070 | +0.530 ± 0.004 |

The prediction held: isotropic A4 is -0.06 ± 0.07 (zero) and central +0.47 ± 0.07,
inside the expected +0.3 to +0.85; c0 and B4 are zero within errors, A8 below
0.01. The term is not a property of VGT, the window or the stencil, so it is
not a bias to calibrate. The k^-11/3 rows without the cut describe the
generator, not VGT; the paired stencil term (now +0.53) is what applies.

CATS b1p1, deg:

| S/N | central A4 | Sobel A4 | isotropic A4 | central - isotropic A4 | c0 (snapshot offset) | AM |
|---|---|---|---|---|---|---|
| no noise | +0.185 ± 0.003 | +0.001 ± 0.003 | +0.061 ± 0.003 | +0.124 ± 0.001 | -3.89 | 0.836 |
| 46 | +0.006 ± 0.110 | +0.022 ± 0.101 | +0.014 ± 0.106 | -0.007 ± 0.018 | -2.33 | 0.866 |
| 8 | -1.74 ± 2.81 | -1.48 ± 2.81 | -1.68 ± 2.85 | -0.06 ± 0.36 | -3.10 | 0.318 |

Reading:

- Central differences have a grid-locked stencil term towards the diagonals
  of +0.30 to +0.85 deg (paired against `isotropic_checked` on the same cubes)
  and +0.12 deg on the noiseless CATS cube. The absolute A4 of the default is
  +0.59 to +0.89 deg on the k^-8/3 random fields, -0.76 ± 0.24 deg on the
  k^-11/3 fields (generator artefact, +0.47 ± 0.07 with the low-k modes cut) and +0.19 deg on CATS (diluted: the
  disk's field angles are fixed and have <cos 4 phi> of their own, so the
  noiseless cube also fits B4 = -0.082 ± 0.003 and the amplitude is 0.20): below 1 deg
  in every case with power. CATS at S/N 8 (AM 0.32) has no power. The stencil
  term follows the wavenumbers that carry the gradients: steeper spectra and
  CATS have less beam-scale power, and noise dilutes it (r = 1: paired +0.50).
- `isotropic_checked` removes the stencil term: A4 = -0.02 to +0.09 deg on the
  k^-8/3 fields and +0.06 on CATS, at the same alignment (AM equal to
  0.001-0.003), and -0.06 ± 0.07 on the k^-11/3 fields with the generator's
  low-k modes cut.
  `sobel_checked` overshoots to the other side (-0.16 to -0.47 deg), as the
  plane-wave expansion predicts (half the central term, opposite sign).
- The window size does not matter (dn 10, 20, 40: paired +0.74, +0.77, +0.84).
- On CATS the noise realisation moves the mean angle of a 190 px disk by about
  0.4 deg at S/N 46 and by 10-20 deg at S/N 8 (AM 0.32), and the snapshot's own
  offset c0 is -2.3 to -3.9 deg. These are not calibratable: they are the
  sample variance of the field estimate, not a property of the stencil.
- Decision: no bias above 1 deg is shown, so the default stays
  `central_checked`. `isotropic_checked` is an opt-in comparison mode for
  analyses where the stencil term matters (a field running near 22.5 deg mod 45
  to the pixel grid over the whole map, where it does not average out). The
  -0.9 deg term of the uncut k^-11/3 fields comes from the generator, so it does
  not weigh against the new stencil.

Not done: other CATS runs (sub-Alfvenic, subsonic), other beams than 3 px, a
direct test on observed data with an independent field tracer (that belongs in
the experiment repositories).
