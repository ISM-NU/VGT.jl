# Faint-gas VGT on MHD cubes: gradient-amplitude weighting and channel pre-filter

Status: validation evidence (2026-11-17). Synthetic (CATS MHD) data only. Neither variant passed the declared rule;
the default (`gradient_weight = "unit"`, no pre-filter) stays. `gradient_weight = "amplitude2"` is kept as an opt-in
option for a declared follow-up (see "Open"), not for production runs. Scripts: `scripts/validate_faint_gas.jl`,
`scripts/figures/faint_gas_summary.py`.

## Question

Where a foreground layer's gradients are near the noise (gradient S/N `rho_g` about 1), can VGT's direction be made
more accurate by

- V1, `gradient_weight = "amplitude2"`: each gradient's doubled angle weighted by `|g|^2` in the window fit (the
  structure tensor). Under isotropic noise a gradient's doubled-angle variance is about `2 s^2 / |G|^2`, so this is
  the inverse-variance weight, and the noise adds nothing to the expected sums, so the angle stays unbiased.
- V2, a Gaussian pre-filter of every channel (FWHM 4 px on a 3 px beam) before the gradients, with the noise mask and
  the Stokes weights still on the unfiltered channel (`input_file` = filtered cube, `weight_file` = raw cube). This is
  the Lazarian group's noise recipe (Lazarian et al. 2017, ApJ 842, 30; Lazarian et al. 2018, ApJ 865, 46: Gaussian
  smoothing of sigma 2-4 px before the gradients).

Literature (full texts read): none of the Lazarian-group papers read weights gradients by amplitude inside a block
(Hu et al. 2019, Nature Astronomy 3, 776, Lu et al. 2020 and Hu & Lazarian 2023 were not read); Lazarian et al.
2018 (footnote 13) suggest filtering the smallest-amplitude gradients as noise but do not test it, and Lazarian &
Yuen 2018 (ApJ 853, 96) link the largest amplitudes to shocks. PCA filtering (Hu, Yuen & Lazarian 2018, MNRAS 480,
1333) was not tested here: it ranks velocity structure by variance, which favours the bright layer, and Hu et al.
2020 (ApJ 888, 96) kept all components in their foreground runs. Whitening (`gradient_whiten`) is not part of this
test: the injected noise is isotropic, where whitening is the identity (`validation-gradient-whitening.md`).

## Setup (declared in the script header before any result)

- Cubes: CATS (Burkhart et al. 2020) b.1p.32, b.1p1, b1p.32, b1p1 at t800, lines of sight 2 and 3: 8 cubes.
  Optically thin PPV, 48 channels over +-4 std(v_los), beam 3 px, beam-correlated noise, 3 noise seeds per cube, the
  same noise for every variant.
- VGT as a foreground run on a 2' cube with 40" pixels: 4-channel windows sliding by one channel over |v| <= 2.5
  std(v_los), dn 20, min_valid 100, central stencil, scalar 3 sigma mask, ker 4.
- Template as the two-layer foreground template: per window, column-weighted Stokes at FWHM 15 px (10' on a 40" grid)
  and the unit Stokes at 180 px (120'). Layers: "all" windows, and "wing" (window centre below -1 std(v_los), a
  fainter layer on one side of the line). Truth: the velocity-resolved, density-weighted Stokes of the sky-plane
  field over the layer's channels, smoothed the same way.
- Score: AM = <cos 2(psi - psi_true)>; paired difference to plain VGT per 32 px block on the pixels finite for every
  variant; error from a hierarchical bootstrap (cubes, then blocks within cubes; seeds averaged per block).
- Gradient S/N: per 20 px block, `rho_g^2 = <|g|^2> / <2 sigma_g^2> - 1` over masked-in pixels.
- Amendment before any variant score: at peak S/N 1-4 the 3 sigma mask leaves almost no pixel, so the grid became peak
  S/N 8, 12, 16, 23, 46 (wing `rho_g` 0.30, 0.78, 1.13, 1.62, 3.18).
- Pass rule: (1) wing layer, pooled over S/N 8-16 and at the S/N levels matched to the target's faint layer (here 16
  and 23): dAM at 10' > 0 at >= 2 sigma, and not worse at 120'; (2) both layers at S/N
  23 and 46: dAM >= -0.005; (3) the sign of (1) in at least 6 of 8 cubes.

## Results

Paired dAM at 10' (variant minus plain, +- bootstrap std; cubes with dAM > 0 in brackets):

| layer, S/N | plain AM | V1 amplitude2 | V2 pre-filter 4 px | V1 + V2 |
|---|---|---|---|---|
| wing, 8 | 0.564 | +0.017 +- 0.016 (5/8) | +0.005 +- 0.034 (4/8) | +0.031 +- 0.029 (5/8) |
| wing, 12 | 0.637 | +0.010 +- 0.012 (4/8) | -0.036 +- 0.014 (2/8) | -0.010 +- 0.012 (2/8) |
| wing, 16 | 0.641 | +0.008 +- 0.009 (3/8) | -0.036 +- 0.010 (0/8) | -0.012 +- 0.012 (1/8) |
| wing, 23 | 0.660 | -0.002 +- 0.008 (2/8) | -0.042 +- 0.009 (0/8) | -0.026 +- 0.012 (1/8) |
| wing, 46 | 0.679 | -0.012 +- 0.008 (3/8) | -0.042 +- 0.008 (0/8) | -0.039 +- 0.011 (0/8) |
| wing, pooled 8-16 | | +0.010 +- 0.006 | -0.029 +- 0.009 | -0.005 +- 0.009 |
| all, 16 | 0.665 | +0.005 +- 0.015 (4/8) | -0.073 +- 0.026 (1/8) | -0.055 +- 0.021 (1/8) |
| all, 23 | 0.657 | -0.006 +- 0.009 (4/8) | -0.058 +- 0.018 (0/8) | -0.060 +- 0.018 (0/8) |
| all, 46 | 0.703 | -0.012 +- 0.006 (2/8) | -0.051 +- 0.013 (0/8) | -0.055 +- 0.014 (0/8) |

Window-level AM at 10' against each window's own truth, by `rho_g` (all windows, cubes, seeds, S/N), and the paired
difference to plain (bootstrap over cubes; cubes with a positive difference):

| rho_g | n | plain | V1 | V2 | V1 + V2 |
|---|---|---|---|---|---|
| < 0.5 | 172 | 0.428 | 0.463 | 0.450 | 0.502 |
| 0.5-0.8 | 411 | 0.507 | 0.544 | 0.506 | 0.545 |
| 0.8-1.2 | 606 | 0.551 | 0.573 | 0.519 | 0.551 |
| 1.2-1.7 | 693 | 0.605 | 0.615 | 0.553 | 0.574 |
| 1.7-2.5 | 515 | 0.632 | 0.633 | 0.569 | 0.581 |
| 2.5-4 | 388 | 0.633 | 0.628 | 0.575 | 0.584 |
| > 4 | 309 | 0.640 | 0.631 | 0.565 | 0.570 |

| rho_g | V1 - plain | V2 - plain |
|---|---|---|
| < 0.5 | +0.035 +- 0.006 (8/8) | +0.022 +- 0.013 (5/8) |
| 0.5-0.8 | +0.037 +- 0.003 (8/8) | -0.001 +- 0.006 (4/8) |
| 0.8-1.2 | +0.022 +- 0.001 (8/8) | -0.032 +- 0.004 (0/8) |
| 1.2-1.7 | +0.010 +- 0.001 (8/8) | -0.052 +- 0.005 (0/8) |
| 1.7-2.5 | +0.002 +- 0.002 (6/8) | -0.063 +- 0.007 (0/8) |
| 2.5-4 | -0.005 +- 0.002 (1/8) | -0.058 +- 0.006 (0/8) |
| > 4 | -0.008 +- 0.004 (2/8) | -0.074 +- 0.007 (0/8) |

Rule: V1 fails (1) (+0.010 +- 0.006 pooled, 1.6 sigma; 3-5 of 8 cubes), (2) (-0.012 at S/N 46) and (3). V2 fails
all three. The layer test is scored only where every window of the layer is covered (159-310 of 512 blocks for the
wing layer at S/N 8-16), i.e. on brighter pixels than the `rho_g` label describes, which works against V1; V2 loses
at both the layer and the window level. The 120' scale could not be tested: the layer sum is NaN wherever one window has no coverage at 180 px,
which in these patchy cubes is almost everywhere below S/N 46 (a harness limitation; at S/N 46 the dAM at 120' is
-0.005 +- 0.041 for V1 and -0.086 +- 0.055 for V2).

## Reading

- Plain VGT's 10' agreement saturates near AM 0.63 above `rho_g` about 2; below that the noise costs 0.03-0.2.
  At `rho_g` 0.8-1.7 the loss is 0.03-0.08 (window level).
- V1 recovers part of that loss where `rho_g` < 1.7 (+0.01 to +0.04 per window, every cube) and costs 0.005-0.008
  above `rho_g` 2.5. On the declared layer test the two do not add up to a 2 sigma gain, so V1 is not a general
  improvement; its window-level gain at low `rho_g` is what the follow-up below would test.
- V2 lowers the ceiling (0.57 instead of 0.63): a 4 px pre-filter on a 3 px beam appears to remove the small-scale
  structure that VGT's gradients read, and only at `rho_g` < 0.5 does the noise reduction win. The maps of the four variants
  share the same pattern of direction errors, set by the gas and field structure, not by the noise.
- The CATS wing windows at matched `rho_g` pass the 3 sigma mask on 30-60 % of pixels; observed faint layers that
  are smooth and bright (high coverage, low `rho_g`) are a different regime that this test does not reproduce.

## Open

- The 120' scale: on a 256 px box each cube gives about one independent 120' sample, so a fixed harness adds only
  8 points.
- An S/N-switched weighting (V1 only in windows with `rho_g` below about 1.7) is suggested by the window-level table but
  was not declared; it needs its own declared test, with the 120' layer sum counting empty windows as zero.
- Noise anisotropy and non-Gaussian (scan-stripe) noise are not in these cubes.
