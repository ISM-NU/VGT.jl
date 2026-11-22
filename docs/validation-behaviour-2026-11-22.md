# VGT behaviour with noise, smoothing, NaNs, edges and seams: pipeline map and synthetic measurements

Status: validation evidence (synthetic data only), 2026-11-22; independently reviewed, corrections applied. Script:
`scripts/validate_behaviour.jl` (metrics and factors declared in its header before the
runs); Slurm driver `hpc/irgetas/behaviour.slurm`.

## 1. Pipeline map

Each stage with what it assumes about noise, NaNs and resolution, and what happens when
the assumption fails. File references are to this branch (it carries the opt-in
`edge_mask` of branch `edge-mask`). Literature claims are marked FT (full text read) or
abs (abstract only).

| Stage | Code | Assumes | When the assumption fails |
|---|---|---|---|
| Input smoothing and binning (outside the package) | `scripts/smooth_cube_to_beam.jl:127-128`: NaN-aware normalised Gaussian convolution, kept where the kernel weight on finite pixels is > 0.5; then block means | smoothing only lowers noise; normalisation removes the blank | near a blank the smoothed values average fewer pixels; measured (section 2, survey edges): masking before smoothing turns directions across the edge by about 0.03 in T at S/N 8 and not measurably at S/N 64, and a stricter cut does not remove it; mechanism not isolated |
| Noise threshold | `src/science/vgt.jl:344-367` (central: all 5 stencil pixels >= threshold), `:368-381` (3 x 3 stencils), `:334-343` (backward: centre only); `src/science/noise.jl:19-31` robust sigma from negative pixels | one noise level per cube (scalar mode); an intensity cut makes gradients signal dominated | a 1 sigma cut passes about a sixth of pure-noise pixels (more as clusters when the noise is correlated); where the noise rises (mosaic edges) a scalar cut passes noise; the cut acts on intensity, not on gradient significance, so bright smooth gas with weak gradients is not removed |
| min_valid | `vgt.jl:426-435` (counts finite gradients in the window) | samples are many and roughly independent | with noise correlated over a grain of FWHM F the gradient field is correlated over about the same scale, so 100 samples in a 20 x 20 window hold only about 400 / (2 pi sigma_c^2), 20 to 40, independent ones for F = 3-4 px |
| Gradient stencil | `vgt.jl:141-168` | noise is white or isotropic | for Gaussian-correlated noise the central-difference noise is `sigma sqrt((1 - rho(2)) / 2)`, rho(2) the noise correlation at 2 px: 0.71 sigma (white), 0.48 sigma (F = 3 px), 0.27 sigma (F = 6 px); an elliptical grain gives the gradient noise a preferred direction (`validation-gradient-whitening.md`); grid-locked bias (`validation-angle-calibration.md`) |
| Whitening (opt-in) | `vgt.jl:267-274`, `:357-361`, `:451-455` | noise gradient covariance known from emission-free channels, one signal direction per window | back-mapping is exact only for one direction per window; covariance from the same channels over-fits (earlier note) |
| Window fit (axial von Mises) | `src/science/fitters.jl:118-146`; window `vgt.jl:419-423` | each finite gradient has unit weight whatever its size; the window mean is one direction | weak (noise) gradients count as much as strong ones, so low-S/N windows are diluted or pulled to the noise direction; the window `i-dn/2+1 : i+dn/2` is centred half a pixel off the pixel centre in both axes |
| Stokes accumulation over channels | `vgt.jl:456-472` | weight = intensity of the window-centre pixel in that channel; positive and finite | a negative centre value enters Q/U rotated by 90 deg (negative weight); a NaN centre value makes that pixel's Q/U NaN for good (any channel); the weight is one noisy pixel, not the window |
| Intensity normaliser | `vgt.jl:683-684` (chunk), `products.jl:57-73` | spectrum finite in every channel | one NaN channel makes the channel sum NaN and the pixel is blanked in every product, whatever `edge_mask` says |
| Output Q/U smoothing (ker) | `src/science/products.jl:51-67` (`imfilter(..., Kernel.gaussian(ker), NA())`: ker is the Gaussian sigma in px, FWHM 2.355 ker) | Q/U are smooth on the ker scale | with `nan_before_smoothing` pixels without any direct fit are filled from neighbours (normalised convolution); `smoothed_stokes_support` keeps them, `direct_fit_support` drops them (`:85-92`) |
| Reliability | `products.jl:18-29`, `:96` | I and Q/U cover the same channels | I sums every channel while Q/U only the channels with a fit, so the coherence also falls where channels lack fits |
| Edges | array border: windows skipped within dn/2 (`vgt.jl:413`); `edge_mask` (`vgt.jl:197-236`, `:416`, `:501-508`; `products.jl:74-84`) | blanks are NaN | zero-filled blanks are data (a hard step); `edge_mask` reacts only to non-finite pixels |
| Chunks | `src/config.jl:361-369` padding dn/2 + 1; smoothing after stitching (`src/pipeline/run.jl:93-97`) | | chunked and one-piece runs agree (`test/chunk_seams.jl`, `test/edge_mask.jl`) |

### What the literature says per stage

- Input smoothing: the Lazarian-group recipe smooths before the gradients, sigma = 2 px
  (Yuen & Lazarian 2017, ApJL 837, L24; Lazarian & Yuen 2018, ApJ 853, 96; FT) or 4 px on
  GALFA (Lazarian et al. 2018, ApJ 865, 46; FT). With white noise the alignment stays good
  up to noise/signal about 1 and is near 0 at about 3 (Lazarian & Yuen 2018, FT); for
  synchrotron gradients a wider filter holds up better at high noise (Lazarian et al. 2017,
  ApJ 842, 30; FT); a 3 px kernel was best at 10-30 % noise (Hu & Lazarian 2024, ApJ 972, 17;
  FT). Every test uses white Gaussian noise.
- Threshold: "mask out the raw gradient whose corresponding intensity is less than three
  times of RMS noise" (Hu 2022, MNRAS 511, 829; Hu et al. 2023, MNRAS 524, 2994; FT). The cut
  is on intensity, not on gradient significance, and its level is not tested.
- Block size: a Gaussian orientation histogram needs blocks of about 1/8 to 1/4 of the map
  side (Yuen & Lazarian 2017, FT); errors converge near 100 px for H I (Lazarian & Yuen 2018,
  FT); 16 x 16 px is called sufficient and 15-60 px typical (Hu & Lazarian 2023, FT); the LMC
  foreground paper used 20 x 20 px of 16' Parkes data (Alina et al. 2023, MNRAS 518, 4466;
  FT). No paper gives a minimum number of unmasked gradients per block.
- Channels: a channel is thin if its width is below the velocity dispersion on the scale
  studied; thermal broadening wider than the channel washes out the velocity term (Lazarian
  & Yuen 2018; Hu et al. 2023; FT). Channels are combined as pseudo-Stokes weighted by the
  channel intensity (Hu & Lazarian 2023, FT), as VGT.jl does; no paper weighs channels by
  noise.
- Output smoothing: Q/U of the gradients smoothed to about 10' FWHM (Hu 2022; Hu & Lazarian
  2024; FT); no error propagation through it.
- Edges and NaNs: failed blocks are NaN (Lazarian et al. 2018, FT); mosaic edges named as a
  sensitivity issue, not tested (Hu & Lazarian 2024, FT).
- Uncertainty: the fit standard error tracks the angle error (Lazarian & Yuen 2018, FT). No
  paper found uses split-data reproducibility.
- Silent: correlated or anisotropic (interferometric, mosaic) noise, the noise bias of the
  fitted angle at a given gradient S/N, stencil choice under noise, per-block NaN fraction,
  noise-weighted channel combination.

## 2. Synthetic measurements

### Setup

- Fields: CATS MHD snapshots (Burkhart et al. 2020), optically thin PPV with thermal width
  1.0 (code units), 64 channels, 3 px beam. Sub-Alfvenic views: b1p1 and b1p.32, each seen
  along map axes 3 and 2 (mean field in the plane). Super-Alfvenic views: b.1p.032 along
  axes 3 and 2 (strong angle structure; VGT reaches AM 0.45 there even without noise, so
  these rows show whether a change generalises, not absolute accuracy).
- Chain and settings as in the script header: smooth to 2 beams with NaN-aware
  normalised convolution, bin 2 x 2, dn 20, ker 4, min_valid 100, 1 sigma threshold.
  Without noise the default reaches AM 0.874 on the sub-Alfvenic views.
- Noise: unit-rms white noise convolved with a Gaussian grain of FWHM 3 px (the beam) unless
  stated, scaled so that S/N = median peak brightness / noise rms on the native grid. The
  gradient S/N on the analysis grid (median signal gradient over the rms of one noise
  gradient component) is 0.15, 0.25, 0.44, 0.77, 2.4 at S/N 2, 4, 8, 16, 64.
- Statistics: mean over views x 3 noise seeds (12 units sub-Alfvenic, 6 super-Alfvenic);
  "±" is the standard error over those units. Differences are paired (same view, seed and
  pixels). The four sub-Alfvenic views come from two MHD runs, and at S/N 16 and above the
  seeds nearly duplicate, so this "±" understates the spread over fields; where a result
  rests on few views the per-view values are given. The truth smoothed further to 16 px (T_prod) gives the same ranking everywhere and
  is listed only where it differs.
- Not covered: coverage-related policies. CATS boxes have emission in every pixel, so
  coverage stays 1 in every non-edge experiment and the support and zero-Stokes policies
  are identical here (section 3).

### Results table

AM is against the true field; dAM is paired against the row's reference.

| Factor | Setting | S/N 4 | S/N 8 | S/N 16 | S/N 64 | Reading |
|---|---|---|---|---|---|---|
| Noise type (ref: grain 3 px) | grain 3 px, AM | 0.587 ± 0.033 | 0.830 ± 0.013 | 0.876 ± 0.012 | 0.877 ± 0.012 | |
| | white, AM (same per-pixel rms) | 0.877 | 0.881 | 0.878 | 0.875 | white noise is mostly removed by the 2-beam smoothing; AM follows the gradient S/N (white at S/N 2: gradient S/N 0.40, AM 0.826; grain at S/N 8: 0.44, 0.830) |
| | grain 6 px, dAM | -0.210 ± 0.015 | -0.094 ± 0.014 | -0.014 ± 0.002 | +0.000 | coarser grain is worse |
| | + channel lag-1 0.25, dAM | -0.010 ± 0.005 | -0.005 ± 0.002 | -0.001 | -0.000 | small |
| | + channel lag-1 0.5, dAM | -0.042 ± 0.006 | -0.013 ± 0.005 | -0.002 | -0.001 | |
| | anisotropic grain (axis ratio 1.05, A = 0.047 native), bias (deg) | 9.7 ± 1.3 | 4.0 ± 0.7 | 1.7 ± 0.5 | 0.3 ± 0.6 | pull toward the noise direction (grain 3 px isotropic: 1.3, 1.4, 0.9, 0.3) |
| Threshold (ref: 1 sigma) | none / 0.5 / 2 / 3 sigma, dAM | +0.000 / +0.001 / -0.011 ± 0.003 / -0.030 ± 0.006 | | -0.000 / 0.000 / -0.000 / -0.000 | | 1 sigma or lower is best in filled fields; 3 sigma costs at low S/N. The cut's role of rejecting emission-free pixels is not tested (filled boxes) |
| min_valid (ref: 100) | 50 / 140 / 200, dAM | -0.000 / -0.000 / -0.002 ± 0.001 | | 0 | | no effect in filled maps |
| Stokes weight, sub-Alfvenic (ref: centre intensity) | clipped at 0, dAM | +0.000 | +0.000 | +0.000 | -0.000 | negative centre weights do not matter |
| | squared, dAM | -0.017 ± 0.004 | -0.016 ± 0.001 | -0.015 ± 0.003 | -0.021 ± 0.003 | |
| | square root, dAM | -0.000 ± 0.002 | +0.005 ± 0.002 | +0.008 ± 0.002 | +0.014 ± 0.003 | at S/N 64: b1p.32 views +0.017, +0.030; b1p1 views +0.005, +0.006 |
| | uniform, dAM | -0.014 ± 0.005 | +0.002 ± 0.004 | +0.013 ± 0.005 | +0.025 ± 0.006 | at S/N 64: b1p.32 +0.035, +0.052; b1p1 +0.004, +0.007. The S/N 4 loss does not replicate in `wsplit` (+0.009 ± 0.007) |
| Stokes weight, super-Alfvenic | square root / uniform, dAM | -0.002 / -0.006 | -0.002 / -0.007 ± 0.001 | -0.002 / -0.008 ± 0.001 | -0.002 / -0.015 ± 0.004 | the gain does not generalise |
| Grid and smoothing, sub-Alfvenic (ref: chain2 dn 20 ker 4, anisotropic grain) | chain2, AM | 0.549 ± 0.039 | | 0.873 ± 0.012 | 0.877 ± 0.012 | |
| | chain2 + whitening, AM | 0.561 (bias 14.5 -> 9.1 deg) | | 0.874 | 0.877 | |
| | native grid dn 40 ker 9 (Q/U smoothing only), dAM | -0.474 ± 0.025 (bias 43 deg) | | -0.252 ± 0.039 (bias 23 deg) | +0.041 ± 0.005 | pulled by the noise direction at low S/N |
| | native grid + whitening, dAM | -0.394 ± 0.020 | | -0.034 ± 0.013 | +0.045 ± 0.006 | per view: S/N 64 +0.025 to +0.069 (all 4 positive, view-level ± 0.011); S/N 16 b1p1 -0.076, -0.073, b1p.32 +0.003, +0.010 |
| | smoothed, not binned, dn 40 ker 8, dAM | -0.048 ± 0.006 | | +0.006 ± 0.001 | +0.008 ± 0.001 | binning costs little |
| Same, super-Alfvenic | native + whitening, dAM | -0.077 ± 0.014 | | +0.019 ± 0.012 | +0.043 ± 0.013 | two views: S/N 64 +0.072 and +0.015 |
| Window size (ref: dn 20) | dn 12 / 32, AM | 0.457 / 0.735 | | 0.803 / 0.944 | | larger windows win on these near-uniform fields; see section 3 |

Whitening control: with C = identity the native-grid result equals the unwhitened one in
every unit.

### Split-data reproducibility

R12 = AM between the two half-cube maps. Halves: odd/even channels, alternating 2-channel
blocks, every 4th channel (channels 1, 5, 9, ... against 3, 7, 11, ...).

| Channel lag-1 | Split | R12 at S/N 4 | half-map AM vs truth at S/N 4 | R12 at S/N 16 |
|---|---|---|---|---|
| 0 | odd/even | 0.382 ± 0.052 | 0.515 | 0.968 |
| 0.25 | odd/even | 0.427 ± 0.050 | 0.523 | 0.973 |
| 0.5 | odd/even | 0.557 ± 0.038 | 0.526 | 0.982 |
| 0 | 2-channel blocks | 0.377 ± 0.054 | 0.517 | 0.967 |
| 0.25 | 2-channel blocks | 0.385 ± 0.054 | 0.514 | 0.967 |
| 0.5 | 2-channel blocks | 0.429 ± 0.047 | 0.498 | 0.969 |
| 0.25 | every 4th | 0.240 ± 0.051 | 0.421 | 0.938 |

- Channel-correlated noise inflates odd/even reproducibility (+0.045 at lag-1 0.25, +0.175
  at 0.5, S/N 4) while the half maps are barely more accurate (+0.008 and +0.011).
  Alternating 2-channel blocks inflate by +0.008 ± 0.006 at lag-1 0.25 (consistent with 0)
  and +0.052 ± 0.010 at 0.5.
- R12 tends to 1 at high S/N (0.996 at S/N 64) while AM against the truth stays at 0.877: it
  measures the noise term only, not the method's own offset from the field.
- R12 does not rank the Stokes weightings by accuracy (`wsplit`): at S/N 64 uniform weights
  are the most accurate (AM 0.903 against 0.877) and the least reproducible (gap R12 0.976
  against 0.990). Reproducibility cannot select between weightings; smoothing choices were
  not tested this way.

### Survey edges

Half-plane footprint, edge across the mean field (n1) or along it (n2), distance in
analysis pixels. T = <cos 2(psi - edge normal)>, so T > 0 means directions across the
edge. Because the truth crosses the edge in n1 and runs along it in n2, a pure loss of
accuracy changes T with opposite signs in the two, while a turn toward one edge direction
changes it with the same sign. The table therefore gives the directional part
dir = (dT_n1 + dT_n2) / 2 and the accuracy part acc = (dAM_n1 + dAM_n2) / 2, both paired
against the same pixels without the edge.

| Case | S/N | 0-5 px dir / acc | 5-10 px dir / acc | 10-15 px dir / acc | 20-30 px acc |
|---|---|---|---|---|---|
| NaN applied after smoothing (VGT only) | 8 | -0.007 / -0.076 | +0.006 / -0.036 | +0.004 / -0.011 | 0.000 |
| | 64 | -0.024 / -0.067 | -0.008 / -0.034 | -0.002 / -0.009 | 0.000 |
| NaN before smoothing, cut 0.5 (production chain) | 8 | +0.035 / -0.082 | +0.032 / -0.040 | +0.016 / -0.014 | 0.000 |
| | 64 | -0.003 / -0.070 | +0.005 / -0.035 | +0.005 / -0.010 | 0.000 |
| cut 0.99 | 8 | +0.017 / -0.082 (coverage 0.60) | +0.024 / -0.052 | +0.017 / -0.020 | 0.000 |
| NaN cases with `edge_mask` | 8, 64 | masked | masked | +0.008 / -0.004 (coverage 0.8) | 0.000 |
| zero outside | 8 | -0.976 / -0.784 | -0.845 / -0.684 | -0.395 / -0.270 | 0.000 |
| | 64 | -0.672 / -0.546 | -0.377 / -0.282 | -0.116 / -0.058 | 0.000 |
| noise x 13 rising toward the edge (e-fold 25 native px), scalar threshold | 8 | -0.131 / -0.824 | -0.107 / -0.736 | -0.134 / -0.682 | -0.401 (30-50 px: -0.124) |
| | 64 | +0.010 / -0.076 | +0.006 / -0.033 | -0.006 / +0.002 | +0.001 |
| same with `edge_mask` | 8 | masked | masked | -0.160 / -0.679 | -0.401 |

Standard errors: dir ±0.002 to ±0.02, acc ±0.001 to ±0.04 (see the CSVs).

- Near a NaN edge VGT loses 0.03 to 0.08 in AM within 10 px, without a preferred
  direction. The loss does not shrink from S/N 8 to 64, so it is not a noise effect of
  partial windows alone; the mechanism was not isolated.
- Masking before smoothing adds a turn across the edge of about +0.03 in T within 10 px at
  S/N 8 and none measurable at S/N 64. It persists at cut 0.99 (+0.017 to +0.024), so a
  stricter smoothing cut is not a fix.
- `edge_mask` removes the affected band; what remains at 10-15 px is within ±0.01.
- Zero-filled outside is a hard step in every channel: directions run along the edge within
  10-15 px. `edge_mask` does not see zeros. Blank the outside with NaN.
- Noise rising toward the edge costs accuracy wherever the local S/N drops below about 8,
  with directions nearly random by a local S/N of about 1.5, and leans directions along the
  edge (dir -0.09 to -0.13 at S/N 8). A per-pixel threshold changes nothing (not shown in
  the table: within ±0.08 of the scalar threshold), because a 1 sigma cut passes the same
  fraction of noise at any noise level. `edge_mask` does not reach this loss beyond 10 px.

### Mosaic seams

One half of the map has a block of channels scaled by a gain g, as a step (L = 0) or a
linear ramp over L = 175 native px. Wing = the blue channels holding 18 % of the flux; core
= a quarter of the channels around the peak, 81 % of the flux. T_seam = <cos 2(psi - seam
tangent)>; dir and acc as for the edges (dir > 0 = turned along the seam), pooled over a
seam across (n1) and along (n2) the mean field.

| Channels, fade | S/N | 0-5 px dir / acc | 5-10 px dir / acc | 20-30 px dir / acc | 30-50 px dir / acc |
|---|---|---|---|---|---|
| wing, to 0, step | 8 | +0.050 / -0.027 | +0.031 / -0.027 | -0.001 / -0.012 | +0.004 / -0.007 |
| | 64 | +0.041 / -0.016 | +0.023 / -0.018 | +0.003 / -0.005 | +0.001 / -0.004 |
| wing, to 0, ramp 175 px | 8 | +0.007 / -0.015 | +0.003 / -0.013 | +0.003 / -0.017 | +0.008 / -0.010 |
| | 64 | +0.029 / -0.010 | +0.020 / -0.012 | +0.014 / -0.012 | +0.007 / -0.009 |
| wing, to 0.5, ramp 175 px | 8 | +0.003 / -0.010 | -0.000 / -0.008 | +0.003 / -0.009 | +0.004 / -0.007 |
| | 64 | +0.015 / -0.006 | +0.008 / -0.005 | +0.007 / -0.005 | +0.003 / -0.004 |
| core, to 0, step | 64 | +0.482 / -0.406 | +0.166 / -0.153 | +0.009 / +0.025 | +0.011 / +0.004 |
| core, to 0, ramp 175 px | 64 | +0.425 / -0.369 | +0.432 / -0.371 | +0.370 / -0.296 | +0.128 / -0.081 |
| wing NaN in one half, NaN replaced by 0 after smoothing | 8 / 64 | +0.001 / +0.018 dir | +0.000 / +0.009 dir | | |
| wing or core NaN in one half, as is | 8, 64 | that half blanked (coverage 0.5) | | | |

Standard errors: ±0.001 to ±0.006 for the wing rows, up to ±0.04 for the core rows.

- A fade of channels that carry little flux turns directions along the seam by a few
  thousandths to 0.03 in T, spread over the ramp, and costs about 0.01 in AM. A fade of the
  line core dominates the gradients over the whole ramp.
- A NaN in some channels blanks every pixel of that region (section 1, intensity
  normaliser). Replacing those NaNs with 0 after smoothing restores coverage and near the
  seam is much milder than a fade to 0 (wing, S/N 8, 0-5 px: dir +0.001 against +0.050),
  because smoothing has not spread the step; beyond 20 px the two agree.
- The synthetic fades scale the signal only; a real one-block fade also changes the noise.

## 3. Reading and limits

- VGT accuracy follows the gradient S/N on the analysis grid: AM about 0.6 at 0.25, 0.83 at
  0.44, and the noise-free value above about 0.8. Beam-scale grain noise keeps much more
  gradient noise through the 2-beam smoothing than white noise of the same per-pixel rms
  (a factor of about 4 in peak S/N at 3 px per beam), which matters because published
  noise tests use white noise.
- At the production resolution the default settings are close to the best tested
  everywhere: threshold at 1 sigma, min_valid irrelevant in filled maps, channel correlation
  costs 0.01 at S/N 4 and nothing above 16.
- A finer grid with Q/U smoothing only (no input smoothing) reaches higher AM at S/N 64 in
  every view (+0.015 to +0.072), and is pulled hard toward the noise direction at lower S/N
  unless the gradients are whitened. Whitened, at S/N 16 it ties or wins in the b1p.32 and
  super-Alfvenic views and loses 0.07 in both b1p1 views; at S/N 4 it loses in all. The crossover lies between S/N 16
  and 64 and was not located. Whitening on the smoothed grid removes about a third of the
  remaining low-S/N bias (with the noise axis at 45 deg to the field, the worst case).
- Square-root and uniform channel weights gain at high S/N in one of the two sub-Alfvenic
  runs (b1p.32: +0.017 to +0.030 and +0.035 to +0.052 at S/N 64), under +0.01 in the
  other, and lose a little in the super-Alfvenic views; the result does not yet
  generalise.
  Square-root weights can be tried without code changes: `weight_file` = sqrt(max(I, 0)) with
  `noise` = sqrt(threshold) gives exactly the square-root weights and the default mask
  (the reliability map then uses the weight cube's channel sum).
- Window size: dn 32 beats dn 20 here, but the sub-Alfvenic truth is nearly uniform, so
  larger windows win by construction; this does not settle the resolution trade-off.
- Limits: two sub-Alfvenic MHD runs and one super-Alfvenic run, one thermal width, one channel count, Gaussian grain
  noise with one anisotropy; filled maps, so coverage policies (support, zero-Stokes,
  min_valid at low fill) are not tested; edges are straight; seams are one straight line.
