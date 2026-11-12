# Synthetic Fitter Validation Notes

Status: active validation note. The Gaussian histogram baseline (`:lsqfit_reference`) and `separable_shrinkage` noise compared here were removed after this validation; the suite now evaluates `:axial_vonmises` only.

Last updated: 2026-06-11.

## Question

Can we justify `:axial_vonmises` as the production fitter over
`:lsqfit_reference` on synthetic cubes with known orientations?

## Current Result

The improved `scripts/validate_synthetic.jl` suite supports
`:axial_vonmises` as the better default for synthetic regimes that resemble
noisy, sparse, spatially structured spectral-line cubes.

Final run summary:

- Overall lower-error cases: `:axial_vonmises` 34/44,
  `:lsqfit_reference` 10/44.
- Multi-seed paired effect-size rows: `:axial_vonmises` 5/6,
  `:lsqfit_reference` 1/6.
- Sparse patchy emission: `:axial_vonmises` wins all tested fill fractions.
- Spatially varying finite-difference truth: `:axial_vonmises` wins all tested
  noise levels.

This is not a proof that `:axial_vonmises` is always lower error. In very clean,
high-SNR, uniform-gradient cases, `:lsqfit_reference` can be slightly better.
The measured high-SNR advantage is small compared with the broader sparse and
moderate-noise gains from the axial estimator.

## Adopted Accuracy Change

`scripts/validate_synthetic.jl` now includes a method-improvement search across
noisy uniform, weak-gradient, sparse-fill, mildly curved, and sharply curved
synthetic cubes.

The adopted config value is `gradient_mode = "central_checked"`:

- Uniform noisy: median error `4.066 deg -> 0.446 deg`, no coverage loss.
- Weak gradient: `5.061 deg -> 0.721 deg`, no coverage loss.
- Sparse fill: `0.569 deg -> 0.120 deg`, coverage ratio `0.95`.
- Mildly varying field: `0.535 deg -> 0.105 deg`, no coverage loss.
- Sharply varying field: `0.673 deg -> 0.204 deg`, no coverage loss.

This replaces the previous production default, `:nan_edge_backward`, with a
central finite-difference stencil that requires the center pixel and the four
derivative-neighbor pixels to be above threshold. The old backward path remains
available for historical comparisons.

Literature context: recent gradient-technique theory
(`arXiv:2406.01545`) emphasizes that measured gradient distributions are
affected by plane-of-sky and line-of-sight averaging, so candidate smoothing
operations must be treated as physical choices, not just denoisers. Hu, Yuen,
and Lazarian (2018, `arXiv:1802.08772`) show that PCA filtering can improve VGT
alignment and acts similarly to filtering along the velocity axis, but that is
kept out of production here because foreground work can contain multiple
velocity components that should not be blended without component validation.

## Accuracy Experiments Tried

The validation script now tests two non-production candidates:

- `rbar2_weight`: multiply each channel's Stokes contribution by the window
  doubled-angle coherence diagnostic.
- `rbar2_gate_0.50`: skip channel-window fits with `rbar2 < 0.50`.

The real-cube parameter validator can now compare these with:

```bash
julia scripts/validate_vgt_parameters.jl \
  --input path/to/cube.fits \
  --noise 0.16 \
  --fitter axial_vonmises \
  --min-valid auto \
  --rbar2-mode baseline,weight,gate0.50
```

`--min-valid auto` expands to fractional candidate thresholds and reports a
recommendation. It does not silently choose or apply a production value.

Both improved the synthetic scenarios tested in the final run:

- Uniform noisy: baseline 4.177 deg; `rbar2_weight` 3.975 deg;
  `rbar2_gate_0.50` 2.510 deg with some coverage loss.
- Sparse fill: baseline 0.690 deg; `rbar2_weight` 0.506 deg;
  `rbar2_gate_0.50` 0.482 deg.
- Varying angle: baseline 0.530 deg; `rbar2_weight` 0.390 deg;
  `rbar2_gate_0.50` 0.270 deg.

These are follow-up candidates only. They should not become production defaults
until checked on component-matched real-cube validation, foreground starlight
agreement, and Planck foreground-subtraction residual behavior.

Additional method candidates from the synthetic search:

- `:sobel_checked` also improved all five scenarios and was often slightly more
  accurate than `:central_checked`, but it uses a wider 3x3 stencil and is kept
  as a validation candidate because it can smooth sharper physical structure.
- Three-channel Hanning filtering improved all tested single-component
  synthetic cases, but it can mix velocity components, so it is not a production
  default.
- Larger `dn=14` windows improved some structured cases but lost coverage in
  noisy and weak-gradient cases; do not adopt larger windows globally.

## Edge-Detail Follow-Up

After adopting `:central_checked`, we ran an additional targeted check for the
main concern: whether the new gradient mode merely oversmooths sharper
structure.

The check compared raw accumulated gradient orientation before final Gaussian
product smoothing, using `:axial_vonmises` for all modes. It focused on
two-region boundaries and thin Gaussian filaments.

Boundary transition results:

- For 30, 45, and 60 degree two-region boundaries, all three modes had
  zero median interior error in the deterministic no-noise probe.
- The estimated transition width was not larger for `:central_checked`:
  `6, 6, 5` columns for central versus `7, 6, 6` for
  `:nan_edge_backward`.
- Boundary-zone errors remained ambiguous by construction because boundary
  windows contain two valid orientations.

Thin-filament results:

- Width 1.5 px: backward median error `13.64 deg`, central `1.21 deg`,
  Sobel `0.45 deg`.
- Width 2.5 px: backward `7.41 deg`, central `0.38 deg`, Sobel `0.19 deg`.
- Width 4.0 px: backward `2.37 deg`, central `0.19 deg`, Sobel `0.10 deg`.

Interpretation: this targeted synthetic probe does not support the idea that
`:central_checked` is only winning by washing out detail. It preserved boundary
transition width and substantially improved narrow-filament orientation
relative to the old backward stencil. `:sobel_checked` remains attractive in
idealized tests but still carries more smoothing risk because it uses the full
3x3 neighborhood.

## Gradient-Mode Ordering Check

Follow-up checks on 2026-06-11 tested whether the stronger statement
"Sobel is always better than central, and central is always better than
backward" is supportable.

Decision: do not make that universal claim.

Synthetic evidence:

- The main synthetic method search still supports `:central_checked` over the
  old backward path in the representative noisy, weak-gradient, sparse-fill,
  mildly curved, and sharply curved scenarios.
- `:sobel_checked` is not always lower-error than `:central_checked`: in the
  sharply varying synthetic field, central was slightly lower median error
  (`0.204 deg`) than Sobel (`0.207 deg`).
- In a deterministic clean linear-gradient cube, all three modes were exact
  within the tested precision. That alone falsifies any strict "always better"
  ordering.
- In an adversarial low-SNR linear-gradient probe with very weak per-pixel
  gradient and independent white noise, the old backward stencil beat central
  in all 12 seeds at noise levels `0.05` and `0.08`. This does not make
  backward the better production default, but it marks a real failure regime:
  central is not universally lower-error when the physical gradient is far
  below channel noise.

Real-cube internal comparison:

- No Planck 353 GHz Q/U maps or starlight-polarization catalogue files were
  present under local `data/`, so external accuracy validation could not be run
  in this workspace.
- On representative local H I/CO crops, central and Sobel produced similar maps
  in dense H I but diverged more in sparse or high-threshold regions. Median
  central-vs-Sobel final-psi residuals were about `0.65 deg` on the GASS crop,
  `0.93 deg` on the GASKAP crop, `8.16 deg` on the HNV signal crop, and
  `3.95 deg` on the CO crop.
- Sobel generally had slightly higher doubled-angle coherence (`rbar2`) but
  lower support or fewer fit attempts. On the sparse CO crop, central produced
  finite psi over `4.67%` of pixels with `1021` fit attempts, while Sobel
  produced `3.48%` finite psi with only `192` fit attempts.

Interpretation: Sobel remains a useful validation candidate, especially for
idealized smooth or filamentary structures, but it should not replace central
as a default without external starlight/Planck validation and an explicit
coverage policy. Central remains the conservative default because it improves
the representative synthetic suite relative to backward while using a narrower
stencil than Sobel.

## Milky-Way-Like Morphology Suite

`scripts/validate_mw_like_synthetic.jl` adds a dedicated synthetic foreground
morphology suite. It does not replace external Planck or starlight validation;
it tests whether the synthetic evidence covers more realistic Galactic
structures than smooth gradients and circular patches.

The script currently includes seven morphology scenarios:

- many elongated filaments;
- crossing filament families;
- curved fibers;
- two separated velocity components with different orientations;
- shells/arcs;
- a Planck-aligned foreground plus a misaligned contaminating component;
- sparse, broken CO-like filaments.

Each scenario reports `:nan_edge_backward`, `:central_checked`, and
`:sobel_checked` separately for clean high-truth-coherence pixels, all signal
pixels, and low-coherence ambiguous pixels. It now scans `min_valid` at
`72`, `50`, and `36` valid samples for the `dn=12` window, corresponding to
roughly 50%, 35%, and 25% of the local window. Ambiguous regions are not
treated as ordinary accuracy truth because crossings and component mixtures may
not have one physically correct orientation.

Representative clean-region median errors from the first run:

- Many filaments: backward `6.22 deg`, central `1.84 deg`, Sobel `2.04 deg`.
- Crossing filaments: backward `5.26 deg`, central `1.87 deg`,
  Sobel `1.48 deg`.
- Curved fibers: backward `5.57 deg`, central `1.17 deg`, Sobel `0.82 deg`.
- Velocity components: backward `6.83 deg`, central `0.75 deg`,
  Sobel `0.63 deg`.
- Shells/arcs: backward `3.94 deg`, central `1.45 deg`, Sobel `1.37 deg`.
- Planck-aligned plus misaligned component, scored against total mixture:
  backward `15.33 deg`, central `1.05 deg`, Sobel `0.92 deg`.
- Sparse CO-like filaments: backward `7.32 deg`, central `5.31 deg`,
  Sobel `11.62 deg`.

Coverage is as important as angle error. Sobel usually has the highest
fitter-level coherence but lower support than central, and this is severe in
the sparse CO-like case: clean-region support was `0.614` for backward,
`0.185` for central, and only `0.084` for Sobel. This supports keeping Sobel
as a validation candidate rather than making it the production default.

The Planck-aligned contaminant case also shows a real limitation: when a
misaligned component is present in other velocity channels, the integrated VGT
orientation tracks the total mixture well but does not reliably recover the
foreground-only orientation everywhere. That argues for component-matched
foreground validation before using any gradient-mode change for Planck
foreground subtraction.

Follow-up support scan result: lowering `min_valid` is the most promising
support-improvement candidate, but it should not be applied blindly. In the
sparse CO-like clean region with `:central_checked`, changing `min_valid` from
`72` to `50` increased support from `0.185` to `0.422` and improved the median
error from `5.31 deg` to `4.08 deg`; lowering further to `36` increased support
to `0.611`, but the 90th percentile error worsened from `32.83 deg` to
`49.62 deg`. Dense filament scenarios also gained support at `50`, with only
small median-error changes for central/Sobel, while backward differences tended
to add high-tail errors. The current recommendation is to validate
`min_valid ≈ 0.35 * dn^2` on real component-matched crops before changing the
production default.

Focused fitter comparison on Milky-Way-like regions: a full 128x128 grid over
all gradient modes and `min_valid` settings is impractical with
`:lsqfit_reference`, so the controlled comparison was run at 96x96 with
`:central_checked` and `min_valid = 50`, across all seven morphology scenarios.
Support was identical between fitters because support is determined before the
fitter-specific angle estimate.

Clean-region median errors:

| Scenario | Von Mises | LsqFit | Lower median |
|---|---:|---:|---|
| many filaments | `2.63 deg` | `2.31 deg` | LsqFit |
| crossing filaments | `2.18 deg` | `1.25 deg` | LsqFit |
| curved fibers | `1.39 deg` | `1.22 deg` | LsqFit |
| velocity components | `0.81 deg` | `1.06 deg` | Von Mises |
| shells/arcs | `1.36 deg` | `1.51 deg` | Von Mises |
| Planck-aligned plus contaminant, total mixture | `0.71 deg` | `0.92 deg` | Von Mises |
| sparse CO-like filaments | `6.42 deg` | `6.53 deg` | Von Mises |

Interpretation: this MW-like focused comparison does not show a decisive
fitter-only accuracy improvement from switching back to LsqFit. LsqFit is
competitive and wins several clean, coherent filament medians, especially the
crossing-filament clean subset. Von Mises is better in the sparse CO-like,
velocity-component, shells/arcs, and Planck-contaminant cases, and it also has
better sparse-CO tail behavior (`p90 30.04 deg` versus `36.20 deg`). Keep
`:axial_vonmises` as the production default; LsqFit is not a rollback target
(it was subsequently removed).

## Removed Ideas

The following ideas were tested and removed from the active candidate list:

- `min_valid_70`: mostly produced no accuracy gain. It improved one sparse-fill
  synthetic case slightly but did not help uniform noisy or varying-angle cases.
- `rbar2_gate_0.35`: inconsistent. It helped some cases but failed the sparse
  fill case and was weaker than the simpler `rbar2_gate_0.50` candidate.

These should stay out of production unless future validation gives a clearer
reason to revisit them.

## Noise-Threshold Experiment

`scripts/evaluate_noise_thresholds.jl` now compares scalar, channel, spatial,
full separable, `separable_shrinkage` with `a = b = 0.5`, and an oracle
threshold on deterministic foreground-style synthetic cubes.

Final run summary:

- Stationary noise: scalar wins.
- Channel-varying noise: channel wins.
- Spatially varying noise: spatial wins.
- Mixed channel/spatial noise: full separable wins.
- `separable_shrinkage` wins 0/4 synthetic noise scenarios in this run.

Decision: keep `noise_mode = "scalar"` as the reproducible default. Keep
`"channel"`, `"spatial"`, `"separable"`, and `"auto"` as validation modes.
`"separable_shrinkage"` was not adopted and has since been removed; use the
simplest model whose coverage and residuals validate on the target cube.

## Boundary Caveat

Two-region boundary windows are not ordinary accuracy tests. They can contain
two valid orientation populations, and the circular mean is only a diagnostic
reference, not a unique physical truth. Boundary rows should remain separate
from claims about fitter correctness.

## Recommended Next Validation

Before adopting any accuracy candidate:

- Run `scripts/validate_vgt_parameters.jl` on representative H I and CO crops
  with `:axial_vonmises`.
- Compare reliability-stratified residuals against foreground starlight
  polarization.
- Check Planck 353 GHz foreground-subtraction residuals by reliability bin.
- Require that any `rbar2` weighting or gate improves accuracy without hiding
  scientifically important low-coherence regions.
