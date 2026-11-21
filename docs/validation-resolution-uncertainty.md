# Effective resolution, window size and noise uncertainty of the VGT angle map

Status: validation evidence, in progress. The rules below were written and committed before any scoring run.
Script: `scripts/validate_resolution_uncertainty.jl` (generator `scripts/synthetic_orientation_fields.jl`);
uncertainty code: `src/science/uncertainty.jl` (`psi_noise_mc`, opt-in, not called by the pipeline).

Labels: [data] measured here, [lit] from a paper, [assumption] taken as given, [interpretation] our reading.

## What VGT does to resolution (from the code)

- [data] The window is evaluated at every pixel: for each channel, the gradient orientations of the `dn` x `dn`
  box centred on the pixel are averaged as doubled-angle unit vectors (`src/science/vgt.jl`), so the windows
  already slide with stride 1. Non-overlapping sub-blocks (the 20 x 20 blocks of Yuen & Lazarian 2017, used by
  Alina et al. 2023) give the same box response on a grid `dn` times coarser; the sliding form keeps the map
  Nyquist-sampled. A smaller stride is not available as an improvement because it is already the production form.
- [data] The window means are weighted by the channel intensity at the centre, summed over channels and smoothed
  with `Kernel.gaussian(ker)` (sigma `ker` px, truncated at +-2 sigma, effective std 0.88 `ker`).
- [assumption, tested below] In the linear regime the angle map is the orientation field convolved with the box,
  the truncated Gaussian and the beam, so FWHM_eff = 2.355 sqrt(dn^2 / 12 + (0.88 ker)^2 + sigma_beam^2) px:
  16.2 px for dn 20, ker 4, a 3 px beam.

## Declared tests and rules (2026-11-21, before scoring)

Synthetic data: 512^2 maps, 8 channels, Gaussian random fields with power k^-8/3 and anisotropy zeta 3 following a
prescribed orientation field, 3 px FWHM beam, beam-correlated noise; S/N r = high-passed (20 px) signal rms over
noise rms, r in {Inf, 2, 1, 0.5}; the 3 sigma mask does not act on these fields (constant level added). MHD: CATS
b1p1 and b.1p1 (Burkhart et al. 2020), line of sight 3, 48 channels, thermal width 1.0, beam 3 px, noise added at the
same r, 3 sigma mask. Windows dn in {8, 12, 16, 20, 32, 48}, ker in {1, 2, 4}, min_valid = dn^2 / 4;
`rbar2_weight = true` at dn 12, 20, 32.

1. Effective resolution: the transfer function of a sinusoidal orientation modulation (15 deg amplitude) against
   wavenumber; FWHM from a Gaussian fit to T between 0.1 and 0.95, and the half-transfer wavelength. Reported for
   every (dn, ker) at every r. The formula above is adopted for the facts tables if it is within 15 % of the
   measured FWHM at r = Inf; otherwise the measured value is used.
2. Window choice. Scores: AM = <cos 2(psi - truth)> over the common interior, (a) against the pixel-scale truth,
   (b) at target resolutions F in {10.5, 15, 21, 30} px (7' and 10' on a 40" grid, 7' and 10' on a 20" grid), where
   the VGT Stokes map is smoothed from its measured FWHM to F and the truth to F.
   - R1 (native map): per regime (field, r), the config with the highest AM (a); configs within max(0.005, 2 standard
     errors over seeds) of the best count as tied and the one with the smallest FWHM is chosen.
   - R2 (map for a comparison at F): among configs with FWHM <= F, the highest AM (b); ties as in R1, with
     dn 20 / ker 4 kept if it is among the tied.
   - Adoption: a config replaces dn 20 / ker 4 only if it beats it by more than 0.01 in AM in the same regime on the
     random fields and on at least one CATS run, without losing more than 0.01 on the other. `rbar2_weight` is
     judged the same way. The choice is made on synthetic and MHD data only; on the LMC the effect of the chosen
     config is reported, never tuned against Planck.
3. Uncertainty: `psi_noise_mc` with K = 16 realisations of the noise model on signal + one noise realisation;
   error = psi(signal + noise) - psi(signal), the noise-free VGT map (method error is not noise). Calibrated if the
   fraction of pixels with |error| <= sigma_psi is within 0.63-0.73 and with |error| <= 2 sigma_psi within
   0.92-0.98, for r >= 1 at dn 20 / ker 4, on the random fields and CATS. Reported for dn 12 and 32 and r = 0.5;
   scale 1 (added noise as large as the cube's) and 0.5 both run, the one closer to calibration is the default.

## Results

Pending (cluster runs).
