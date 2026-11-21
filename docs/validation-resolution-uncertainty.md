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

Amendment (2026-11-21, after the smoke test and before any full scoring run): sigma_psi is estimated from K = 16
maps, so even a perfectly calibrated sigma gives error / sigma_psi distributed as Student's t with K - 1 = 15 degrees
of freedom (after the sqrt(K / (K - 1)) correction now applied in `psi_noise_mc`): 0.667 within 1 sigma, 0.936 within
2 sigma, rms z 1.074. The calibration windows of rule 3 are therefore centred on these values: 0.62-0.72 within
1 sigma, 0.91-0.96 within 2 sigma, rms z 0.95-1.25.

Amendment 2 (2026-11-21, after the full transfer-function run and before the sweep was rerun): the effective FWHM
of test 1 is the half-transfer one, FWHM = lambda_half * 2 ln 2 / pi (the Gaussian with the same half-transfer
wavelength), replacing the Gaussian fit. Reason, seen in the data: the box window has a sinc response, and beyond its
first zero the sidelobes give T = 0.1-0.2 (dn 32 at lambda 21 px: 0.19, dn 48 at 32 px: 0.21), which fall inside the
0.1-0.95 fit range and pull the fit low. The fit was non-monotonic in dn (dn 20 / ker 1 12.0 px, below dn 16 / ker 1
12.8 px; dn 32 / ker 4 16.4 px without noise, 26.5 px at r = 0.5), while lambda_half is monotonic and stable with
S/N. The half-transfer wavelength was one of the two declared outputs; only the choice of primary changed. The
Gaussian fit stays in `tf_fit.csv` as `fwhm_gaussfit_px`. Because part sweep smooths each map to the targets F from
this FWHM and rule R2 selects on it, the sweep was rerun with the half-transfer values; the first sweep's R2 lines
are superseded (its native AM scores do not depend on the FWHM and stand).

## What the noise uncertainty depends on (first-order chain)

[interpretation, standard error propagation; the Monte Carlo is the measurement] Each step adds its own noise term:

1. Gradient: central differences of a channel with noise rms sigma_n and lag-2 noise correlation rho_2 have noise
   rms sigma_g = sigma_n sqrt((1 - rho_2) / 2) per component (0.48 sigma_n for a 3 px beam, rho_2 = 0.54). The
   gradient angle of pixel i then scatters by about sigma_g / |g_i|, saturating at uniform angles when |g_i| is
   comparable to sigma_g.
2. Window: the window angle moves by sum_i cos 2(phi_i - Phi) d(phi_i) / (n Rbar), so its noise falls as
   1 / (Rbar sqrt(n_eff)), with n_eff = n divided by the number of pixels a noise blob spans (about 10 for a 3 px
   beam, so a 20 x 20 window holds about 40 independent noise samples), and grows when the gradients in the window
   disagree (low Rbar).
3. Channels: the intensity-weighted sum over channels averages independent channel noise; adjacent-channel noise
   correlation reduces the gain.
4. Smoothing: the Gaussian `ker` averages over neighbouring windows, which share most of their pixels, so it adds
   little beyond the box (variance ratio about dn^2/12 / (dn^2/12 + (0.88 ker)^2), 0.73 for dn 20, ker 4).
5. Mask: at the 3 sigma threshold, noise switches faint voxels in and out of the fit; this is non-linear and is
   why the error is measured by Monte Carlo rather than from the formula.

So sigma_psi is largest where the gas is faint (small |g| / sigma_n), where gradients in a window disagree (low
coherence), and where few channels carry signal. It shrinks with dn at fixed S/N, which is the resolution cost:
the first-order chain suggests 1/dn, but the measured noise error falls about as dn^-0.6 (random fields, r = 1: rms
12.8, 9.4, 7.1 deg for dn 12, 20, 32), probably because the fixed `ker` smoothing adds relatively more averaging to
small windows [interpretation].

## Results

[data] Full runs finished 2026-11-21: transfer function 3 seeds, sweep 4 seeds (rerun after amendment 2), Monte
Carlo 3 seeds, K = 16. Synthetic and MHD only; the LMC application is in the experiment repository.

- Effective FWHM (half-transfer, r = Inf; lower value at half the low-k plateau of about 0.94): dn 8 / ker 4
  10.0-10.4 px, dn 12 / ker 4 11.8-12.2, dn 16 / ker 4 13.6-14.1, dn 20 / ker 1 14.1-14.7, dn 20 / ker 4
  16.1-16.7, dn 32 / ker 4 23.9-25.0, dn 48 / ker 4 34.2-36.4. The formula is within 10 % everywhere (within 3 % at
  dn 20 / ker 4), so it passes the 15 % test. S/N changes the FWHM by less than 1 px except dn 48 at r = 0.5.
- Native AM (ker 4, dn 8/12/16/20/32/48): random fields r = 1 0.859/0.890/0.911/0.923/0.934/0.928, r = 2
  0.934/0.940/0.945/0.946/0.940/0.925; CATS b1p1 r = 1 0.738/0.787/0.830/0.865/0.928/0.959. On native AM dn 32 /
  ker 4 passes the adoption rule at r = 1 (random fields +0.0102, CATS b1p1 +0.063) and dn 32 / dn 48 at r = 0.5;
  nothing passes at r >= 2. The trade is resolution (25 against 16.7 px).
- R2 at 21 and 30 px: dn 20 / ker 4 is chosen or tied on the random fields at every r and on CATS b1p1 at
  r = 0.5-2; the adoption rule finds nothing at r >= 1. At 10.5 and 15 px dn 20 / ker 4 is not eligible; the declared
  tie-break (smallest FWHM) picks dn 8-12 settings on the random fields and dn 16 / ker 4 on CATS b1p1 at 15 px (r >= 1).
- Calibration (dn 20, scale 1, within 1 / 2 sigma): random fields r = 4 0.668/0.938, r = 2 0.747/0.968, r = 1
  0.845/0.988; CATS b1p1 r = 4 0.484/0.776, r = 2 0.482/0.787, r = 1 0.663/0.917. The amended window is met in 2
  of the 6 regimes (1 of 6 under the original 0.92-0.98 two-sigma window, which CATS r = 1 misses at 0.917); rms z
  is in 0.95-1.25 in the same 2. sigma_psi is conservative on the random fields at r <= 2 and about 1.5x too small
  on CATS at r >= 2. Scale 0.5 is closer by the summed distance from 0.667 (0.61 against 0.63) and in rms z, so by
  rule 3 it is the default.


Smoke test (1 seed, K = 4 realisations, two wavelengths; not a result, a check that the pipeline works):
- dn 20, ker 4: FWHM 16.7 px measured against 16.3 px from the formula; dn 20, ker 1: 14.5 against 14.1; low-k
  transfer at lambda = 128 px 0.92.
- CATS b1p1 at r = 1: AM 0.88 for dn 20 / ker 4, 0.80 for dn 12 / ker 4, 0.71 for dn 12 / ker 1.
- Coverage at r = 1, dn 20 (K = 4, before the sqrt(K / (K - 1)) correction, where a calibrated sigma gives 0.55
  within 1 sigma and rms z is undefined): random fields 0.71, CATS 0.56-0.61. These are pipeline checks only; the
  CATS value is what a calibrated sigma predicts at K = 4. The full run also scores sqrt(sigma^2 + offset^2)
  (Monte Carlo mean shift) as a post hoc variant.
