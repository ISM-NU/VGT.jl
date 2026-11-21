# Scripts

Validation, benchmarking, profiling and cube-preparation scripts that sit
outside the package API. Experiment scripts live in other repositories (see the
end).

Julia scripts assume they are run from the repository root, for example
`julia --project=. scripts/validate_synthetic.jl`. Raw input cubes belong under
`data/exp_raw/`; generated outputs default to `data/sims/...`. Explicit
`--output-dir` values are accepted for one-off runs.

Scripts that compute scalar validation, comparison, or profiling metrics write
compact DrWatson-tagged JLD2 records under `<output-dir>/summaries/` and refresh
`<output-dir>/collected_results.jld2`. These aggregates are for scanning runs;
they do not replace the script-specific FITS, text, or CSV outputs.

## Pipeline tools

`analyze_inputs.jl`
: Inspects a FITS cube, reports header metadata, estimates 1σ noise from negative
  pixels or a labeled no-negative fallback, summarizes channel occupancy, and
  recommends `noise` and masking choices. `--noise-vrange=LO,HI` also measures
  the gradient noise (central and Sobel)
  in emission-free channels: the gradient covariance `gradient_noise_cov`, the
  anisotropy A and its direction. `--cov-map=OUT.fits [--cov-fwhm=PX]` also writes
  the smoothed per-pixel covariance cube for `gradient_noise_cov_file`.

`channel_quality_report.jl`
: Writes a per-channel CSV with robust noise, above-threshold fill at 1/2/3σ,
  positive-vs-negative summed brightness, and a conservative likely-noise-only
  flag. Use it before VGT on broad H I cubes so line-free channels are handled
  explicitly.

`make_velocity_slabs.jl`
: Splits a FITS cube into velocity slab cubes using header-derived
  `CRVAL3`/`CRPIX3`/`CDELT3` coordinates.

`regrid_cube_to_reference.jl`
: Reprojects a 3-D source cube onto a reference cube's spatial WCS and
  overlapping linear velocity channels (bilinear in space, linear in velocity).
  Row-blocked and memory-mapped, so large cubes do not need a second copy in
  RAM. It aligns sampling; it does not match beams.

`smooth_cube_to_beam.jl`
: Smooths a cube from one Gaussian beam to a larger one and resamples it to a new
  pixel size on the same sky window. NaN pixels are excluded with normalised
  convolution.

`profile_pipeline.jl`
: Profiles the VGT pipeline on a configurable crop and writes a stage-by-stage
  timing report plus a compact summary record.

## Method validation

`validate_synthetic.jl`
: Bounded synthetic sweeps of the `:axial_vonmises` fitter with known
  orientations (noise, sparsity, amplitude, curvature, boundaries) plus
  gradient-stencil and `rbar2` candidates, to check changes before touching
  production outputs.

`validate_mw_like_synthetic.jl`
: Gradient-mode validation on Milky-Way-like synthetic morphologies (filament
  networks, crossings, curved fibers, component mixtures, shells, sparse
  CO-like filaments). Compares `:nan_edge_backward`, `:central_checked` and
  `:sobel_checked` across fractional `min_valid` thresholds using angle error,
  support and `rbar2`. Focused sweeps can be filtered with `VGT_MW_MODES`,
  `VGT_MW_FITTERS`, `VGT_MW_MIN_VALID` and `VGT_MW_SIZE`.

`validate_angle_bias.jl`
: Grid-locked angle bias (VGT minus true angle against the true angle) of the
  central, Sobel and isotropic stencils on anisotropic Gaussian random fields
  and on rotated CATS MHD cubes, with the harmonic fit A4 sin(4t) + A8 sin(8t);
  `--quick 1` is a smoke test. Results: `docs/validation-angle-calibration.md`.

`validate_coherence_debias.jl`
: Noise in the beam coherence |<exp(2i psi)>| of VGT angle maps: raw coherence,
  floor subtraction (`debias_coherence`), cross product and the split estimator
  (`split_coherence`, two independent noise realizations corrected for angle
  attenuation) against the noiseless VGT map, on a ramp-to-random-field mixture
  with isotropic and anisotropic correlated noise. Run on the cluster.

`validate_gradient_whitening.jl`
: Opt-in gradient whitening (`gradient_whiten`) against correlated, slightly
  anisotropic noise: angle bias toward the noise direction, scatter and coverage on
  a ramp and on filament morphologies, with Q/U smoothing as in production;
  also isotropic-noise, mis-estimated-covariance and spatially varying noise
  (global covariance against a per-pixel map) cases. Results in
  `docs/validation-gradient-whitening.md`.

`validate_vgt_parameters.jl`
: Candidate scan for `dn`, `min_valid` and experimental `rbar2` accumulation
  choices. Writes `summary.txt`, `summary.csv` and `pairwise_residuals.csv`.
  Its recommendation is an internal-stability gate, not external validation.

`evaluate_noise_thresholds.jl`
: Synthetic harness for scalar, per-channel, spatial, separable
  and oracle noise thresholds. Reports estimator error, VGT angle
  error and valid-window coverage on cubes with known noise and true angles.

`compare_vgt_runs.jl`
: Compares existing VGT run directories: per-run coverage, reliability and
  `rbar2` summaries, pairwise axial residual CSVs, and
  `support_partitions.csv` (common, gained and lost support against the first
  run). `--write-residual-maps true` adds residual FITS maps.

## MHD-truth tests

These score VGT against a known field instead of starlight.

`synthetic_ppv_from_mhd.jl`
: Builds an optically thin 21 cm cube from an MHD snapshot (CATS) plus the true
  plane-of-sky field angle `psi_true`. Writes `ppv.fits` and `psi_true.fits`.

`score_vgt_vs_true.jl`
: Scores VGT angle maps against `psi_true.fits`: alignment measure
  `AM = <cos 2 dpsi>` and median `|dpsi|` over pixels where both maps are finite.

`vda_velocity_caustics.jl`
: Velocity Decomposition Algorithm: splits every channel map into the part that
  follows the integrated intensity and the velocity-caustics part, so VGT can
  run on either (or on the unchanged crop as baseline).

`gradient_magnetization.jl`
: Magnetization from the anisotropy of velocity-centroid gradients: per-block
  `J2` and the Alfvén Mach number estimate, written with the velocity
  dispersion for the plane-of-sky field estimate.

`figures/figY_vgt_sim_test.py`
: Figure Y: VGT on CATS simulations against the true field (intensity with true
  segments, VGT field, angle error).

`figures/figAA_lmc_j2.py`
: Figure AA: `J2` maps of the LMC cube next to two CATS simulations.

Experiment scripts (foreground, starlight, Planck, ROHSA, cold gas, shells) are
in `ISM-NU/LMC-foreground`, `ISM-NU/ROHSA-LMC`, `ISM-NU/LMC-cold-gas` and
`ISM-NU/LMC-shells`.
