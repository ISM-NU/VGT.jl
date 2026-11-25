# AGENTS.md

## What This Project Is

This repository is a Julia implementation of the Velocity Gradient Technique
(VGT) for spectral-line FITS cubes on the LMC sightline. It estimates
plane-of-sky magnetic-field orientation proxies from H I and CO gas kinematics.

Primary science outputs:

- `psi_*.fits` - 2D angle map in radians.
- `psi_reliability_*.fits` - smoothed Stokes coherence map,
  `P = sqrt(Q^2 + U^2) / W` in [0, 1], where `W` is the summed Stokes weight of
  the fitted windows (definition of 2026-11-25, `docs/decisions.md`).
- `psi_fit_rbar2_*.fits` - fitter-level doubled-angle coherence diagnostic. It
  stores the mean resultant length R̄, not R̄²; the name is kept.
- `psi_support_*.fits` - fraction of velocity channels with successful local
  VGT fits per sky pixel.
- `psi_plot_*.png` and `run_metadata.toml` - quick-look plot (skipped for an
  all-NaN map) and run provenance.

Primary local data locations:

- `data/exp_raw/Tb_combined_LMC_182_339_res1_3pperb_reprojected_on_gass_masked.fits`
  - H I cube on the LMC footprint.
- `data/exp_raw/lmc_v11.co.vradlsr.base.fits` - MAGMA CO cube.
- `data/sims/vgt/...` - default DrWatson output tree for pipeline products.

The large FITS cubes are local inputs and are not expected to be tracked by Git.

## Research Standard

This is a research codebase. The angle map is a scientific result, not just a
software artifact.

Critical thinking comes first:

- Do not change algorithms, defaults, masks, interpolation, fitter behavior, or
  numerical shortcuts unless the change is scientifically motivated and
  validated.
- Do not do actions for their own sake. Every code change, doc change, benchmark,
  or cleanup must have a clear research or maintainability reason.
- Do not optimize runtime at the cost of accuracy, interpretability, or
  reproducibility.
- Do not assume a prettier map is a better map.
- Do not treat Planck or foreground starlight polarization as ground truth.
  They are comparison tracers with their own assumptions.
- When uncertain, preserve the current scientific behavior and add validation
  before changing production outputs.

The right question is not "can we do this?" but:

> Does this improve the reliability, interpretability, or validated accuracy of
> the Milky Way foreground estimate, Planck foreground subtraction, or
> starlight-polarization validation workflow?

## Current Scientific Interpretation

- `:axial_vonmises` is the production default. It is an unbinned axial circular
  estimator on doubled angles.
  It is the only fitter. The Gaussian-histogram fitter was removed in 0.1.0.
  The 2026-11 estimator audit (branch `estimator-audit`) found it worse than
  `:axial_vonmises` in windows with a single structure but better within about
  9 px of filament crossings. The legacy-matching modes and
  `separable_shrinkage` noise were removed because validation found them worse.
  All of them live in the archived private repository (see
  `docs/provenance.md`). Config files that still set them fail with an
  explanatory error.

The production CPU path uses exact per-pixel windows, checked central-difference
gradients, and Stokes Q/U accumulation. The former backward-difference path is
retained as `gradient_mode = "nan_edge_backward"` in TOML config for comparisons. Preserve these
semantics unless validation justifies a documented change.

## Scope And Language

- This repository holds the VGT method only: the package, its tests, and
  validation that runs on synthetic data. Experiments on LMC cubes belong in
  ISM-NU/LMC-foreground, ISM-NU/LMC-cold-gas or ISM-NU/LMC-shells; ROHSA belongs
  in ISM-NU/ROHSA-LMC.
- Use Julia where performance needs it; simple scripts can be Python.
- `data/` is a symlink to the shared `~/Projects/LMC/lmc-data` folder.

## Repository Map

Start here:

- `README.md` - setup, configuration, run commands, output files.
- `docs/README.md` - documentation map.
- `docs/EXPLAINER.md` - beginner-oriented VGT workflow and physics.
- `docs/dependency-guidance.md` - dependency adoption and helper-package notes.

Core code:

- `main.jl` - thin entrypoint that loads config and runs the package pipeline.
- `examples/quickstart.jl` - synthetic end-to-end example; regenerates `docs/images/quickstart.png`.
- `src/VGT.jl` - package module and includes.
- `src/config.jl` - typed `VGTConfig` and config normalization.
- `src/science/vgt.jl` - core VGT gradients, windows, Q/U accumulation, chunk
  reads.
- `src/science/fitters.jl` - axial window estimator.
- `src/science/noise.jl` - scalar/channel/spatial/separable noise models.
- `src/science/products.jl` - final `psi`, reliability, and axial residual
  helpers.
- `src/pipeline/run.jl` - full run orchestration.
- `src/pipeline/chunks.jl` - chunk stitching.
- `src/pipeline/outputs.jl` - FITS/PNG output writing and metadata.
- `src/fits_io.jl`, `src/paths.jl` - FITS wrappers and DrWatson path helpers.

Validation and analysis (full catalog in `scripts/README.md`):

- `scripts/analyze_inputs.jl` - input metadata, noise, and fill analysis.
- `scripts/evaluate_noise_thresholds.jl` - synthetic noise-threshold validation.
- `scripts/validate_synthetic.jl` - bounded synthetic VGT/fitter validation.
- `scripts/validate_mw_like_synthetic.jl` - gradient-mode validation on
  Milky-Way-like synthetic morphologies.
- `scripts/validate_vgt_parameters.jl` - `dn`/`min_valid` candidate scans.
- `scripts/compare_vgt_runs.jl` - compare existing run directories.
- `scripts/profile_pipeline.jl` - representative crop profiling.
- `scripts/README.md` - current script catalog.

Historical code:

- The legacy prototype is not in this repository; see `docs/provenance.md`.

## Commands

Run commands from the repository root.

```bash
# Run the configured pipeline
julia --project=. main.jl

# Estimate noise and fill statistics
julia --project=. scripts/analyze_inputs.jl path/to/cube.fits

# Run tests
julia --project=. test/runtests.jl

# Run bounded synthetic validation
julia --project=. scripts/validate_synthetic.jl
```

The checked-in `.envrc` sets `JULIA_PROJECT=@.` and defaults
`JULIA_NUM_THREADS` to `16` for `direnv` users. Override thread count per
command when needed, for example:

```bash
JULIA_NUM_THREADS=8 julia --project=. main.jl
```

## Configuration Rules

- Start from `config.example.toml`; keep local `config.toml` out of Git.
- Put raw FITS cubes under `data/exp_raw/`.
- Leave `output_dir = ""` to use the DrWatson parameter-named run directory
  under `data/sims/vgt/`.
- Keep `noise_mode = "scalar"` unless cube-specific validation supports
  `"channel"`, `"spatial"`, `"separable"`, or `"auto"`.
- Use `ny_chunks = "auto"` unless there is a measured reason to pin chunk count.
  Julia thread count is controlled at launch, not by config.
- Do not enable `rbar2_weight` by default; it needs component-matched validation.

## Accuracy-Critical Rules

- Use stable axial-angle math. Prefer `atan(y, x)` and doubled-angle Stokes
  statistics. Do not use fragile `atan(y / x)` ratios in production paths.
- Average headless orientations with Stokes Q/U. Never average angles directly.
- Preserve CPU gradient semantics: checked central differences are part of the
  current production path; `:nan_edge_backward` is retained for comparisons.
- Keep `psi_reliability` semantics unchanged unless a documented decision says
  otherwise. Before changing it or any other output's meaning, read and extend
  `docs/decisions.md`.
- Treat changes to `dn`, `ker`, `min_valid`, masks, noise thresholds, fitter
  modes, chunk padding, and output support as science changes.
- Any performance change must pass parity or validation tests on
  component-matched project benchmark crops.

## Documentation Rules

Docs must be useful, short, and connected.

- Add or update cross-file pointers whenever a new doc becomes part of active
  guidance.
- Do not add long exploratory text without marking its status: active guidance,
  validation evidence, historical note, or archive.
- Prefer one authoritative document over multiple overlapping explanations.
- If a doc becomes stale, replace it with a short superseded note that points to
  the current source of truth.
- Keep PDFs and papers under `docs/references/` and summarize their relevance in
  project docs instead of scattering notes.
- Do not link to deleted roadmap files. Keep active research direction either in
  `README.md`, this file, or a current doc linked from `docs/README.md`.

## Cleanup Rules

- Generated local products belong in ignored output directories, normally
  `data/sims/...`. Explicit one-off outputs may use ignored scratch directories
  such as `validation_output/` when a script or test expects that shape.
- Durable validation evidence belongs in `docs/` with clear status and links
  from `docs/README.md`.
- Do not delete legacy code, validation evidence, raw data, or user-provided
  papers without explicit permission.
- Do not remove generated intermediates unless the run completed successfully or
  the user explicitly asks for cleanup. The pipeline can remove chunk JLD2 files
  via `keep_intermediate_chunks = false`.

## Research Direction

This repository is the method only. It serves the foreground-subtraction,
starlight-validation and cold-gas experiments, which live in
`ISM-NU/LMC-foreground`, `ISM-NU/LMC-cold-gas` and `ISM-NU/ROHSA-LMC`. Method
changes still need the validation described above.
