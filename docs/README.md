# Documentation Map

Each entry gives the purpose and the status (active guidance, validation
evidence, or reference).

## Method docs

- `../README.md` — setup, configuration, run commands, outputs. Active guidance.
- `EXPLAINER.md` — plain-language VGT workflow and physics. Active guidance.
- `dependency-guidance.md` — Julia package adoption guidance and optional
  library candidates. Active guidance.
- `synthetic-fitter-validation-notes.md` — synthetic-cube evidence for fitter
  defaults (including the removed Gaussian histogram baseline) and follow-up
  accuracy candidates. Validation evidence.
- `validation-gradient-whitening.md` — synthetic test of opt-in gradient
  whitening against correlated, anisotropic noise (and the rejected gradient
  S/N mask). Validation evidence.
- `validation-angle-calibration.md` — grid-locked angle bias of the gradient
  stencils on synthetic fields with known angles, and the opt-in
  `isotropic_checked` stencil. Validation evidence.
- `validation-resolution-uncertainty.md` — effective resolution of the angle
  map against `dn` and `ker`, window choice by declared rules, and the opt-in
  Monte Carlo noise uncertainty `psi_noise_mc`. Validation evidence, in progress.
- `provenance.md` — where the method and this implementation come from, and
  why the legacy prototype and its matching modes are not included. Active
  reference.
- `../hpc/irgetas/README.md` — NU Irgetas cluster setup and the VGT job.
  Active guidance.

## Code layout

- `../main.jl` — thin command-line entrypoint.
- `../src/config.jl` — typed VGT configuration and execution resolution.
- `../src/science/vgt.jl` — core VGT gradients, windows, Q/U accumulation, chunk reads.
- `../src/science/fitters.jl` — axial von Mises-style window estimator.
- `../src/science/noise.jl` — scalar/channel/spatial/separable noise models.
- `../src/science/products.jl` — final `psi`, reliability, and axial residual helpers.
- `../src/pipeline/run.jl` — full run orchestration.
- `../src/pipeline/chunks.jl` — chunk stitching.
- `../src/pipeline/outputs.jl` — FITS/PNG output writing.
- `../scripts/README.md` — script catalog.

Experiments that use VGT (foreground, starlight, Planck, ROHSA, cold gas,
shells) are documented in their own repositories: `ISM-NU/LMC-foreground`,
`ISM-NU/ROHSA-LMC`, `ISM-NU/LMC-cold-gas`, `ISM-NU/LMC-shells`.

Do not add long exploratory notes without linking them here and stating their
status.
