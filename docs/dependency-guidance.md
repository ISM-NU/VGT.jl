# Julia Dependency Guidance

Status: active guidance for dependency choices.

This repo is a research codebase. Add Julia packages when they reduce
maintenance burden, improve reproducibility, or make validation more explicit.
Do not add packages merely to make the code look modern, and do not replace
scientific behavior without component-matched validation evidence.

## Installed Now

### `WCS.jl`

Use for FITS World Coordinate System handling and coordinate metadata checks.
`getcoor` now uses `WCS.from_header` and `WCS.pix_to_world` with a linear
fallback. Use this for coordinate transforms, footprint checks, and explicit
pixel/world assumptions.

Guardrails:

- preserve FITS pixel-index conventions carefully; `WCS.jl` uses WCSLIB-style
  coordinates, so test one-indexed Julia pixel assumptions explicitly;
- do not change output map values while adopting WCS metadata helpers;
- add small tests using synthetic FITS headers before touching production
  foreground products.

### `DrWatson.jl`

Use for parameter-named run directories and experiment organization. The main
pipeline now resolves legacy `subcube`/empty output directories to
`data/sims/vgt/<parameter-name>/`, validation/report scripts default to
`data/sims/...`, and `validate_vgt_parameters.jl` writes parameter-named
candidate directories under `data/sims/vgt_parameter_scan/` by default. Scripts
that produce scalar validation or profiling metrics also write tagged summary
records under `<output-dir>/summaries/` and refresh `<output-dir>/collected_results.jld2`.
The collector uses DrWatson's result collection API when available and otherwise
falls back to a small JLD2 aggregate compatible with the pinned dependency set.

Guardrails:

- keep generated `data/` outputs ignored by Git;
- commit durable conclusions under `docs/validation/`, not full generated run
  directories;
- include every science-affecting `VGTConfig` switch in run parameters so
  incompatible outputs do not share a directory;
- keep writing a per-run metadata sidecar with full config, diagnostics, output
  paths, and DrWatson git provenance.
- keep summary records scalar and compact; do not make DrWatson aggregation load
  full FITS-scale arrays or candidate matrix caches.

### `AstroFITS.jl`

Use through the `VGT.read_fits_header`, `VGT.read_fits_image`, and
`VGT.write_fits` wrappers. The wrappers prefer AstroFITS for full-image
reads, plain header reads, and simple image writes, while retaining FITSIO for
sliced cube reads and fallback writes.

Guardrails:

- do not call AstroFITS directly from pipeline code; use the wrappers;
- keep FITSIO until sliced reads and header-preserving writes are fully covered;
- add tests for every new FITS wrapper behavior before migrating call sites.

### `BenchmarkTools.jl`

Use for repeatable performance experiments and regression checks. Prefer it for
new benchmark scripts over ad hoc timing when measuring candidate changes to
window scanning, fitters, interpolation, or chunking.

## Install Only When Needed

### `Reproject.jl`

Candidate for Planck/VGT footprint matching and component-labelled comparison
workflows. Use only when a script needs explicit astronomical reprojection
between FITS/WCS grids.

### `Distributions.jl`

Candidate for reference checks involving `VonMises` distributions or synthetic
sampling. Do not replace the production axial estimator unless doubled-angle,
headless-orientation behavior is explicit and validated.

### `DirectionalStatistics.jl`

Candidate for comparison tests against independent circular-statistics
implementations. Treat it as a validation aid, not a production dependency,
unless adoption is justified by tests and maintenance value.

## Dependency Review Checklist

Before adding another package:

1. State whether it supports metadata, validation, reproducibility, performance,
   or a scientifically motivated algorithm change.
2. Add or update a doc pointer so later agents know why it exists.
3. Add tests or validation commands that cover the package boundary.
4. Keep optional experiment packages out of production code until their value is
   demonstrated on component-matched foreground workflows.
