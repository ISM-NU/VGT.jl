# Provenance of the VGT implementation

Status: active reference.

## The method

The Velocity Gradient Technique (VGT) was developed by A. Lazarian's group.
It traces the plane-of-sky magnetic-field orientation with gradients of
spectral-line velocity information, rotated by 90 degrees:

- González-Casanova & Lazarian 2017, ApJ 835, 41 (`2017ApJ...835...41G`):
  velocity-centroid gradients as a magnetic-field tracer.
- Yuen & Lazarian 2017, ApJL 837, L24 (`2017ApJ...837L..24Y`): sub-block
  averaging of gradient orientations and its application to H I data.
- Lazarian & Yuen 2018, ApJ 853, 96 (`2018ApJ...853...96L`): gradients of
  thin velocity channel maps.
- Lazarian, Yuen, Ho et al. 2018, ApJ 865, 46 (`2018ApJ...865...46L`):
  pseudo-Stokes parameters and the distribution of gradient orientations.

Cite these papers for the method. This repository is an implementation of it,
not a new method.

## This code

This package is a new implementation written at Nazarbayev University for the
LMC sightline. It does not contain code from the Lazarian group.

- It was written from the papers above and from the behaviour of the group's
  earlier prototype, which combined Lazarian-group modules (`LazVGT.jl`,
  `LazCore.jl`) with team scripts.
- That prototype is not part of this repository and is not distributed here,
  because the Lazarian-group modules are not ours to publish. It is kept in a
  private repository for parity checks only.
- Function names that the prototype also used (for example `getcoor`,
  `imfilter_gaussian`) are kept so that older analysis scripts keep
  working through `VGT.getcoor` or `using VGT: getcoor`; their
  implementations are new.

## What differs from the prototype

The production path differs on purpose. Each difference was validated before
it became a default (see `synthetic-fitter-validation-notes.md`):

- an unbinned axial von Mises window estimator instead of a Gaussian fit to an
  angle histogram;
- checked central-difference gradients instead of backward differences;
- exact per-pixel windows, chunked and threaded;
- Stokes Q/U accumulation with doubled angles throughout, `atan(y, x)` instead
  of ratio angles;
- scalar, channel, spatial and separable noise models;
- 2D product headers that record the VGT parameters and the angle convention.

## Reproducing the prototype

The legacy-matching modes and the Gaussian histogram baseline were removed for
the public release. The legacy-matching modes validated worse than the
defaults. The Gaussian histogram fitter is not uniformly worse: the 2026-11
estimator audit (branch `estimator-audit`) found it worse than the von Mises
estimator in windows with a single structure but better within about 9 px of
filament crossings, in line with the crossing-filament rows of
`synthetic-fitter-validation-notes.md`. The removed code, and the validation
runs that compared it, live in the archived private repository.

## Third-party code

None is vendored here. Dependencies are Julia packages listed in
`Project.toml`.
