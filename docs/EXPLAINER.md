# How This Code Works — A Beginner's Guide

This document explains what `main.jl` does, step by step, in plain language.
No prior knowledge of programming or astrophysics is assumed.

---

## What Is the Goal?

The code analyzes a radio observation of a molecular cloud — a giant cloud of gas
in space, mostly made of carbon monoxide (CO). The observation is stored as a 3D
data file (a FITS "data cube").

The goal is to figure out **which direction the magnetic field is pointing** across
the face of that cloud, and save the result as a 2D map.

The technique used is called **VGT — Velocity Gradient Technique**. It is based on
the idea that turbulence in a magnetized gas tends to produce velocity variations
that are mostly *perpendicular* to the local magnetic field. By measuring how the
gas velocity changes from pixel to pixel, you can infer the field's orientation.

The final output is a 2D image called **ψ (psi)**, where each pixel contains an
angle (in radians) representing the inferred magnetic field direction at that
sky position.

---

## The Input Data: What Is a FITS Cube?

The input file (`CO_data_reproject.fits`) is a **3D array** of numbers:

```
Axis 1 (x): sky position — left/right (Right Ascension)
Axis 2 (y): sky position — up/down (Declination)
Axis 3 (v): velocity — how fast the gas is moving toward or away from us
```

Each "slice" along the velocity axis is called a **channel map**. Each pixel in a
channel map tells you how much CO gas is moving at that particular velocity at that
sky location. This is how radio telescopes observe gas — not just where it is, but
at what velocity, which traces the motion and structure of the cloud.

---

## The Physics in Brief

### Magnetic Fields and the VGT

In turbulent magnetized gas, the fluid tends to flow more freely *along* field
lines than across them. This creates velocity gradients (changes in velocity from
one point to the next) that are statistically stronger *perpendicular* to the
magnetic field.

So if you measure the direction in which velocity changes the fastest — the
gradient direction — you get a direction that is **90° rotated from the magnetic
field**. The code perfonoise this rotation at the end (the `+ π/2` in the `psi`
formula).

### Stokes Q and U

The magnetic field orientation is a *headless* direction — a field pointing left
is the same as one pointing right. Angles like this cannot be averaged directly
because, for example, the average of 1° and 179° should be 90°, not 90°, but naive
arithmetic gives 90° for the wrong reasons and breaks down near 0°/180°.

The standard solution (used in polarimetry) is to encode the angle `θ` as two
components:

```
Q = cos(2θ)
U = sin(2θ)
```

These are the **Stokes parameters**. They can be averaged and smoothed normally.
The angle is recovered at the end with `atan(U, Q) / 2`.

---

## Step-by-Step Walkthrough

### Step 0: Start Julia and Load Packages

```julia
using DrWatson
quickactivate(@__DIR__, "VGT")
```

DrWatson finds the project root and activates the `Project.toml` environment.
Run the first-time setup command in the README once to download missing
packages.

---

### Step 1: Start Julia With Threads

```bash
julia main.jl
```

This starts one Julia process with a fixed thread pool. The pipeline reads one
chunk of the cube at a time to control memory use, and the per-pixel fitting
inside that chunk uses the thread pool. With the repository `.envrc` loaded,
Julia defaults to 16 threads; override with `JULIA_NUM_THREADS=N` when needed.

---

### Step 2: Configure the Run

Before choosing `noise`, estimate it from your cube:

```bash
julia scripts/analyze_inputs.jl path/to/cube.fits
```

The script prints:

- `Suggested config: noise = ...` for `config.toml`
- channel-to-channel Noise spread
- spatial patch Noise spread
- a recommendation on whether a scalar `noise` is adequate or a dynamic threshold should be investigated

If there are no negative pixels, the script falls back to a line-free/blank-region
robust scatter estimate and prints a warning. Prefer the original unclipped
calibrated cube when available.

```toml
[noise]
noise = 0.27
noise_mode = "scalar"

[vgt]
dn = 20
ker = 4

[execution]
ny_chunks = 32
```

These are the tuning knobs:

- **`noise`**: the configured noise threshold. The VGT mask uses this value as the noise
  threshold and ignores pixels below it.
- **`noise_mode`**: how the Noise field is estimated. `"scalar"` uses one value for
  the whole cube. `"channel"`, `"spatial"`, `"separable"`,
  and `"auto"` are validation modes for cubes whose Noise varies strongly with
  velocity or position; use them only if validation on the target cube shows
  they improve coverage and residual behavior.
- **`dn`**: the width of the local square window used when estimating a
  dominant gradient angle. A larger window averages over more pixels (smoother
  but lower resolution).
- **`ker`**: how much to blur the final Q/U maps before computing the angle. This
  reduces pixel-to-pixel noise at the cost of angular resolution.
- **`ny_chunks`**: the data is too large to hold in memory all at once, so it is
  divided into horizontal strips (chunks) that are processed one by one.

---

### Step 3: Read the File Header

```julia
f = FITS(input_file)
dims = size(f[1])
nx, ny, nv = dims[1], dims[2], dims[3]
header = read_header(f[1])
```

The FITS file has two parts: the **header** (metadata — telescope name, sky
coordinates, units, etc.) and the **data** (the actual numbers). Here we only read
the header to learn the dimensions of the cube: `nx` pixels wide, `ny` pixels tall,
`nv` velocity channels deep.

We do not load the actual data yet — it is too large to fit in RAM all at once.

---

### Step 4: Process in Chunks With Threads

```julia
for j in 1:ny_chunks
    process_chunk(j, ...)
end
```

Each chunk number `j` (from 1 to `ny_chunks`) is processed in order. Inside
`process_chunk`, the expensive local-window fits are spread over Julia threads:

#### Inside `process_chunk`:

1. **Open the file and read only one strip**
   The pipeline reads just the rows for the current strip — plus a small overlap
   ("padding") on each side. The overlap is needed because the gradient
   calculation at the edge of a strip needs to see a few pixels from the
   neighboring strip.

2. **Compute the velocity gradients**
   `VGT_QU_mv_2d` loops over every velocity channel. For each channel, it computes
   how fast the intensity changes from pixel to pixel in the x and y directions:
   ```
   grad_x[i,j] = intensity[i-1, j] - intensity[i, j]   <- x-gradient
   grad_y[i,j] = intensity[i, j-1] - intensity[i, j]   <- y-gradient
   ```
   The angle of the gradient vector `atan(grad_y, grad_x)` is the direction of steepest
   change in velocity at that pixel. The first row/column are marked invalid in
   the production mode instead of wrapping around the image edge.

3. **Estimate the dominant local angle**
   Within each local window (a `dn x dn` pixel box), the gradient angles from all
   valid pixels are collected. The production default, `:axial_vonmises`, uses
   doubled-angle circular statistics to estimate one headless orientation and
   the `rbar2` coherence diagnostic. It is the only fitter; the older Gaussian
   histogram fit was removed after validation found it worse.

4. **Convert to Stokes Q and U**
   The dominant gradient angle `θ` is converted:
   ```
   Q += intensity × cos(2θ)
   U += intensity × sin(2θ)
   ```
   The intensity weighting means brighter pixels contribute more to the average,
   which is physically appropriate.

5. **Trim the padding and save to disk**
   The overlap rows are trimmed off, and the result — a 2D Q map and a 2D U map
   for this strip — is saved to a small file (`chunk_2d_j.jld2`).

---

### Step 5: Stitch the Chunks Back Together

```julia
for j in 1:ny_chunks
    q_stitched[:, y_start:y_end] = q_chunk_accum
    u_stitched[:, y_start:y_end] = file["u"]
    intensity_stitched[:, y_start:y_end] = file["i"]
end
```

The saved strip files are loaded one by one and placed into the correct rows of
three full-size 2D arrays: `q_stitched`, `u_stitched`, and
`intensity_stitched` (total intensity).

---

### Step 6: Smooth and Mask

```julia
q_smooth = imfilter_gaussian(q_stitched, ker)
u_smooth = imfilter_gaussian(u_stitched, ker)

q_smooth[intensity_stitched .== 0] .= NaN
u_smooth[intensity_stitched .== 0] .= NaN
```

**Gaussian smoothing** replaces each pixel with a weighted average of its
neighbors, where nearby neighbors count more than distant ones (a bell-curve
weight). This reduces pixel-level noise. The radius of the smoothing is `ker = 4`
pixels.

**Masking**: pixels where the total intensity is zero have no signal — they are
outside the cloud or purely noise. Those pixels are set to `NaN` ("Not a Number"),
which causes them to appear as blank/transparent in plots and be ignored in further
calculations.

---

### Step 7: Compute the Polarization Angle ψ

```julia
psi = 0.5 .* atan.(u_smooth, q_smooth) .+ (pi / 2)
```

This is the key formula:

1. `atan(U, Q)` recovers the doubled angle `2θ` from the Stokes parameters.
2. Dividing by 2 (`0.5 ×`) gives back the gradient angle `θ`.
3. Adding `π/2` (90°) rotates from the gradient direction to the magnetic field
   direction — because the field runs *perpendicular* to the gradient.

The result `psi` is a 2D array where each pixel holds the inferred magnetic field
angle (in radians) at that sky position.

---

### Step 8: Save Results

```julia
write_pipeline_outputs(run)
```

The `psi` map, reliability map, and fitter diagnostics are written as FITS
files. The WCS (World Coordinate System) header from the original data cube is
copied in so that the output images know where they are on the sky — this lets
you overlay them on the original data in astronomy visualization tools like DS9.

---

### Step 9: Plot

```julia
fig, ax, hm = heatmap(psi, colormap=:jet)
```

A quick sanity-check heatmap is saved as `psi_plot.png`. Each pixel is colored by
its angle value. Makie `heatmap` uses the first matrix dimension as x and the
second as y, matching the in-memory FITS map orientation used by the pipeline.

---

## Output Files

By default, results are saved in a parameter-named directory under
`data/sims/vgt/`. If `output_dir` is set explicitly, they are written there.

| File | Contents |
|---|---|
| `chunk_2d_1.jld2` ... `chunk_2d_<ny_chunks>.jld2` | Intermediate per-strip Q, U, I arrays |
| `psi_<input-stem>.fits` | Final 2D magnetic field angle map (in radians) |
| `psi_reliability_<input-stem>.fits` | Stokes coherence/reliability map |
| `psi_fit_rbar2_<input-stem>.fits` | Fitter coherence diagnostic |
| `psi_support_<input-stem>.fits` | Fraction of velocity channels with successful local fits |
| `psi_plot_<input-stem>.png` | Quick-look color image of the angle map |

---

## Summary Diagram

```
data/exp_raw/CO_data_reproject.fits   (3D cube: x, y, velocity)
         │
         ▼
   Divide into y-strips
         │
         ▼ (one chunk at a time; threaded fits inside each chunk)
   For each strip:
     • Read only that strip from disk
     • For each velocity channel:
         – Compute x/y intensity gradients
         – Estimate the dominant local gradient angle
         – Accumulate into Q, U (Stokes parameters)
     • Save strip Q, U, I to disk
         │
         ▼
   Stitch strips → q_stitched, u_stitched, intensity_stitched
         │
         ▼
   Gaussian smooth Q, U, I  (reduce noise)
         │
         ▼
   Mask zero-intensity pixels → NaN
         │
         ▼
   psi = 0.5 × atan(U, Q) + π/2   (magnetic field angle)
         │
         ▼
   Save psi_<input>.fits, reliability, diagnostics, and psi_plot_<input>.png
```
