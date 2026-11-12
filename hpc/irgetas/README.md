# VGT on NU Irgetas

How we run the VGT pipeline on the Irgetas cluster. Login node: use it only to
copy files and submit jobs.

| Script | What it does |
|---|---|
| `vgt.slurm` | Run the VGT pipeline (`main.jl`) on one config with `--cpus-per-task` Julia threads |

## Access

Connect to the NU GlobalProtect VPN when off campus, then
`ssh YOUR_USERNAME@<irgetas login node>` (address in the NU HPC user guide). Our copy of the repository is
`/shared/home/YOUR_USERNAME/vgt`. It is not a git checkout: sync tracked files
from a laptop checkout with

```bash
rsync -az --files-from=<(git ls-files) \
  ./ YOUR_USERNAME@<irgetas login node>:/shared/home/YOUR_USERNAME/vgt/
```

Large cubes go in `data/exp_raw/` with `rsync -az --partial`, which resumes
after a VPN drop.

`module` only exists in a login shell: use `bash -lc "module load ..."`
over ssh, and `#!/bin/bash -l` in ad-hoc batch scripts. Julia packages live in
the depot `.julia-irgetas` inside the repository (`JULIA_DEPOT_PATH`); set it up
once with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

## VGT runs

```bash
sbatch --job-name=vgt_<name> --export=ALL,VGT_REPO=$PWD,VGT_CONFIG=<config.toml> \
  hpc/irgetas/vgt.slurm
squeue --me
```

VGT threads come from `--cpus-per-task` (default 32) and its chunk memory budget
from `--mem` (default 64 GB). Give each config its own `output_dir`: chunk files
are not run-specific. Logs: `validation_output/vgt_irgetas_logs/`.

ROHSA and experiment jobs live in the `hpc/` folders of their repositories
(`ISM-NU/ROHSA-LMC`, `ISM-NU/LMC-foreground`, `ISM-NU/LMC-cold-gas`).

NU asks users to acknowledge Nazarbayev University Research Computing when its
HPC facilities contribute materially to published results.
