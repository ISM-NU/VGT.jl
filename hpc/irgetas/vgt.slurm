#!/bin/bash
#SBATCH --job-name=vgt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=validation_output/vgt_irgetas_logs/%x-%j.out
#SBATCH --error=validation_output/vgt_irgetas_logs/%x-%j.err
#SBATCH --mail-type=END,FAIL

# Run the VGT pipeline on one cube (or several, see below). Pass the config with
#   sbatch --export=ALL,VGT_REPO=$PWD,VGT_CONFIG=<config.toml> hpc/irgetas/vgt.slurm
# Threads come from --cpus-per-task; the chunk memory budget follows --mem
# (src/config.jl memory_budget_bytes). Give each run its own output_dir:
# chunk files are not run-specific.

set -euo pipefail

: "${VGT_REPO:?set VGT_REPO to the repository root}"
: "${VGT_CONFIG:?set VGT_CONFIG to a VGT config.toml}"

cd "$VGT_REPO"
mkdir -p validation_output/vgt_irgetas_logs
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$VGT_REPO/.julia-irgetas}"
threads="${SLURM_CPUS_PER_TASK:-1}"
export JULIA_NUM_THREADS="$threads"

echo "job_id=$SLURM_JOB_ID host=$(hostname) threads=$threads mem_mb=${SLURM_MEM_PER_NODE:-unset}"
echo "config=$VGT_CONFIG"
date --iso-8601=seconds

command -v julia >/dev/null || {
    echo "Julia is not on PATH. Follow hpc/irgetas/README.md before submitting." >&2
    exit 2
}

# VGT_CONFIG may list several configs separated by spaces; they run in one Julia process.
read -ra configs <<< "$VGT_CONFIG"
srun --ntasks=1 julia --project=. --threads="$threads" main.jl "${configs[@]}"

date --iso-8601=seconds
