#!/bin/bash
#SBATCH --job-name=vgt_behaviour
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=03:30:00
#SBATCH --output=validation_output/behaviour/logs/%x-%A_%a.out
#SBATCH --error=validation_output/behaviour/logs/%x-%A_%a.err

# VGT behaviour study (scripts/validate_behaviour.jl), one experiment per array task:
#   sbatch --array=0-9 --export=ALL,VGT_REPO=$PWD,CATS_DIR=<dir with b1p1/, b1p.32/, b.1p.032/> hpc/irgetas/behaviour.slurm
# QUICK=1 runs every experiment in its smoke-test size in one task. Results are small CSV
# files in validation_output/behaviour/.
# Measured with 16 threads (2026-11-22, 3 seeds): peak memory 1.9-3.3 GB per task; wall time
# noise 20 min, split 18, chains 2 h 29, thr 10, weights 17, support 4, dn 5, edges 26,
# seams 36, wsplit 19.

set -euo pipefail
: "${VGT_REPO:?set VGT_REPO to the repository root}"
: "${CATS_DIR:?set CATS_DIR to the CATS simulation directory}"
cd "$VGT_REPO"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$VGT_REPO/.julia-irgetas}"
threads="${SLURM_CPUS_PER_TASK:-1}"
exps=(noise split chains thr weights support dn edges seams wsplit)
out=validation_output/behaviour
if [ "${QUICK:-0}" = 1 ]; then
    for e in "${exps[@]}"; do
        julia --project=. --threads="$threads" scripts/validate_behaviour.jl --exp "$e" --cats-dir "$CATS_DIR" \
            --out "$out/quick" --quick 1
    done
else
    e=${exps[${SLURM_ARRAY_TASK_ID:-0}]}
    julia --project=. --threads="$threads" scripts/validate_behaviour.jl --exp "$e" --cats-dir "$CATS_DIR" \
        --out "$out" --seeds "${SEEDS:-3}"
fi
