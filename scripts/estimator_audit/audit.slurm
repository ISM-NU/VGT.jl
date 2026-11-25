#!/bin/bash -l
#SBATCH --job-name=est_audit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=0-13
#SBATCH --output=validation_output/estimator_audit/logs/%x-%A_%a.out
#SBATCH --error=validation_output/estimator_audit/logs/%x-%A_%a.err

# Estimator audit grid: array task -> (suite, scenario), seeds AUDIT_SEEDS.
#   sbatch --export=ALL,AUDIT_OLD_FITTERS=<private fitters.jl> scripts/estimator_audit/audit.slurm
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mw=(many_filaments crossing_filaments curved_fibers velocity_components shells_arcs planck_misaligned sparse_co_like)
uni=(uniform_n0.05 uniform_n0.12 uniform_n0.20 weak_amp0.05 sparse_fill0.35 varying_c8e-4 varying_c2e-3)
i=${SLURM_ARRAY_TASK_ID}
if (( i < 7 )); then suite=mw; scen=${mw[$i]}; else suite=uni; scen=${uni[$((i - 7))]}; fi
export AUDIT_SUITE=$suite AUDIT_SCEN=$scen
export AUDIT_SEEDS="${AUDIT_SEEDS:-0:9}"
export AUDIT_EST="${AUDIT_EST:-vm,gauss,st1,st2,ms2,ms4,st2ms2,vmu,hg100,hg100n,hg36}"
export AUDIT_OUT="validation_output/estimator_audit/${suite}_${scen}"
echo "job=$SLURM_JOB_ID task=$i host=$(hostname) suite=$suite scen=$scen seeds=$AUDIT_SEEDS"
date --iso-8601=seconds
srun --ntasks=1 julia --project=scripts/estimator_audit --threads="${SLURM_CPUS_PER_TASK}" \
    scripts/estimator_audit/run_audit.jl
date --iso-8601=seconds
