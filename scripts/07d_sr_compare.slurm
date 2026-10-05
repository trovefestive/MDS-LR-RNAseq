#!/usr/bin/env bash
#SBATCH --job-name=07d_sr_compare
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=results/logs/07d_sr_compare_%j.log

# =============================================================================
# 07d_sr_compare.slurm — Phase 7: short-read vs long-read comparison (runs 17_phase7_compare.R; Q157R − WT)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/17_phase7_compare.R"
echo "[$(date)] 07d complete"
