#!/usr/bin/env bash
#SBATCH --job-name=01b_qcsum
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=results/logs/01b_qc_summary_%j.log

# =============================================================================
# 01b_qc_summary.slurm — Phase 1: run 10_phase1_qc_summary.R after the 01a array
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
echo "[$(date)] $(R --version | head -1)"
Rscript "${PROJ}/scripts/10_phase1_qc_summary.R"
echo "[$(date)] 01b complete"
