#!/usr/bin/env bash
#SBATCH --job-name=06b_diffiso
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=results/logs/06b_diff_isoforms_%j.log

# =============================================================================
# 06b_diff_isoforms.slurm — Phase 6: edgeR DGE/DTE, DRIMSeq DTU, switches, LRP2 comparison (negated), Ybx1,
# dot plots. Runs 15_phase6_diff_isoforms.R (Q157R − WT; rule-8 assertions inside).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/15_phase6_diff_isoforms.R"
echo "[$(date)] 06b complete"
