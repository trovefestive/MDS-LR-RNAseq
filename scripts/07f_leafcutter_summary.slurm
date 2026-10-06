#!/usr/bin/env bash
#SBATCH --job-name=07f_lc_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=results/logs/07f_leafcutter_summary_%j.log

# =============================================================================
# 07f_leafcutter_summary.slurm — Phase 7b: LeafCutter summary, rule-8 checks, 3′SS tests, comparisons
#   Runs 19_phase7b_leafcutter.R (needs 07e outputs). deltapsi = Q157R − WT.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/19_phase7b_leafcutter.R"
echo "[$(date)] 07f complete"
