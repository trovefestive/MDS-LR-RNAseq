#!/usr/bin/env bash
#SBATCH --job-name=06c_rmatssum
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=results/logs/06c_rmats_summary_%j.log

# =============================================================================
# 06c_rmats_summary.slurm — Phase 6: rMATS-long import, rule-8 direction checks, DRIMSeq agreement,
# isoform-difference classes. Runs 16_phase6_rmats_summary.R (needs 06a and 06b outputs).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/16_phase6_rmats_summary.R"
echo "[$(date)] 06c complete"
