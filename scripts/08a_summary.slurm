#!/usr/bin/env bash
#SBATCH --job-name=08a_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=results/logs/08a_summary_%j.log

# =============================================================================
# 08a_summary.slurm — Phase 8: isoform-switch candidates with evidence levels (runs 18_phase8_candidates.R)
#   Needs Phase 3 (final set, IsoQuant comparison), Phase 6 (rMATS-long, DRIMSeq, edgeR) and Phase 7 (Salmon) outputs.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/18_phase8_candidates.R"
echo "[$(date)] 08a complete"
