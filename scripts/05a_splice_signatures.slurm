#!/usr/bin/env bash
#SBATCH --job-name=05a_splice
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=results/logs/05a_splice_signatures_%j.log

# =============================================================================
# 05a_splice_signatures.slurm — Phase 5: genotype-associated novel isoforms, junction classes,
# U2AF1 3'SS position tests (Q157 +1 G/A; S34 -3 C/T negative control), logos, ORA.
# Runs 14_phase5_splice_signatures.R. Needs Phase 3 final tables, SQANTI3 junctions, Phase 4 gene names.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/14_phase5_splice_signatures.R"
echo "[$(date)] 05a complete"
