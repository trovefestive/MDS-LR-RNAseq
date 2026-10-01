#!/usr/bin/env bash
# =============================================================================
# submit.sh — submit a pipeline script from the project root with site settings.
#   usage: scripts/submit.sh [sbatch options] scripts/NN_step.slurm
# Loads config.sh (+ untracked config.local.sh -> SBATCH_ACCOUNT / SBATCH_PARTITION),
# so the tracked scripts carry no account name or absolute path.
# Scripts must be submitted from the project root: they source ${SLURM_SUBMIT_DIR}/config.sh
# and write logs to results/logs/ relative to it.
# =============================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source ./config.sh
: "${SBATCH_ACCOUNT:?set SBATCH_ACCOUNT in config.local.sh}"
mkdir -p results/logs
sbatch "$@"
