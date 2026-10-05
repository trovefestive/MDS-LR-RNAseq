#!/usr/bin/env bash
#SBATCH --job-name=07c_rmats_prep
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-6
#SBATCH --output=results/logs/07c_sr_rmats_prep_%A_%a.log

# =============================================================================
# 07c_sr_rmats_prep.slurm — Phase 7: rMATS-turbo 4.1.2 "prep" per BAM (task N = row N of samples_shortread.tsv)
#   The single all-BAM run (--nthread 8) segfaulted after 33 s with no message, while every reduced test passed
#   (full GTF, 8 threads, unplaced contigs). Running prep separately per BAM is rMATS's documented route; it also
#   shows which BAM fails if one does. 07c_sr_rmats.slurm then runs "post" on the combined .rmats files.
#   Parameters must match the post step exactly (GTF, -t, --readLength, --variable-read-length, --libType, --novelSS).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_RMATS_SR"
S=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print $1}' "$SR_SAMPLES")
R=${RES}/07_shortread/rmats; P=${R}/prep/${S}; mkdir -p "$P"
BAM=${RES}/07_shortread/star/${S}.Aligned.sortedByCoord.out.bam
[ -s "$BAM" ] || { echo "missing $BAM"; exit 1; }
echo "$BAM" > "${P}/b1.txt"
echo "[$(date)] prep ${S}"

if [ ! -e "${P}/.done" ]; then
    rm -rf "${P}/tmp" "${P}/out"
    rmats.py --b1 "${P}/b1.txt" --gtf "${SR_REF}/expanded_final.gtf" \
        -t paired --readLength 150 --variable-read-length --libType "$SR_RMATS_LIBTYPE" --novelSS \
        --nthread 1 --od "${P}/out" --tmp "${P}/tmp" --task prep
    touch "${P}/.done"
fi
ls -la "${P}/tmp"
echo "[$(date)] prep ${S} complete"
