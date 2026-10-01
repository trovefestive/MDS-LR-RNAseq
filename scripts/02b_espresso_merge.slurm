#!/usr/bin/env bash
#SBATCH --job-name=02b_espmerge
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=results/logs/02b_espresso_merge_%j.log

# =============================================================================
# 02b_espresso_merge.slurm — Phase 2: merge per-chromosome ESPRESSO_Q outputs
# (as in paper assemble_transcripts.sh "merging outputs"): one header + all rows.
# Output: results/02_espresso/merged/merged_N2_R0_{abundance.esp,updated.gtf}
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
ESP=${RES}/02_espresso; OUT=${ESP}/merged; mkdir -p "$OUT"
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)

# === check all chromosomes finished ===
for c in "${CHRS[@]}"; do
    ls "${ESP}/work/${c}/${c}_N2_R0_abundance.esp" "${ESP}/work/${c}/${c}_N2_R0_updated.gtf" >/dev/null
done

# === abundance ===
head -n 1 "${ESP}/work/chr1/chr1_N2_R0_abundance.esp" > "${OUT}/merged_N2_R0_abundance.esp"
for c in "${CHRS[@]}"; do tail -n +2 "${ESP}/work/${c}/${c}_N2_R0_abundance.esp"; done >> "${OUT}/merged_N2_R0_abundance.esp"

# === GTF (drop features with start > end, as in paper) ===
head -n 1 "${ESP}/work/chr1/chr1_N2_R0_updated.gtf" > "${OUT}/merged_N2_R0_updated.gtf"
for c in "${CHRS[@]}"; do tail -n +2 "${ESP}/work/${c}/${c}_N2_R0_updated.gtf"; done \
  | awk -F'\t' '!($4 > $5)' | sort -k1,1 -k4,4n >> "${OUT}/merged_N2_R0_updated.gtf"

# === summary ===
echo "[$(date)] header: $(head -1 "${OUT}/merged_N2_R0_abundance.esp" | cut -f1-9)"
echo "[$(date)] isoforms in abundance: $(tail -n +2 "${OUT}/merged_N2_R0_abundance.esp" | wc -l)"
echo "[$(date)] transcripts in GTF: $(awk -F'\t' '$3=="transcript"' "${OUT}/merged_N2_R0_updated.gtf" | wc -l)"
echo "[$(date)] novel (ESPRESSO) transcripts: $(awk -F'\t' '$3=="transcript"' "${OUT}/merged_N2_R0_updated.gtf" | grep -c 'ESPRESSO' || true)"
echo "[$(date)] per-sample total counts:"
awk -F'\t' 'NR==1{for(i=4;i<=NF;i++)h[i]=$i; n=NF; next} {for(i=4;i<=n;i++)s[i]+=$i} END{for(i=4;i<=n;i++) printf "  %s\t%.0f\n", h[i], s[i]}' "${OUT}/merged_N2_R0_abundance.esp"
echo "[$(date)] 02b complete"
