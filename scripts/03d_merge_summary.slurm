#!/usr/bin/env bash
#SBATCH --job-name=03d_p3summary
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/03d_merge_summary_%j.log

# =============================================================================
# 03d_merge_summary.slurm — Phase 3: merge requantified matrices, compare, summarize
#   1) merge per-chr ESPRESSO outputs: filtered/ (R2, all reads) and fsm/ (FSM reads only)
#   2) gffcompare of the filtered transcriptome against IsoQuant models and the old LRP2 set
#      (intron-chain match = class code "=")
#   3) 12_phase3_summary.R: read loss, primary-matrix decision, novel-isoform rule, TMM CPM,
#      category/length/expression summaries, plots
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
W=${RES}/02_espresso/work; Q=${RES}/03_sqanti/quant; mkdir -p "$Q"
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)

# === 1) merge ===
for m in filtered fsm; do
    tag=$([ "$m" = filtered ] && echo R2 || echo fsm)
    for c in "${CHRS[@]}"; do ls "${W}/${c}/${m}/${c}_N2_R2_abundance.esp" >/dev/null; done
    head -n 1 "${W}/chr1/${m}/chr1_N2_R2_abundance.esp" > "${Q}/${tag}_abundance.esp"
    for c in "${CHRS[@]}"; do tail -n +2 "${W}/${c}/${m}/${c}_N2_R2_abundance.esp"; done >> "${Q}/${tag}_abundance.esp"
    echo "[$(date)] ${tag}: $(tail -n +2 "${Q}/${tag}_abundance.esp" | wc -l) isoforms"
done
for c in "${CHRS[@]}"; do tail -n +2 "${W}/${c}/filtered/${c}_N2_R2_updated.gtf"; done \
  | awk -F'\t' '!($4 > $5)' | sort -k1,1 -k4,4n > "${Q}/R2_updated.gtf"
echo "[$(date)] R2 GTF transcripts: $(awk -F'\t' '$3=="transcript"' "${Q}/R2_updated.gtf" | wc -l)"

# === 2) gffcompare vs IsoQuant and old LRP2 ===
module load gffcompare >/dev/null 2>&1 || true
command -v gffcompare >/dev/null || { activate_env "$ENV_ALIGN"; }
CMP=${RES}/03_sqanti/compare; mkdir -p "$CMP"
gffcompare -r "${RES}/02_isoquant/LR/LR.transcript_models.gtf" -o "${CMP}/vs_isoquant" "${Q}/R2_updated.gtf"
OLD=${PROJ}/2026_04/260327_MDS_lrp2/S2_TRANSCRIPTOME/M3_FILTER_TRANSCRIPTOME/MDS.transcriptome.filtered.gtf
if [ -s "$OLD" ]; then
    gffcompare -r "$OLD" -o "${CMP}/vs_lrp2" "${Q}/R2_updated.gtf"
fi
mv -f "${Q}"/vs_*.tmap "${Q}"/vs_*.refmap "$CMP"/ 2>/dev/null || true
ls -la "$CMP"

# === 3) summary ===
activate_env "$ENV_R"
Rscript "${PROJ}/scripts/12_phase3_summary.R"
echo "[$(date)] 03d complete"
