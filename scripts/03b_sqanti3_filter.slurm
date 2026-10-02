#!/usr/bin/env bash
#SBATCH --job-name=03b_sqanti_filt
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=results/logs/03b_sqanti3_filter_%j.log

# =============================================================================
# 03b_sqanti3_filter.slurm — Phase 3: SQANTI3 filtering + GENCODE restore + per-chr GTFs
#   1) ML filter (default settings, as in the paper) and 2) rules filter (default JSON);
#      the primary one is P3_FILTER in config.sh (rules since 2026-10-01; see phase3 report)
#   3) keep = primary-filter "Isoform" + GENCODE (ENSMUST) transcripts flagged "Artifact" (restored, as the
#      paper restored ENST: reference treated as truth set)
#   4) subset the ESPRESSO GTF (original ESPRESSO formatting, needed by ESPRESSO_Q) to keep-list
#      and split per chromosome -> results/02_espresso/work/<chr>/<chr>_N2_R0_updated.filtered.gtf
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_SQANTI3"
OUT=${RES}/03_sqanti; QC=${OUT}/qc
CLASS=${QC}/espresso_classification.txt
[ -s "$CLASS" ] || { echo "missing $CLASS (run 03a first)"; exit 1; }
CORR_GTF=$(ls "${QC}"/espresso_corrected.gtf)

# === 1) ML filter ===
if [ ! -s "${OUT}/filter_ml/espresso_ML_classification.txt" ]; then
    sqanti3_filter.py ml --sqanti_class "$CLASS" --filter_gtf "$CORR_GTF" \
        -o espresso -d "${OUT}/filter_ml" -c "$SLURM_CPUS_PER_TASK" --skip_report
fi
# === 2) rules filter (comparison) ===
if [ ! -s "${OUT}/filter_rules/espresso_RulesFilter_classification.txt" ]; then
    sqanti3_filter.py rules --sqanti_class "$CLASS" --filter_gtf "$CORR_GTF" \
        -o espresso -d "${OUT}/filter_rules" -c "$SLURM_CPUS_PER_TASK" --skip_report
fi

# === 3) keep list: primary filter (P3_FILTER in config.sh) pass + restored GENCODE artifacts ===
case "$P3_FILTER" in
    rules)    ML=${OUT}/filter_rules/espresso_RulesFilter_classification.txt ;;
    ml_noORF) ML=${OUT}/filter_ml_noORF/espresso_ML_classification.txt ;;
    ml)       ML=${OUT}/filter_ml/espresso_ML_classification.txt ;;
    *) echo "unknown P3_FILTER=$P3_FILTER"; exit 1 ;;
esac
[ -s "$ML" ] || { echo "missing $ML"; exit 1; }
echo "[$(date)] primary filter: ${P3_FILTER} (${ML})"
ln -sfn "$ML" "${OUT}/primary_filter_classification.txt"
KEEP=${OUT}/keep_isoforms.txt
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="isoform")a=i; if($i=="filter_result")f=i}; next}
            $f=="Isoform" || ($f=="Artifact" && $a ~ /^ENSMUST/) {print $a}' "$ML" | sort -u > "$KEEP"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="isoform")a=i; if($i=="filter_result")f=i}; next}
            {k=($a ~ /^ENSMUST/ ? "gencode" : "novel"); n[k"\t"$f]++} END{for(x in n) print x"\t"n[x]}' "$ML" \
  | sort > "${OUT}/filter_${P3_FILTER}_counts.tsv"
echo "[$(date)] ${P3_FILTER} filter result (origin, result, n):"; cat "${OUT}/filter_${P3_FILTER}_counts.tsv"
echo "[$(date)] keep list (${P3_FILTER} pass + restored GENCODE): $(wc -l < "$KEEP")"

# === 4) filtered ESPRESSO GTF + per-chromosome split ===
FILT=${OUT}/espresso_filtered.ensrestore.gtf
# gene_id "NA" (unassigned novel isoforms) -> ESPRESSO cluster ID, as for the SQANTI3 input (03a): ESPRESSO_Q
# otherwise pools all NA-gene isoforms of a chromosome into one pseudo-gene and reports some twice (Phase 3 v1).
awk 'BEGIN{FS=OFS="\t"} NR==FNR{k[$1]=1; next} /^#/{print; next}
     { if (match($9, /transcript_id "[^"]+"/)) { t=substr($9, RSTART+15, RLENGTH-16); if (!(t in k)) next
         if ($9 ~ /gene_id "NA";/) { n=split(t, p, ":"); g=p[1] ":" p[2] ":" p[3]
                                     sub(/gene_id "NA";/, "gene_id \"" g "\";", $9) }
         print } }' \
    "$KEEP" "${RES}/02_espresso/merged/merged_N2_R0_updated.gtf" > "${FILT}.tmp"
mv "${FILT}.tmp" "$FILT"
echo "[$(date)] filtered GTF transcripts: $(awk -F'\t' '$3=="transcript"' "$FILT" | wc -l); gene_id NA lines left: $(grep -c 'gene_id "NA"' "$FILT" || true)"
for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY; do
    awk -v c="$c" -F'\t' '$1==c' "$FILT" > "${RES}/02_espresso/work/${c}/${c}_N2_R0_updated.filtered.gtf"
done
wc -l "${RES}"/02_espresso/work/*/*_N2_R0_updated.filtered.gtf | tail -1
echo "[$(date)] 03b complete"
