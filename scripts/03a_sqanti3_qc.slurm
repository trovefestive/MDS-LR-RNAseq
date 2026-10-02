#!/usr/bin/env bash
#SBATCH --job-name=03a_sqanti_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --output=results/logs/03a_sqanti3_qc_%j.log

# =============================================================================
# 03a_sqanti3_qc.slurm — Phase 3: SQANTI3 QC of the merged ESPRESSO transcriptome
#   input : results/02_espresso/merged/merged_N2_R0_{updated.gtf,abundance.esp}
#   prep  : novel isoforms with gene_id "NA" get their ESPRESSO cluster as gene_id
#           (otherwise SQANTI3 pools them into one pseudo-gene for iso/gene ratios);
#           per-sample ESPRESSO counts (rounded) -> multi-sample FL count file (-fl)
#   QC    : GENCODE vM39 + GRCm39, refTSS CAGE peaks, PolyASite (lifted mm39), polyA motifs,
#           --include_ORF (needed for predicted_NMD in Phase 4). No short reads (none matched).
# Paper (v5.2.1) used short-read SJ/BAM support instead of CAGE/polyA; see phase0 report.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_SQANTI3"
IN=${RES}/02_espresso/merged
OUT=${RES}/03_sqanti; QC=${OUT}/qc
mkdir -p "$QC"
echo "[$(date)] $(sqanti3_qc.py --version 2>&1 | tail -1)"

# === input GTF with per-cluster gene_id for unassigned novel isoforms ===
GTF_IN=${OUT}/espresso_input.gtf
if [ ! -s "$GTF_IN" ]; then
    awk 'BEGIN{FS=OFS="\t"} /^#/{print; next}
         { if ($9 ~ /gene_id "NA";/) {
               match($9, /transcript_id "ESPRESSO:[^"]+"/)
               tid = substr($9, RSTART+15, RLENGTH-16); n = split(tid, p, ":")
               g = p[1] ":" p[2] ":" p[3]; sub(/gene_id "NA";/, "gene_id \"" g "\";", $9) }
           print }' "${IN}/merged_N2_R0_updated.gtf" > "${GTF_IN}.tmp"
    mv "${GTF_IN}.tmp" "$GTF_IN"
fi
echo "[$(date)] remaining gene_id NA lines: $(grep -c 'gene_id "NA"' "$GTF_IN" || true)"

# === FL count file: isoform + one integer column per sample ===
FL=${OUT}/espresso_fl_counts.tsv
if [ ! -s "$FL" ]; then
    awk 'BEGIN{FS=OFS="\t"} NR==1{printf "isoform"; for(i=4;i<=NF;i++) printf "\t%s",$i; print ""; next}
         {printf "%s",$1; for(i=4;i<=NF;i++) printf "\t%d",int($i+0.5); print ""}' \
        "${IN}/merged_N2_R0_abundance.esp" > "$FL"
fi
head -2 "$FL"

# === SQANTI3 QC ===
if [ ! -s "${QC}/espresso_classification.txt" ]; then
    sqanti3_qc.py \
        --isoforms "$GTF_IN" --refGTF "$GTF" --refFasta "$GENOME_FA" \
        --CAGE_peak "$CAGE_BED" --polyA_peak "$POLYA_BED" --polyA_motif_list "$POLYA_MOTIFS" \
        -fl "$FL" --include_ORF \
        -o espresso -d "$QC" -t "$SLURM_CPUS_PER_TASK" -n "$SLURM_CPUS_PER_TASK" --report skip
fi

ls -la "$QC" | head -30
echo "[$(date)] structural categories:"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="structural_category") c=i; next} {n[$c]++} END{for(k in n) print n[k]"\t"k}' \
    "${QC}/espresso_classification.txt" | sort -k1,1nr
echo "[$(date)] 03a complete"
