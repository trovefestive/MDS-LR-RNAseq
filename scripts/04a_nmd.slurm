#!/usr/bin/env bash
#SBATCH --job-name=04a_nmd
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=results/logs/04a_nmd_%j.log

# =============================================================================
# 04a_nmd.slurm — Phase 4: ORF / NMD characterization (runs 13_phase4_nmd.R)
#   inputs : results/03_sqanti/final/{final_counts,final_isoform_annotation}.tsv.gz (R2 primary matrix)
#   prep   : GENCODE vM39 gene_id -> gene_name table (results/04_nmd/gene_names.tsv)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
OUT=${RES}/04_nmd; mkdir -p "$OUT"

# === gene names from GENCODE (gene lines only) ===
if [ ! -s "${OUT}/gene_names.tsv" ]; then
    awk -F'\t' 'BEGIN{print "gene_id\tgene_name\tgene_type"} $3=="gene" {
        match($9, /gene_id "[^"]+"/);   id = substr($9, RSTART+9,  RLENGTH-10)
        match($9, /gene_name "[^"]+"/); nm = substr($9, RSTART+11, RLENGTH-12)
        match($9, /gene_type "[^"]+"/); ty = substr($9, RSTART+11, RLENGTH-12)
        print id "\t" nm "\t" ty }' "$GTF" > "${OUT}/gene_names.tsv"
fi
echo "[$(date)] gene names: $(($(wc -l < "${OUT}/gene_names.tsv") - 1))"

activate_env "$ENV_R"
Rscript "${PROJ}/scripts/13_phase4_nmd.R"
echo "[$(date)] 04a complete"
