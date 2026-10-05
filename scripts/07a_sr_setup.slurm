#!/usr/bin/env bash
#SBATCH --job-name=07a_sr_setup
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/07a_sr_setup_%j.log

# =============================================================================
# 07a_sr_setup.slurm — Phase 7: references for the short-read cross-check (all GRCm39 + expanded transcriptome)
#   1) conda env lr-salmon (salmon) — skip if present; version recorded
#   2) expanded GTF = Phase 3 final set (171,224 isoforms; ESPRESSO GTF records) -> transcript FASTA (gffread)
#   3) Salmon decoy-aware index (gentrome = transcripts + GRCm39 genome as decoys), k=31
#   4) STAR GRCm39 index with the expanded GTF as splice-junction database (sjdbOverhang 149 for 2x150)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
module load miniforge; source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir -p "$SR_REF"; T=$SLURM_CPUS_PER_TASK

# === 1) salmon env ===
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_SALMON"; then
    conda create -y -n "$ENV_SALMON" -c conda-forge -c bioconda --strict-channel-priority salmon
fi
conda activate "$ENV_SALMON"; echo "[$(date)] $(salmon --version)"
grep -q "^## ${ENV_SALMON}" "${META}/software_versions.txt" || printf "## %s\nsalmon\t%s\n" "$ENV_SALMON" "$(salmon --version | awk '{print $2}')" >> "${META}/software_versions.txt"

# === 2) expanded GTF + transcript FASTA ===
GTF_X=${SR_REF}/expanded_final.gtf
if [ ! -s "$GTF_X" ]; then
    awk 'NR==FNR{k[$1]=1; next} /^#/{next}
         { if (match($0, /transcript_id "[^"]+"/)) { t=substr($0, RSTART+15, RLENGTH-16); if (t in k) print } }' \
        "${RES}/03_sqanti/final/final_isoform_ids.txt" "${RES}/03_sqanti/quant/R2_updated.gtf" > "$GTF_X"
fi
echo "[$(date)] expanded GTF transcripts: $(awk -F'\t' '$3=="transcript"' "$GTF_X" | wc -l)"
conda activate "$ENV_ALIGN"
[ -s "${SR_REF}/expanded_transcripts.fa" ] || gffread -w "${SR_REF}/expanded_transcripts.fa" -g "$GENOME_FA" "$GTF_X"
echo "[$(date)] transcript FASTA records: $(grep -c '^>' "${SR_REF}/expanded_transcripts.fa")"

# === 3) Salmon decoy-aware index ===
conda activate "$ENV_SALMON"
if [ ! -s "${SR_REF}/salmon_index/versionInfo.json" ]; then
    grep '^>' "$GENOME_FA" | cut -d' ' -f1 | sed 's/^>//' > "${SR_REF}/decoys.txt"
    cat "${SR_REF}/expanded_transcripts.fa" "$GENOME_FA" > "${TMPDIR_PROJ}/gentrome.fa"
    salmon index -t "${TMPDIR_PROJ}/gentrome.fa" -d "${SR_REF}/decoys.txt" -i "${SR_REF}/salmon_index" -k 31 -p "$T"
    rm -f "${TMPDIR_PROJ}/gentrome.fa"
fi

# === 4) STAR GRCm39 index with expanded GTF ===
conda activate "$ENV_RMATS_SR"; echo "[$(date)] STAR $(STAR --version)"
if [ ! -s "${SR_REF}/star_index/SA" ]; then
    mkdir -p "${SR_REF}/star_index"; rm -rf "${TMPDIR_PROJ}/star_genomegen_tmp"
    STAR --runMode genomeGenerate --runThreadN "$T" --genomeDir "${SR_REF}/star_index" \
         --genomeFastaFiles "$GENOME_FA" --sjdbGTFfile "$GTF_X" --sjdbOverhang 149 \
         --outTmpDir "${TMPDIR_PROJ}/star_genomegen_tmp" --limitGenomeGenerateRAM 90000000000
fi
ls -la "${SR_REF}" "${SR_REF}/star_index" | head -30
echo "[$(date)] 07a complete"
