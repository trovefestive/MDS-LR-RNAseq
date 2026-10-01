#!/usr/bin/env bash
#SBATCH --job-name=01a_align
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --array=1-6
#SBATCH --output=results/logs/01a_qc_align_%A_%a.log

# =============================================================================
# 01a_qc_align.slurm — Phase 1: per-sample read QC, alignment, genotype check
#   array task N = data row N of samples.tsv
#   1) seqkit stats + read-length histogram (100 bp bins) of FLNC fastq
#   2) minimap2 -ax splice:hq -uf --secondary=no --junc-bed (GENCODE vM39)
#      -> coordinate-sorted, indexed BAM with @RG (ESPRESSO needs sorted BAM + NM tag)
#   3) flagstat / samtools stats
#   4) U2af1 Q157R pileup at chr17:31867169 (GRCm39; minus-strand gene,
#      CAG>CGG on transcript = T>C on genome + strand)
# Skip-if-exists on each output.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_ALIGN"

# === sample for this task ===
ROW=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$SAMPLES")
S=$(echo "$ROW" | cut -f1)
GENO=$(echo "$ROW" | cut -f2)
FQ=${PROJ}/$(echo "$ROW" | cut -f5)   # samples.tsv paths are relative to PROJ
[ -s "$FQ" ] || { echo "missing fastq for $S: $FQ"; exit 1; }
T=${SLURM_CPUS_PER_TASK}

QC=${RES}/01_align/qc; BAMDIR=${RES}/01_align/bam; GT=${RES}/01_align/genotype
mkdir -p "$QC" "$BAMDIR" "$GT"
BAM=${BAMDIR}/${S}.sorted.bam
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${S} (${GENO}) fastq=${FQ}"
echo "[$(date)] minimap2 $(minimap2 --version); samtools $(samtools --version | head -1 | cut -d' ' -f2); seqkit $(seqkit version | cut -d' ' -f2)"

# === 1) read QC ===
if [ ! -s "${QC}/${S}.seqkit_stats.tsv" ]; then
    echo "[$(date)] seqkit stats"
    seqkit stats -a -T -j 4 "$FQ" | sed "2s/^[^\t]*/${S}/" > "${QC}/${S}.seqkit_stats.tsv"
fi
if [ ! -s "${QC}/${S}.length_hist.tsv" ]; then
    echo "[$(date)] length histogram"
    # header written separately so sort cannot move it (bug fixed 2026-09-30)
    {
      printf "sample\tbin_start\tn\n"
      seqkit fx2tab -n -l -j 4 "$FQ" \
        | awk -v s="$S" 'BEGIN{OFS="\t"} {b=int($NF/100)*100; h[b]++} END{for (k in h) print s,k,h[k]}' \
        | sort -k2,2n
    } > "${QC}/${S}.length_hist.tsv"
fi

# === 2) alignment ===
if [ ! -s "${BAM}.bai" ]; then
    echo "[$(date)] minimap2 splice:hq -> ${BAM}"
    minimap2 -t "$T" -ax splice:hq -uf --secondary=no \
        --junc-bed "$JUNC_BED" \
        -R "@RG\tID:${S}\tSM:${S}\tPL:PACBIO" \
        "$MM2_INDEX" "$FQ" \
      | samtools sort -@ 4 -m 6G -T "${TMPDIR_PROJ}/${S}_sort" -o "${BAM}.tmp.bam" -
    mv "${BAM}.tmp.bam" "$BAM"
    samtools index -@ 4 "$BAM"
    echo "[$(date)] alignment done: $(du -h "$BAM" | cut -f1)"
fi

# === 3) alignment stats ===
if [ ! -s "${QC}/${S}.flagstat.txt" ]; then
    samtools flagstat -@ 4 "$BAM" > "${QC}/${S}.flagstat.txt"
fi
if [ ! -s "${QC}/${S}.samtools_stats_SN.tsv" ]; then
    samtools stats -@ 4 "$BAM" | grep '^SN' | cut -f2- > "${QC}/${S}.samtools_stats_SN.tsv"
fi
cat "${QC}/${S}.flagstat.txt"

# === 4) U2af1 Q157R genotype (chr17:31867169, ref T, alt C) ===
POS=chr17:31867169-31867169
{
  printf "sample\tgenotype_expected\tpos\tref\tdepth\tref_T\talt_C\tother\tdel\talt_C_frac\n"
  samtools mpileup -r "$POS" -f "$GENOME_FA" -q 20 -Q 0 -d 0 \
      --no-output-ins --no-output-del --no-output-ends "$BAM" 2>/dev/null \
    | awk -v s="$S" -v g="$GENO" '
        { b=$5; nref=gsub(/[.,]/,"",b); nC=gsub(/[Cc]/,"",b); nOth=gsub(/[AaGgTtNn]/,"",b);
          nDel=gsub(/[*#]/,"",b); tot=nref+nC+nOth;
          printf "%s\t%s\t%s:%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\n", s,g,$1,$2,$3,tot,nref,nC,nOth,nDel,(tot>0?nC/tot:0) }'
} > "${GT}/${S}.u2af1_Q157R.tsv"
echo "[$(date)] U2af1 Q157R pileup:"; cat "${GT}/${S}.u2af1_Q157R.tsv"

echo "[$(date)] ${S} complete"
