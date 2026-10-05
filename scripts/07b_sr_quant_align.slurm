#!/usr/bin/env bash
#SBATCH --job-name=07b_sr_quant
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-6
#SBATCH --output=results/logs/07b_sr_quant_align_%A_%a.log

# =============================================================================
# 07b_sr_quant_align.slurm — Phase 7: per-sample short-read processing (task N = row N of samples_shortread.tsv)
#   1) Salmon quant against the expanded (long-read) transcriptome, decoy-aware; -l SR_SALMON_LIBTYPE (config.sh)
#   2) STAR 2-pass-free alignment to GRCm39 with the expanded-GTF junction database -> sorted BAM (for rMATS-turbo)
#   3) U2af1 Q157R genotype check at chr17:31867169 (T>C on + strand), as in Phase 1
# FASTQ directory: SR_FASTQ_DIR (untracked config.local.sh). Skip-if-exists per step.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
ROW=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "$SR_SAMPLES")
S=$(cut -f1 <<<"$ROW"); G=$(cut -f2 <<<"$ROW")
FQ1=${SR_FASTQ_DIR}/$(cut -f3 <<<"$ROW"); FQ2=${SR_FASTQ_DIR}/$(cut -f4 <<<"$ROW")
[ -s "$FQ1" ] && [ -s "$FQ2" ] || { echo "missing FASTQ for $S"; exit 1; }
O=${RES}/07_shortread; mkdir -p "${O}/salmon" "${O}/star" "${O}/genotype"
T=$SLURM_CPUS_PER_TASK
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${S} (${G})"

# === 1) Salmon ===
activate_env "$ENV_SALMON"
if [ ! -s "${O}/salmon/${S}/quant.sf" ]; then
    salmon quant -i "${SR_REF}/salmon_index" -l "$SR_SALMON_LIBTYPE" -1 "$FQ1" -2 "$FQ2" --validateMappings --gcBias --seqBias \
        -p "$T" -o "${O}/salmon/${S}"
fi
grep -E '"expected_format"|"compatible_fragment_ratio"' "${O}/salmon/${S}/lib_format_counts.json" || true
grep -E '"percent_mapped"|"num_processed"' "${O}/salmon/${S}/aux_info/meta_info.json" || true

# === 2) STAR to GRCm39 (expanded junction database) ===
activate_env "$ENV_RMATS_SR"
BAM=${O}/star/${S}.Aligned.sortedByCoord.out.bam
if [ ! -s "${BAM}" ]; then
    STAR --runThreadN "$T" --genomeDir "${SR_REF}/star_index" --readFilesIn "$FQ1" "$FQ2" --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM XS \
         --limitBAMsortRAM 30000000000 --outFileNamePrefix "${O}/star/${S}."
fi
grep -E "Uniquely mapped reads %|Number of input reads|% of reads mapped to multiple loci" "${O}/star/${S}.Log.final.out"

# === 3) genotype at U2af1 codon 157 ===
activate_env "$ENV_ALIGN"
[ -s "${BAM}.bai" ] || samtools index -@ 4 "$BAM"
{
  printf "sample\tgenotype_expected\tdepth\tref_T\talt_C\talt_C_frac\n"
  samtools mpileup -r chr17:31867169-31867169 -f "$GENOME_FA" -q 20 -Q 20 -d 0 \
      --no-output-ins --no-output-del --no-output-ends "$BAM" 2>/dev/null \
    | awk -v s="$S" -v g="$G" '{ b=$5; r=gsub(/[.,]/,"",b); c=gsub(/[Cc]/,"",b); o=gsub(/[AaGgTtNn]/,"",b); t=r+c+o
                                 printf "%s\t%s\t%d\t%d\t%d\t%.3f\n", s, g, t, r, c, (t>0?c/t:0) }'
} > "${O}/genotype/${S}.u2af1_Q157R.tsv"
cat "${O}/genotype/${S}.u2af1_Q157R.tsv"
echo "[$(date)] ${S} complete"
