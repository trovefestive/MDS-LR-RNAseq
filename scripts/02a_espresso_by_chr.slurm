#!/usr/bin/env bash
#SBATCH --job-name=02a_espresso
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=47:00:00
#SBATCH --array=1-21
#SBATCH --output=results/logs/02a_espresso_%A_%a.log

# =============================================================================
# 02a_espresso_by_chr.slurm — Phase 2: ESPRESSO assembly + quantification, per chromosome
# Follows Miller et al. 2026 assemble_transcripts.sh (per-chr split, all samples jointly):
#   0) split each sample BAM to this chromosome, primary alignments only (-F 0x900)
#   1) ESPRESSO_S  (all 6 samples, GENCODE vM39 annotation for this chr)
#   2) ESPRESSO_C  (each sample; sequential, -T 8, up to 3 attempts)
#   3) ESPRESSO_Q  (defaults N=2, R=0 as in the paper) -> *_abundance.esp, *_updated.gtf,
#                   compatible-isoform tsv (needed for the later FSM-read filter)
# Array task N -> chromosome N of: chr1..chr19, chrX, chrY (chrM and scaffolds excluded, as in paper)
# Skip-if-exists per step (a timed-out task can be resubmitted and resumes).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_ESPRESSO"
ESPBIN=$(dirname "$(command -v ESPRESSO_S.pl)")   # .pl files have no shebang: call via perl

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
CHR=${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}
T=${SLURM_CPUS_PER_TASK}
C_THREADS=8
ESP=${RES}/02_espresso
IN=${ESP}/inputs/${CHR}
WD=${ESP}/work/${CHR}
mkdir -p "$IN" "$WD"
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${CHR} on $(hostname)"
echo "[$(date)] ESPRESSO $(grep -o "[0-9]\.[0-9]\.[0-9]" "${ESPBIN}/ESPRESSO_Version.pm" | head -1); samtools $(samtools --version | head -1 | cut -d' ' -f2)"

# === 0) per-chromosome inputs ===
GTF_CHR=${IN}/${CHR}.gencode.v${GENCODE_REL}.gtf
[ -s "$GTF_CHR" ] || awk -v c="$CHR" '$1==c' "$GTF" > "$GTF_CHR"
echo "[$(date)] ${CHR} GTF transcripts: $(awk '$3=="transcript"' "$GTF_CHR" | wc -l)"

TSV=${WD}/${CHR}.tsv
if [ ! -s "$TSV" ]; then
    : > "${TSV}.tmp"
    tail -n +2 "$SAMPLES" | cut -f1 | while read -r S; do
        OUT=${IN}/${S}.${CHR}.bam
        samtools view -@ 4 -F 0x900 -b -o "$OUT" "${RES}/01_align/bam/${S}.sorted.bam" "$CHR"
        samtools index "$OUT"
        printf "%s\t%s\n" "$OUT" "$S" >> "${TSV}.tmp"
        echo "[$(date)] ${S} ${CHR}: $(samtools view -c "$OUT") primary alignments"
    done
    mv "${TSV}.tmp" "$TSV"
fi
cat "$TSV"

# === 1) ESPRESSO_S ===
if [ ! -s "${WD}/${CHR}.tsv.updated" ]; then
    echo "[$(date)] ESPRESSO_S"
    perl "${ESPBIN}/ESPRESSO_S.pl" -L "$TSV" -F "$GENOME_FA" -A "$GTF_CHR" -O "$WD" -T 6
fi
echo "[$(date)] S done; tsv.updated:"; cat "${WD}/${CHR}.tsv.updated"

# === 2) ESPRESSO_C (per sample, sequential) ===
for ID in $(cut -f3 "${WD}/${CHR}.tsv.updated" | sort -u); do
    if [ -e "${WD}/.C_done_${ID}" ]; then echo "[$(date)] C ${ID} already done"; continue; fi
    echo "[$(date)] ESPRESSO_C sample ${ID}"
    # C is very verbose -> per-sample log, gzipped on success. Default --sort_buffer_size (2G).
    # Compute nodes use strict overcommit (vm.overcommit_memory=2): each worker thread's
    # fork() for `sort` must commit the whole perl process, and fails transiently (system()=-1)
    # when several tasks share a node (20679746_3/_15, 20680397_14). Use -T 8 (as in the paper)
    # and retry up to 3 times.
    CLOG=${WD}/C_${ID}.log
    ok=0
    for attempt in 1 2 3; do
        if perl "${ESPBIN}/ESPRESSO_C.pl" -I "$WD" -F "$GENOME_FA" -X "$ID" -T "$C_THREADS" > "$CLOG" 2>&1; then
            ok=1; break
        fi
        echo "[$(date)] ESPRESSO_C sample ${ID} attempt ${attempt} FAILED:"
        grep -m3 -E "abnormally|[Ee]rror|[Dd]ied" "$CLOG" | cut -c1-300 || true
        sleep 60
    done
    if [ "$ok" -ne 1 ]; then echo "[$(date)] ESPRESSO_C sample ${ID} failed 3 times"; tail -20 "$CLOG"; exit 1; fi
    gzip -f "$CLOG"
    touch "${WD}/.C_done_${ID}"
done

# === 3) ESPRESSO_Q ===
if ! ls "${WD}"/*_abundance.esp >/dev/null 2>&1; then
    echo "[$(date)] ESPRESSO_Q"
    perl "${ESPBIN}/ESPRESSO_Q.pl" -L "${WD}/${CHR}.tsv.updated" -A "$GTF_CHR" \
        -V "${WD}/${CHR}_compatible_isoform.tsv" -T 1
fi
ls -la "$WD" | grep -E "abundance|updated.gtf|compatible"
echo "[$(date)] ${CHR}: $(tail -n +2 "${WD}"/*_abundance.esp | wc -l) isoforms quantified;" \
     "$(grep -c -P '\ttranscript\t' "${WD}"/*_updated.gtf) in updated GTF"
echo "[$(date)] ${CHR} complete"
