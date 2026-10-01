#!/usr/bin/env bash
#SBATCH --job-name=02c_isoquant
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=47:00:00
#SBATCH --output=results/logs/02c_isoquant_%j.log

# =============================================================================
# 02c_isoquant.slurm — Phase 2: IsoQuant transcript discovery + quantification (cross-check)
#   all 6 sorted BAMs jointly (one column per sample), GENCODE vM39
#   --data_type pacbio_ccs  HiFi reads
#   --fl_data               FLNC reads are full-length (both ends reliable)
#   --stranded forward      FLNC reads are oriented 5'->3' by isoseq refine
#   --polya_trimmed all     refine --require-polya trimmed polyA tails
# Resumable: if the output folder exists, rerun with --resume.
# --check_canonical removed: IsoQuant 4.0.0 crashes in that code path
#   (AttributeError: 'GeneInfo' object has no attribute 'reference_region', job 20679627).
#   SQANTI3 reports canonical/non-canonical junctions in Phase 3 instead.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_ISOQUANT"
OUT=${RES}/02_isoquant
echo "[$(date)] $(isoquant --version 2>&1)"

mapfile -t S < <(tail -n +2 "$SAMPLES" | cut -f1)
BAMS=(); for s in "${S[@]}"; do BAMS+=("${RES}/01_align/bam/${s}.sorted.bam"); done
echo "[$(date)] samples: ${S[*]}"

if [ -d "${OUT}/LR" ] && [ ! -e "${OUT}/.done" ]; then
    echo "[$(date)] resuming previous IsoQuant run"
    isoquant --resume --output "$OUT" --threads "$SLURM_CPUS_PER_TASK"
elif [ ! -e "${OUT}/.done" ]; then
    isoquant \
        --reference "$GENOME_FA" \
        --genedb "$GTF" --complete_genedb \
        --bam "${BAMS[@]}" --labels "${S[@]}" \
        --data_type pacbio_ccs --fl_data --stranded forward --polya_trimmed all \
        --sqanti_output \
        --prefix LR --output "$OUT" --threads "$SLURM_CPUS_PER_TASK"
fi
touch "${OUT}/.done"

ls -la "${OUT}/LR" | head -40
echo "[$(date)] 02c complete"
