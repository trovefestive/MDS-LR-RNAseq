#!/usr/bin/env bash
#SBATCH --job-name=03c_requant
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-21
#SBATCH --output=results/logs/03c_requant_%A_%a.log

# =============================================================================
# 03c_requant_by_chr.slurm — Phase 3: requantify on the filtered GTF, then FSM-only (paper steps 4–5)
#   per chromosome (array task N -> chr1..chr19, chrX, chrY):
#   1) ESPRESSO_Q on <chr>_N2_R0_updated.filtered.gtf with --read_ratio_cutoff 2 (only GTF isoforms
#      quantified)                     -> work/<chr>/filtered/<chr>_N2_R2_{abundance.esp,updated.gtf}
#   2) keep only reads that are FSM to a filtered isoform (paper remove_nonfsm_data.pl; fixed copy in
#      scripts/paper_derived/)       -> work/<chr>/fsm/{0..5,bams,<chr>.tsv.updated}
#   3) ESPRESSO_Q again on the FSM-only inputs -> work/<chr>/fsm/<chr>_N2_R2_abundance.esp
# Both matrices are kept: R2 (all reads) and fsm (FSM reads only). Skip-if-exists per step.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_ESPRESSO"
ESPBIN=$(dirname "$(command -v ESPRESSO_S.pl)")   # .pl files have no shebang: call via perl

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
CHR=${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}
WD=${RES}/02_espresso/work/${CHR}
GTF_F=${WD}/${CHR}_N2_R0_updated.filtered.gtf
[ -s "$GTF_F" ] || { echo "missing $GTF_F (run 03b first)"; exit 1; }
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: ${CHR}; filtered GTF transcripts: $(awk -F'\t' '$3=="transcript"' "$GTF_F" | wc -l)"

# === 1) requant on filtered GTF (R2) ===
mkdir -p "${WD}/filtered"
if [ ! -s "${WD}/filtered/${CHR}_N2_R2_abundance.esp" ]; then
    echo "[$(date)] ESPRESSO_Q (filtered, R2)"
    perl "${ESPBIN}/ESPRESSO_Q.pl" -L "${WD}/${CHR}.tsv.updated" -O "${WD}/filtered" -A "$GTF_F" \
        --read_ratio_cutoff 2 -V "${WD}/filtered/${CHR}_compatible_isoform.tsv" -T 1
fi
echo "[$(date)] read classes vs filtered GTF:"
cut -f3 "${WD}/filtered/${CHR}_compatible_isoform.tsv" | sort | uniq -c | sort -k1,1nr

# === 2) FSM-read filter ===
if [ ! -s "${WD}/fsm/${CHR}.tsv.updated" ]; then
    echo "[$(date)] removing non-FSM reads"
    rm -rf "${WD}/fsm"          # incomplete previous attempt
    ok=0
    for attempt in 1 2 3; do    # samtools forks can fail transiently under strict overcommit
        if perl "${PROJ}/scripts/paper_derived/remove_nonfsm_data.pl" \
              "${WD}/filtered/${CHR}_compatible_isoform.tsv" "$WD" "${WD}/${CHR}.tsv.updated" \
              2> "${WD}/fsm_filter.log"; then ok=1; break; fi
        echo "[$(date)] FSM filter attempt ${attempt} failed:"; tail -5 "${WD}/fsm_filter.log"; rm -rf "${WD}/fsm"; sleep 60
    done
    [ "$ok" -eq 1 ] || exit 1
fi

# === 3) requant FSM-only ===
if [ ! -s "${WD}/fsm/${CHR}_N2_R2_abundance.esp" ]; then
    echo "[$(date)] ESPRESSO_Q (FSM reads only)"
    perl "${ESPBIN}/ESPRESSO_Q.pl" -L "${WD}/fsm/${CHR}.tsv.updated" -O "${WD}/fsm" -A "$GTF_F" \
        --read_ratio_cutoff 2 -V "${WD}/fsm/${CHR}_compatible_isoform.tsv" -T 1
fi

# === summary ===
for m in filtered fsm; do
    echo "[$(date)] ${m}: $(tail -n +2 "${WD}/${m}/${CHR}_N2_R2_abundance.esp" | wc -l) isoforms; per-sample totals:"
    awk -F'\t' 'NR==1{for(i=4;i<=NF;i++)h[i]=$i; n=NF; next} {for(i=4;i<=n;i++)s[i]+=$i}
                END{for(i=4;i<=n;i++) printf "  %s\t%.0f\n", h[i], s[i]}' "${WD}/${m}/${CHR}_N2_R2_abundance.esp"
done
echo "[$(date)] ${CHR} complete"
