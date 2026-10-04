#!/usr/bin/env bash
#SBATCH --job-name=06a_rmatslong
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/logs/06a_rmats_long_%j.log

# =============================================================================
# 06a_rmats_long.slurm — Phase 6: rMATS-long differential isoform analysis (paper espresso_step6_rmats.sh)
#   input : Phase 3 final set — counts (R2 primary, duplicate-collapsed) + matching ESPRESSO GTF records
#   groups: group 1 = Q157R, group 2 = WT. rMATS-long computes
#           delta_isoform_proportion = group_1_avg - group_2_avg (calculate_isoform_proportion.py: calc_delta(a,b)=a-b)
#           so delta = Q157R - WT (rule 8); verified downstream on U2af1-213 (15_phase6_diff_isoforms.R).
#   --num-threads 1: with workers, cleanup_threads() waits only 1 s per process (join(1)) and raises
#     "thread(s) were not joined" if one exits slowly (job 20775824); single-threaded avoids the race.
#   thresholds: --adj-pvalue 0.05 (plan default; paper used 0.1), --delta-proportion 0.1 (as paper)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_RMATSLONG"
FIN=${RES}/03_sqanti/final; OUT=${RES}/06_diff/rmats_long; IN=${RES}/06_diff/inputs
mkdir -p "$OUT" "$IN"

# === inputs: abundance (ESPRESSO layout) + GTF restricted to the final set ===
zcat "${FIN}/final_counts.tsv.gz" > "${IN}/final_abundance.esp"
awk 'NR==FNR{k[$1]=1; next} /^#/{print; next}
     { if (match($0, /transcript_id "[^"]+"/)) { t=substr($0, RSTART+15, RLENGTH-16); if (t in k) print } }' \
    "${FIN}/final_isoform_ids.txt" "${RES}/03_sqanti/quant/R2_updated.gtf" > "${IN}/final_updated.gtf"
paste -sd, <(awk -F'\t' 'NR>1 && $2=="Q157R"{print $1}' "$SAMPLES") > "${IN}/group_Q157R.txt"
paste -sd, <(awk -F'\t' 'NR>1 && $2=="WT"{print $1}' "$SAMPLES")    > "${IN}/group_WT.txt"
echo "[$(date)] group 1 (Q157R): $(cat "${IN}/group_Q157R.txt")  group 2 (WT): $(cat "${IN}/group_WT.txt")"
echo "[$(date)] abundance rows: $(($(wc -l < "${IN}/final_abundance.esp") - 1)); GTF transcripts: $(awk -F'\t' '$3=="transcript"' "${IN}/final_updated.gtf" | wc -l)"

# === rMATS-long ===
if [ ! -e "${OUT}/.done" ]; then
    rmats-long rmats_long.py \
        --abundance "${IN}/final_abundance.esp" --updated-gtf "${IN}/final_updated.gtf" --gencode-gtf "$GTF" \
        --group-1 "${IN}/group_Q157R.txt" --group-2 "${IN}/group_WT.txt" \
        --group-1-name Q157R --group-2-name WT \
        --out-dir "$OUT" --num-threads 1 \
        --adj-pvalue 0.05 --delta-proportion 0.1 --process-top-n 100 --plot-file-type all
    touch "${OUT}/.done"
fi
ls -la "$OUT" | head -30
echo "[$(date)] 06a complete"
