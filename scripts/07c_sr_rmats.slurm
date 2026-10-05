#!/usr/bin/env bash
#SBATCH --job-name=07c_sr_rmats
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=results/logs/07c_sr_rmats_%j.log

# =============================================================================
# 07c_sr_rmats.slurm — Phase 7: rMATS-turbo 4.1.2 "post" on the short-read GRCm39 BAMs with the expanded GTF
#   Inputs: per-BAM .rmats files from 07c_sr_rmats_prep.slurm (same parameters), copied with unique prefixes.
#   --b1 = Q157R, --b2 = WT. rMATS_P/inclusion_level.py: IncLevelDifference = mean(PSI b1) - mean(PSI b2)
#   -> ΔPSI = Q157R − WT (rule 8). --libType fr-firststrand (dUTP; Salmon ISR in all 6 samples, checked below).
#   --novelSS: also detect events at splice sites absent from the GTF.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_RMATS_SR"
O=${RES}/07_shortread; R=${O}/rmats; mkdir -p "$R"

# === library type check: every sample must be mostly ISR (fraction of strand-informative fragments) ===
for s in $(tail -n +2 "$SR_SAMPLES" | cut -f1); do
    f=${O}/salmon/${s}/lib_format_counts.json
    isr=$(grep -o '"ISR": [0-9]*' "$f" | awk '{print $2}'); isf=$(grep -o '"ISF": [0-9]*' "$f" | awk '{print $2}')
    awk -v s="$s" -v a="$isr" -v b="$isf" 'BEGIN{fr=a/(a+b); printf "%s ISR fraction %.3f\n", s, fr; if (fr < 0.8) exit 1}' \
        || { echo "library of $s is not fr-firststrand"; exit 1; }
done

# === group files: b1 = Q157R, b2 = WT ===
bams() { awk -F'\t' -v g="$1" 'NR>1 && $2==g{print $1}' "$SR_SAMPLES" | sed "s|^|${O}/star/|; s|$|.Aligned.sortedByCoord.out.bam|" | paste -sd,; }
bams Q157R > "${R}/b1_Q157R.txt"; bams WT > "${R}/b2_WT.txt"

if [ ! -e "${R}/.done" ]; then
    rm -rf "${R}/tmp" "${R}/out"; mkdir -p "${R}/tmp"
    for s in $(tail -n +2 "$SR_SAMPLES" | cut -f1); do
        [ -e "${R}/prep/${s}/.done" ] || { echo "prep missing for $s"; exit 1; }
        python "$(dirname "$(which rmats.py)")/../rMATS/cp_with_prefix.py" "${s}_" "${R}/tmp" "${R}/prep/${s}"/tmp/*.rmats
    done
    rmats.py --b1 "${R}/b1_Q157R.txt" --b2 "${R}/b2_WT.txt" --gtf "${SR_REF}/expanded_final.gtf" \
        -t paired --readLength 150 --variable-read-length --libType "$SR_RMATS_LIBTYPE" --novelSS \
        --nthread "$SLURM_CPUS_PER_TASK" --od "${R}/out" --tmp "${R}/tmp" --task post
    touch "${R}/.done"
fi
ls "${R}/out" | head -40
for e in SE A3SS A5SS RI MXE; do
    awk -F'\t' -v e="$e" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {n++; if ($h["FDR"]<0.05 && ($h["IncLevelDifference"]>=0.1 || $h["IncLevelDifference"]<=-0.1)) s++}
        END{printf "%s JCEC: %d events, %d FDR<0.05 & |dPSI|>=0.1\n", e, n, s}' "${R}/out/${e}.MATS.JCEC.txt"
done
echo "[$(date)] 07c complete"
