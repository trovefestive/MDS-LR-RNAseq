#!/usr/bin/env bash
#SBATCH --job-name=03b2_ml_noORF
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/03b2_sqanti3_ml_noORF_%j.log

# =============================================================================
# 03b2_sqanti3_ml_noORF.slurm — Phase 3: ML filter rerun without ORF/NMD and per-sample count features
# Why: with default features the random forest's top predictors were predicted_NMD and coding
#   (TP set = reference FSM, mostly coding) and it called 99.7% of novel isoforms artifacts.
#   Filtering on coding/NMD would bias Phase 4 (NMD analysis); per-sample FL columns let it learn
#   sample/genotype-specific patterns. Excluded columns: scripts/03_ml_remove_columns.txt.
#   Total FL (summed reads) is kept as the coverage feature.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_SQANTI3"
OUT=${RES}/03_sqanti
sqanti3_filter.py ml --sqanti_class "${OUT}/qc/espresso_classification.txt" \
    --filter_gtf "${OUT}/qc/espresso_corrected.gtf" -r "${PROJ}/scripts/03_ml_remove_columns.txt" \
    -o espresso -d "${OUT}/filter_ml_noORF" -c "$SLURM_CPUS_PER_TASK" --skip_report
echo "[$(date)] variable importance (top 12):"; head -12 "${OUT}/filter_ml_noORF/classifier_variable-importance_table.txt"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="isoform")a=i; if($i=="filter_result")f=i}; next}
            {k=($a ~ /^ENSMUST/ ? "gencode" : "novel"); n[k"\t"$f]++} END{for(x in n) print x"\t"n[x]}' \
    "${OUT}/filter_ml_noORF/espresso_ML_classification.txt" | sort
echo "[$(date)] 03b2 complete"
