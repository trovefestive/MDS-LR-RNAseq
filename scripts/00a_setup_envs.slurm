#!/usr/bin/env bash
#SBATCH --job-name=00a_envs
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=results/logs/00a_setup_envs_%j.log

# =============================================================================
# 00a_setup_envs.slurm — Phase 0: create conda envs for the long-read pipeline
# Skip-if-exists: an env that already exists is not rebuilt.
# Each env is built independently; failures are reported at the end.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
module load miniforge
source "$(conda info --base)/etc/profile.d/conda.sh"
echo "[$(date)] conda: $(conda --version)"

CH="-c conda-forge -c bioconda --strict-channel-priority"
FAILED=()

# === helper ===
make_env() {
    local name=$1; shift
    if conda env list | awk '{print $1}' | grep -qx "$name"; then
        echo "[$(date)] $name exists — skipping"
        return 0
    fi
    echo "[$(date)] creating $name: $*"
    if conda create -y -n "$name" $CH "$@"; then
        echo "[$(date)] $name OK"
    else
        echo "[$(date)] $name FAILED"
        FAILED+=("$name")
    fi
}

# === envs ===
make_env "$ENV_ALIGN"     minimap2=2.31 samtools seqkit pbmm2 nanoplot bedtools gffread ucsc-liftover k8
make_env "$ENV_ESPRESSO"  espresso=1.6.0 samtools
make_env "$ENV_ISOQUANT"  isoquant=4.0.0
make_env "$ENV_SQANTI3"   sqanti3=6.0.2 "setuptools<81"   # TD2 (ORF calls) imports pkg_resources, removed in setuptools>=81
make_env "$ENV_RMATSLONG" rmats-long=2.1.0
make_env "$ENV_R"         r-base=4.4 r-tidyverse r-data.table r-ggrepel r-ggseqlogo r-patchwork \
                          bioconductor-edger bioconductor-limma bioconductor-deseq2 bioconductor-drimseq \
                          bioconductor-dexseq bioconductor-isoformswitchanalyzer bioconductor-clusterprofiler \
                          bioconductor-org.mm.eg.db bioconductor-bsgenome.mmusculus.ucsc.mm39 \
                          bioconductor-rtracklayer bioconductor-genomicfeatures

# === record versions ===
mkdir -p "$META"
OUT=${META}/software_versions.txt
{
    echo "# software versions — written $(date)"
    for e in "$ENV_ALIGN" "$ENV_ESPRESSO" "$ENV_ISOQUANT" "$ENV_SQANTI3" "$ENV_RMATSLONG" "$ENV_R"; do
        echo "## $e"
        conda list -n "$e" 2>/dev/null | grep -E -i \
          "^(minimap2|samtools|seqkit|pbmm2|nanoplot|bedtools|gffread|ucsc-liftover|espresso|isoquant|sqanti3|rmats-long|r-base|bioconductor-(edger|limma|deseq2|drimseq|dexseq|isoformswitchanalyzer|clusterprofiler)) " \
          | awk '{print $1"\t"$2}' || echo "(env missing)"
    done
} > "$OUT"
cat "$OUT"

if [ ${#FAILED[@]} -gt 0 ]; then
    echo "[$(date)] FAILED envs: ${FAILED[*]}"
    exit 1
fi
echo "[$(date)] all envs ready"
