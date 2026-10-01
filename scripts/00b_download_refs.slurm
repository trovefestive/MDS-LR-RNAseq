#!/usr/bin/env bash
#SBATCH --job-name=00b_refs
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=results/logs/00b_download_refs_%j.log

# =============================================================================
# 00b_download_refs.slurm — Phase 0: download references + paper code
#   GRCm39 primary assembly + GENCODE vM39 primary_assembly GTF
#   refTSS v4.1 mouse CAGE peaks (mm39)
#   PolyASite 2.0 mouse atlas (GRCm38; lifted to mm39 in 00c) + UCSC mm10->mm39 chain
#   Miller et al. 2026 analysis code (zenodo 10.5281/zenodo.20314532)
# Skip-if-exists on every file. No conda env needed.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
mkdir -p "$REF" "$META" "${PROJ}/data/paper_code"
cd "$REF"

# === helper: download (+ optional gunzip) ===
fetch() {
    local url=$1 out=$2 unzip=${3:-no}
    local final=$out; [ "$unzip" = "gunzip" ] && final=${out%.gz}
    if [ -s "$final" ]; then echo "[$(date)] exists: $final"; return 0; fi
    echo "[$(date)] downloading $url"
    wget -q --tries=3 -O "$out" "$url"
    [ "$unzip" = "gunzip" ] && gunzip -f "$out"
    echo "[$(date)] done: $final ($(du -h "$final" | cut -f1))"
}

# === GENCODE vM39 ===
fetch ${GENCODE_URL}/GRCm39.primary_assembly.genome.fa.gz              GRCm39.primary_assembly.genome.fa.gz gunzip
fetch ${GENCODE_URL}/gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf.gz \
      gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf.gz gunzip

# === refTSS v4.1 mouse (mm39) ===
fetch https://reftss.riken.jp/datafiles/current/mouse/refTSS_v4.1_mouse_coordinate.mm39.bed.gz \
      refTSS_v4.1_mouse_coordinate.mm39.bed.gz gunzip

# === PolyASite 2.0 mouse (GRCm38) + liftOver chain ===
fetch https://polyasite.unibas.ch/download/atlas/2.0/GRCm38.96/atlas.clusters.2.0.GRCm38.96.bed.gz \
      atlas.clusters.2.0.GRCm38.96.bed.gz gunzip
fetch https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz \
      mm10ToMm39.over.chain.gz

# === Paper code (Miller et al. 2026) ===
cd "${PROJ}/data/paper_code"
fetch "https://zenodo.org/api/records/20314532/files/chrisamiller/longread-0.1.zip/content" longread-0.1.zip
[ -d longread-0.1 ] || [ -n "$(ls -d */ 2>/dev/null)" ] || unzip -q longread-0.1.zip

# === checksums / provenance ===
{
    echo "# reference downloads — $(date)"
    echo "GENCODE release: ${GENCODE_REL} (${GENCODE_URL})"
    cd "$REF"; md5sum *.fa *.gtf *.bed *.chain.gz 2>/dev/null
} > "${META}/reference_provenance.txt"
cat "${META}/reference_provenance.txt"

# === quick sanity ===
cd "$REF"
echo "[$(date)] genome contigs: $(grep -c '^>' GRCm39.primary_assembly.genome.fa)"
echo "[$(date)] GTF genes: $(awk '$3=="gene"' gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf | wc -l)"
echo "[$(date)] GTF transcripts: $(awk '$3=="transcript"' gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf | wc -l)"
echo "[$(date)] refTSS peaks: $(wc -l < refTSS_v4.1_mouse_coordinate.mm39.bed); first line:"; head -1 refTSS_v4.1_mouse_coordinate.mm39.bed
echo "[$(date)] PolyASite clusters: $(wc -l < atlas.clusters.2.0.GRCm38.96.bed); first line:"; head -1 atlas.clusters.2.0.GRCm38.96.bed
echo "[$(date)] 00b complete"
