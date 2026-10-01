#!/usr/bin/env bash
#SBATCH --job-name=00c_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=results/logs/00c_build_indexes_%j.log

# =============================================================================
# 00c_build_indexes.slurm — Phase 0: build indexes/derived references
#   samtools faidx + chrom sizes
#   minimap2 splice:hq index
#   GENCODE junction BED (paftools.js gff2bed) for minimap2 --junc-bed
#   PolyASite 2.0 GRCm38 -> mm39 liftOver (SQANTI3 --polyA_peak)
# Requires: 00a (lr-align env) and 00b (downloads). Skip-if-exists.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
activate_env "$ENV_ALIGN"
cd "$REF"
echo "[$(date)] minimap2 $(minimap2 --version); samtools $(samtools --version | head -1)"

# === faidx + chrom sizes ===
if [ ! -s "${GENOME_FA}.fai" ]; then
    echo "[$(date)] samtools faidx"
    samtools faidx "$GENOME_FA"
fi
cut -f1,2 "${GENOME_FA}.fai" > GRCm39.chrom.sizes

# === junction BED from GENCODE ===
if [ ! -s "$JUNC_BED" ]; then
    echo "[$(date)] paftools.js gff2bed -> $JUNC_BED"
    paftools.js gff2bed "$GTF" > "$JUNC_BED"
fi
echo "[$(date)] junction BED lines: $(wc -l < "$JUNC_BED")"

# === minimap2 index (splice:hq preset) ===
if [ ! -s "$MM2_INDEX" ]; then
    echo "[$(date)] building minimap2 index"
    minimap2 -x splice:hq -t "$SLURM_CPUS_PER_TASK" -d "$MM2_INDEX" "$GENOME_FA"
fi
echo "[$(date)] index: $(du -h "$MM2_INDEX" | cut -f1)"

# === PolyASite GRCm38 -> mm39 ===
if [ ! -s "$POLYA_BED" ]; then
    echo "[$(date)] lifting PolyASite atlas mm10 -> mm39"
    # Ensembl-style contigs (1, X, MT) -> UCSC (chr1, chrX, chrM); keep BED6
    awk 'BEGIN{OFS="\t"} !/^(track|#)/ {c=$1; if (c !~ /^chr/) c="chr"c; if (c=="chrMT") c="chrM";
         print c,$2,$3,$4,$5,$6}' atlas.clusters.2.0.GRCm38.96.bed > tmp_polyasite_mm10.bed
    liftOver tmp_polyasite_mm10.bed mm10ToMm39.over.chain.gz tmp_polyasite_mm39.bed polyasite_unlifted.bed
    # keep only contigs present in the GRCm39 primary assembly
    awk 'NR==FNR{ok[$1]=1; next} ($1 in ok)' GRCm39.chrom.sizes tmp_polyasite_mm39.bed \
        | sort -k1,1 -k2,2n > "$POLYA_BED"
    echo "[$(date)] PolyASite: input $(wc -l < tmp_polyasite_mm10.bed), lifted $(wc -l < "$POLYA_BED"), unlifted $(grep -vc '^#' polyasite_unlifted.bed || true)"
    rm -f tmp_polyasite_mm10.bed tmp_polyasite_mm39.bed
fi

# === refTSS contig check ===
echo "[$(date)] refTSS contigs not in genome: $(awk 'NR==FNR{ok[$1]=1; next} !($1 in ok){print $1}' GRCm39.chrom.sizes "$CAGE_BED" | sort -u | tr '\n' ' ')"

# === U2af1 sanity check (GRCm39 chr17) ===
echo "[$(date)] U2af1 in GTF:"
awk '$3=="gene" && /gene_name "U2af1";/' "$GTF" | cut -f1,4,5,7

echo "[$(date)] 00c complete"
