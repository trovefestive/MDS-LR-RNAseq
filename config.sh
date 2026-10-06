#!/usr/bin/env bash
# =============================================================================
# config.sh — paths and parameters sourced by every script in scripts/
# Project: PacBio Kinnex long-read RNA-seq, U2af1 Q157R vs WT (mouse LK)
# =============================================================================

# === Project paths (PROJ = directory containing this file) ===
export PROJ=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export SAMPLES=${PROJ}/samples.tsv
export RAW=${PROJ}/data/raw                  # -> raw_files (read-only)
export REF=${PROJ}/data/ref
export RES=${PROJ}/results
export LOGS=${RES}/logs
export META=${RES}/metadata
export TMPDIR_PROJ=${PROJ}/tmp

# === Site-specific settings (untracked): SBATCH_ACCOUNT, SBATCH_PARTITION ===
[ -f "${PROJ}/config.local.sh" ] && source "${PROJ}/config.local.sh"

# === Reference: GRCm39 + GENCODE vM39 (keep one release throughout) ===
export GENCODE_REL=M39
export GENCODE_URL=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${GENCODE_REL}
export GENOME_FA=${REF}/GRCm39.primary_assembly.genome.fa
export GTF=${REF}/gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf
export JUNC_BED=${REF}/gencode.v${GENCODE_REL}.junctions.bed      # minimap2 --junc-bed
export MM2_INDEX=${REF}/GRCm39.primary_assembly.splice.mmi

# SQANTI3 orthogonal support
export CAGE_BED=${REF}/refTSS_v4.1_mouse_coordinate.mm39.noheader.bed   # header line stripped
export POLYA_BED=${REF}/atlas.clusters.2.0.GRCm38.96.lifted_mm39.bed   # PolyASite 2.0, lifted mm10 -> mm39
export POLYA_MOTIFS=${REF}/mouse_and_human.polyA_motif.txt

# === Conda envs (module load miniforge) ===
export ENV_ALIGN=lr-align          # minimap2, samtools, seqkit, pbmm2, nanoplot, bedtools, gffread, liftover
export ENV_ESPRESSO=lr-espresso
export ENV_ISOQUANT=lr-isoquant
export ENV_SQANTI3=lr-sqanti3
export ENV_RMATSLONG=lr-rmatslong
export ENV_R=lr-r

activate_env() {
    module load miniforge >/dev/null 2>&1
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$1"
}

# === Analysis parameters ===
export GROUP_REF=WT                # reference level: all effects are Q157R - WT
export GROUP_TEST=Q157R
export NOVEL_MIN_SAMPLES=2         # paper: novel isoform in >=2 samples
export NOVEL_MIN_READS=5           #        and >=5 total reads
export FDR_CUT=0.05
export DPROP_CUT=0.10

# === Phase 3: primary SQANTI3 filter (rules | ml_noORF | ml) ===
# rules chosen 2026-10-01: default ML filter learned coding/NMD (kept 337 novel); ML without ORF/NMD
# features dropped 63k novel isoforms with CAGE/polyA support as good as kept ones. See phase3 report.
export P3_FILTER=rules

# === Phase 7: short-read cross-check (path of the FASTQs: SR_FASTQ_DIR in untracked config.local.sh) ===
export SR_SAMPLES=${PROJ}/samples_shortread.tsv
export ENV_SALMON=lr-salmon                 # salmon (created by 07a)
export ENV_RMATS_SR=rmatsenv-2.7            # existing env: rMATS 4.1.2 + STAR 2.7.11b
export SR_REF=${PROJ}/data/ref/expanded     # transcriptome FASTA, Salmon index, STAR index (expanded GTF)
# Library type: all 6 libraries are dUTP-stranded (86–92% of Salmon fragments ISR), but Salmon's -l A auto-detect
# called 2/6 "IU" from its early-read sample, so the type is fixed for every sample. ISR = rMATS fr-firststrand.
export SR_SALMON_LIBTYPE=ISR
export SR_RMATS_LIBTYPE=fr-firststrand

# === Phase 7b: LeafCutter (annotation-free intron clustering; R package 0.2.x from github davidaknowles/leafcutter) ===
export ENV_LEAFCUTTER=lr-leafcutter                       # r-rstan 2.32 (Stan models use pre-2.33 array syntax) + regtools
export LEAFCUTTER_DIR=${PROJ}/data/tools/leafcutter       # repo clone (clustering + leafcutter_ds.R scripts)
