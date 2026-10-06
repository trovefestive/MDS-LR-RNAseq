#!/usr/bin/env bash
#SBATCH --job-name=07e_leafcutter
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=results/logs/07e_leafcutter_%j.log

# =============================================================================
# 07e_leafcutter.slurm — Phase 7b: LeafCutter differential intron usage, short reads (LSK) and long reads (LK)
#   1) env lr-leafcutter (R 4.3, rstan 2.32, regtools) + leafcutter R package installed from the repo clone
#   2) exon file from GENCODE vM39 (gene names only; clustering is annotation-free)
#   3) per dataset: regtools junctions extract -> leafcutter_cluster_regtools.py -> leafcutter_ds.R
#      Strand: short reads dUTP (-s RF); long reads get an XS tag from the read flag (FLNC + minimap2 -uf) and -s XS,
#      because regtools RF/FR assume paired reads and mislabelled every long-read junction strand (0% GT-AG; fixed 2026-10-05).
#   Groups file lists WT first: leafcutter_ds.R codes group_names[1] = 0 (baseline), deltapsi = perturbed − baseline
#   -> deltapsi = Q157R − WT (rule 8; checked in 19_phase7b_leafcutter.R with the U2af1-213 junction).
#   n = 3 vs 3: -i 3 (intron in ≥ 3 samples), -g 3 (≥ 3 samples per group with coverage).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"   # submit from project root (scripts/submit.sh)
O=${RES}/07_leafcutter; mkdir -p "$O"; T=$SLURM_CPUS_PER_TASK
module load miniforge >/dev/null 2>&1; source "$(conda info --base)/etc/profile.d/conda.sh"

# === 1) environment + R package ===
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_LEAFCUTTER"; then
    echo "[$(date)] creating $ENV_LEAFCUTTER"
    conda create -y -n "$ENV_LEAFCUTTER" -c conda-forge -c bioconda r-base=4.3 r-rstan=2.32 r-stanheaders=2.32 \
        r-rstantools r-rcpp r-rcppeigen r-bh r-ggplot2 r-r.utils r-gridextra r-reshape2 r-hmisc r-foreach r-dplyr \
        r-magrittr r-stringr r-domc r-optparse r-shiny r-shinyjs r-dt r-intervals r-shinycssloaders r-gtable \
        bioconductor-dirichletmultinomial r-data.table regtools python=3.11
fi
conda activate "$ENV_LEAFCUTTER"
if ! Rscript -e 'suppressMessages(library(leafcutter))' 2>/dev/null; then
    echo "[$(date)] installing leafcutter R package (compiles Stan models)"
    # TailRank is not on conda: CRAN, with Biobase (conda) + oompaBase/oompaData (CRAN) as dependencies
    Rscript -e 'quit(status = !requireNamespace("Biobase", quietly = TRUE))' || \
        conda install -y -n "$ENV_LEAFCUTTER" -c conda-forge -c bioconda bioconductor-biobase
    # current oompaBase/oompaData need R >= 4.4 (env is R 4.3 for rstan 2.32): newest archived versions supporting R 4.3
    A=https://cloud.r-project.org/src/contrib/Archive
    Rscript -e 'quit(status = !requireNamespace("oompaBase", quietly = TRUE))' || Rscript -e "install.packages('${A}/oompaBase/oompaBase_3.2.9.tar.gz', repos = NULL, type = 'source')"
    Rscript -e 'quit(status = !requireNamespace("oompaData", quietly = TRUE))' || Rscript -e "install.packages('${A}/oompaData/oompaData_3.1.4.tar.gz', repos = NULL, type = 'source')"
    Rscript -e 'if (!requireNamespace("TailRank", quietly = TRUE)) install.packages("TailRank", repos = "https://cloud.r-project.org")'
    # local patch (repo clone only): the 2018 rstantools Makevars lacks the TBB/RcppParallel + StanHeaders flags that
    # rstan 2.32 needs, and the shipped stanExports_*.cc come from an old stanc -> regenerate, current Makevars template
    LP=${LEAFCUTTER_DIR}/leafcutter
    if ! grep -q RcppParallel "${LP}/src/Makevars"; then
        rm -f "${LP}"/src/*.o "${LP}"/src/*.so
        ( cd "$LP" && Rscript -e 'rstantools::rstan_config()' )
        cat > "${LP}/src/Makevars" <<'MK'
# local patch: current rstantools Makevars template (TBB via RcppParallel, StanHeaders flags, C++17)
STANHEADERS_SRC = $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "message()" -e "cat(system.file('include', 'src', package = 'StanHeaders', mustWork = TRUE))" -e "message()" | grep "StanHeaders")
STANC_FLAGS = $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "cat(ifelse(utils::packageVersion('rstan') >= 2.26, '-DUSE_STANC3',''))")
PKG_CPPFLAGS = -I"../inst/include" -I"$(STANHEADERS_SRC)" -DBOOST_DISABLE_ASSERTS -DEIGEN_NO_DEBUG -DRCPP_PARALLEL_USE_TBB=1 $(STANC_FLAGS) -D_HAS_AUTO_PTR_ETC=0
PKG_CXXFLAGS = $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "RcppParallel::CxxFlags()") $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "StanHeaders:::CxxFlags()")
PKG_LIBS = $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "RcppParallel::RcppParallelLibs()") $(shell "$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "StanHeaders:::LdFlags()")
CXX_STD = CXX17
MK
        printf "leafcutter local patch: src/Makevars -> current rstantools template (TBB/RcppParallel, StanHeaders flags, C++17); stanExports regenerated with rstantools::rstan_config() under rstan 2.32. CRAN-archived oompaBase 3.2.9 / oompaData 3.1.4 (R 4.3).\n" >> "${META}/software_versions.txt"
    fi
    MAKEFLAGS="-j${T}" R CMD INSTALL "$LP"
fi
grep -q "^## ${ENV_LEAFCUTTER}" "${META}/software_versions.txt" || printf "## %s\nleafcutter\t%s (github %s)\nregtools\t%s\nrstan\t%s\n" \
    "$ENV_LEAFCUTTER" "$(Rscript -e 'cat(as.character(packageVersion("leafcutter")))')" \
    "$(git -C "$LEAFCUTTER_DIR" rev-parse --short HEAD)" "$(regtools --version 2>&1 | grep -i -m1 version | awk '{print $NF}')" \
    "$(Rscript -e 'cat(as.character(packageVersion("rstan")))')" >> "${META}/software_versions.txt"

# === 2) exon file (gene names for clusters) ===
EXONS=${O}/gencode_vM39_exons.txt.gz
if [ ! -s "$EXONS" ]; then   # gtf_to_exons.R reads its input through zcat: give it a gzipped exon-only GTF
    awk -F'\t' '$3=="exon"' "$GTF" | gzip > "${O}/gencode_exon_lines.gtf.gz"
    Rscript "${LEAFCUTTER_DIR}/scripts/gtf_to_exons.R" "${O}/gencode_exon_lines.gtf.gz" "$EXONS"
    rm -f "${O}/gencode_exon_lines.gtf.gz"
fi

# === 3) per dataset ===
SAMTOOLS=$(conda run -n "$ENV_ALIGN" which samtools)   # samtools from the alignment env (not in lr-leafcutter)
[ -x "$SAMTOOLS" ] || { echo "samtools not found in $ENV_ALIGN"; exit 1; }
run_set() {   # name  strand  min_clu_reads  samples_tsv  bam_pattern(%s = sample)
    local N=$1 STR=$2 MINR=$3 TSV=$4 PAT=$5 D=${O}/$1
    mkdir -p "${D}/junc"
    : > "${D}/juncfiles.txt"
    for s in $(awk -F'\t' 'NR>1{print $1}' "$TSV"); do
        local B; B=$(printf "$PAT" "$s"); [ -s "$B" ] || { echo "missing $B"; exit 1; }
        if [ ! -s "${D}/junc/${s}.junc" ]; then
            if [ "$STR" = "XSFLAG" ]; then
                # unpaired long reads: regtools RF/FR assume mates and mislabel strands (all "+" became "-", "-" became "?").
                # FLNC reads are transcript-oriented and aligned with -uf, so read strand = transcript strand: add XS from the flag.
                XB=${D}/junc/${s}.xs.bam
                "$SAMTOOLS" view -h -F 0x904 "$B" | awk 'BEGIN{OFS="\t"} /^@/{print; next} {print $0, (int($2/16)%2==1 ? "XS:A:-" : "XS:A:+")}' \
                    | "$SAMTOOLS" view -b -o "$XB" - && "$SAMTOOLS" index "$XB"
                regtools junctions extract -a 8 -m 50 -M 500000 -s XS -o "${D}/junc/${s}.junc" "$XB"
                rm -f "$XB" "${XB}.bai"
            else
                regtools junctions extract -a 8 -m 50 -M 500000 -s "$STR" -o "${D}/junc/${s}.junc" "$B"
            fi
        fi
        echo "${D}/junc/${s}.junc" >> "${D}/juncfiles.txt"
        echo "[$(date)] ${N} ${s}: $(wc -l < "${D}/junc/${s}.junc") junctions"
    done
    # groups: WT first (baseline), then Q157R
    { awk -F'\t' 'NR>1 && $2=="WT"{print $1"\tWT"}' "$TSV"; awk -F'\t' 'NR>1 && $2=="Q157R"{print $1"\tQ157R"}' "$TSV"; } > "${D}/groups.txt"
    if [ ! -s "${D}/${N}_perind_numers.counts.gz" ]; then
        ( cd "$D" && python "${LEAFCUTTER_DIR}/clustering/leafcutter_cluster_regtools.py" -j juncfiles.txt -m "$MINR" \
              -o "$N" -l 500000 -r "$D" )
    fi
    if [ ! -s "${D}/${N}_ds_cluster_significance.txt" ]; then
        Rscript "${LEAFCUTTER_DIR}/scripts/leafcutter_ds.R" --num_threads "$T" -i 3 -g 3 -e "$EXONS" \
            -o "${D}/${N}_ds" "${D}/${N}_perind_numers.counts.gz" "${D}/groups.txt"
    fi
    awk -F'\t' 'NR>1{n[$2]++; if($2=="Success" && $6<0.05) s++} END{for(k in n) print "  status", k, n[k]; print "  p.adjust<0.05:", s+0}' \
        "${D}/${N}_ds_cluster_significance.txt"
}
run_set shortread RF 50 "$SR_SAMPLES" "${RES}/07_shortread/star/%s.Aligned.sortedByCoord.out.bam"
run_set longread  XSFLAG 30 "${PROJ}/samples.tsv" "${RES}/01_align/bam/%s.sorted.bam"
echo "[$(date)] 07e complete"
