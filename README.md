# PacBio long-read RNA-seq — U2af1 Q157R vs WT (mouse LK)

The detailed analysis plan is kept locally and is not part of this repository. Phase reports are in `results/reports/`.

## Samples
Kinnex full-length RNA, Revio run r84050_20240510; FLNC reads delivered by the core (skera → lima → refine --require-polya).
See [samples.tsv](samples.tsv).

| Genotype | Samples | FLNC reads |
|---|---|---|
| WT | V335, V334, A310 | 10.27 M, 5.68 M, 9.15 M |
| Q157R | X504, A258, A309 | 7.16 M, 9.83 M, 8.30 M |

## Comparison
Q157R vs WT, with WT as the reference. **Every effect size (logFC, ΔPSI, Δproportion) is Q157R − WT.**

## Reference
GRCm39 primary assembly + GENCODE vM39 (481,956 transcripts, much larger than older releases, so known/novel counts are not comparable with older annotations). SQANTI3 support comes from refTSS v4.1 (mm39) and PolyASite 2.0 lifted
mm10 → mm39.

## Tools
Versions are in `results/metadata/software_versions.txt`. Main tools: minimap2, ESPRESSO (primary assembler),
IsoQuant (cross-check), SQANTI3, rMATS-long, and R (edgeR, DRIMSeq, IsoformSwitchAnalyzeR).

## Scripts
| Script | Phase | Purpose |
|---|---|---|
| scripts/submit.sh | all | submit any step from the project root: `scripts/submit.sh scripts/NN_step.slurm` (reads untracked `config.local.sh` for the SLURM account/partition) |
| scripts/00a_setup_envs.slurm | 0 | conda envs + version record |
| scripts/00b_download_refs.slurm | 0 | GENCODE vM39, refTSS, PolyASite, paper code |
| scripts/00c_build_indexes.slurm | 0 | faidx, minimap2 index, junction BED, PolyASite liftOver |
| scripts/00_plot_theme.R | all | shared ggplot theme; WT #3B6FB6, Q157R #C8553D |
| scripts/01a_qc_align.slurm | 1 | 6-task array: seqkit QC, minimap2 splice:hq + junc-bed, sort/index, U2af1 Q157R pileup (chr17:31867169 T>C) |
| scripts/01b_qc_summary.slurm → 10_phase1_qc_summary.R | 1 | QC table, length/depth plots, genotype check (stops on mismatch) |
| scripts/02a_espresso_by_chr.slurm | 2 | 21-task array: per-chromosome ESPRESSO S → C (-T 8, retries) → Q, as in Miller et al. |
| scripts/02b_espresso_merge.slurm | 2 | merge per-chromosome abundance tables and GTFs |
| scripts/02c_isoquant.slurm | 2 | IsoQuant 4.0.0 cross-check (FL, stranded, polyA-trimmed HiFi) |
| scripts/11_phase2_summary.py | 2 | assignment rates, known/novel counts, A310 composition check, U2af1 junction alleles |
| scripts/03a_sqanti3_qc.slurm | 3 | SQANTI3 6.0.2 QC (CAGE, PolyASite, polyA motifs, read counts, TD2 ORFs) |
| scripts/03b_sqanti3_filter.slurm | 3 | ML + rules filters; primary = `P3_FILTER` (rules); restore GENCODE; per-chromosome filtered GTFs |
| scripts/03b2_sqanti3_ml_noORF.slurm + 03_ml_remove_columns.txt | 3 | ML filter without ORF/NMD/per-sample features (comparison) |
| scripts/03c_requant_by_chr.slurm | 3 | 21-task array: ESPRESSO_Q on filtered GTF (R2) → FSM-read filter → requant (fsm) |
| scripts/paper_derived/remove_nonfsm_data.pl | 3 | paper's FSM-read filter (MIT), one-line bug fix |
| scripts/03d_merge_summary.slurm → 12_phase3_summary.R | 3 | merge, read loss + primary-matrix rule, novel rule, TMM CPM, IsoQuant/LRP2 overlap, plots |

## Status
Phase 0 (setup) complete 2026-09-30. See results/reports/phase0_setup.md.
Phase 1 (QC + alignment) complete 2026-09-30. See results/reports/phase1_qc_alignment.md.
Phase 2 (ESPRESSO assembly + IsoQuant) complete 2026-09-30. See results/reports/phase2_assembly.md.
Phase 3 (SQANTI3 curation + requantification) complete 2026-10-02. See results/reports/phase3_curation.md.

## Key results
- Genotypes confirmed from the reads: U2af1 Q157R (chr17:31,867,169 T>C) at 40–42% allele fraction in A258/A309/X504, and 0–0.2% in WT.
- Mapping ≥99.94% in all samples; HiFi mismatch rate 0.21–0.22%.
- ESPRESSO (unfiltered): 253,272 isoforms (126,317 GENCODE, 126,955 novel); 81–85% of reads assigned; novel isoforms carry 6.2–6.8% of reads in every sample.
- The "mutant-exclusive" U2af1 isoform is GENCODE U2af1-213. It is present only in Q157R (UP, not down), and every read using its junction carries the mutant allele: the Q157R base creates a 5′ splice site in U2af1 itself (cis effect).
- Final transcriptome (Phase 3): 171,213 isoforms = 126,317 GENCODE + 44,896 new (rules filter; ≥2 samples, ≥5 reads). New isoforms are 4.2–4.7% of expression in every sample; 36% of new coding isoforms are predicted NMD vs 14% of GENCODE. Primary counts = all reads (R2); the paper's FSM-read filter removes 15–20% of reads.
- A310 (WT) has about 2× Mpo/Elane, which suggests more promyelocyte-like cells in that sort. Treat granule-gene differences with caution.

## Earlier analysis (do not reuse its signs)
`2026_04/260327_MDS_lrp2` holds the LRP2 pipeline outputs from April 2026. Its edgeR logFC values are WT − Q157R even
though the files are named "Q157R_vs_WT", so every direction call in them is reversed.
