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

## Status
Phase 0 (setup) complete 2026-09-30. See results/reports/phase0_setup.md.
Phase 1 (QC + alignment) complete 2026-09-30. See results/reports/phase1_qc_alignment.md.

## Key results
- Genotypes confirmed from the reads: U2af1 Q157R (chr17:31,867,169 T>C) at 40–42% allele fraction in A258/A309/X504, and 0–0.2% in WT.
- Mapping ≥99.94% in all samples; HiFi mismatch rate 0.21–0.22%.

## Earlier analysis (do not reuse its signs)
`2026_04/260327_MDS_lrp2` holds the LRP2 pipeline outputs from April 2026. Its edgeR logFC values are WT − Q157R even
though the files are named "Q157R_vs_WT", so every direction call in them is reversed.
