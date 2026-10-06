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
IsoQuant (cross-check), SQANTI3, rMATS-long, and R (edgeR, DRIMSeq, IsoformSwitchAnalyzeR). Phase 7 short reads: Salmon,
STAR, rMATS-turbo 4.1.2. Phase 7b: LeafCutter 0.2.9 + regtools (annotation-free).

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
| scripts/04a_nmd.slurm → 13_phase4_nmd.R | 4 | NMD known vs novel (Fisher), per-sample NMD expression share (Welch t), per-gene NMD share (edgeR diffSpliceDGE), rule-8 sign checks |
| scripts/05a_splice_signatures.slurm → 14_phase5_splice_signatures.R | 5 | genotype-specific/enriched novel isoforms (edgeR diffSpliceDGE), junction classes, U2AF1 3′SS tests (+1 G/A; −3 C/T control), logos, GO ORA |
| scripts/06a_rmats_long.slurm | 6 | rMATS-long (group 1 = Q157R, group 2 = WT → Δ = Q157R − WT), single-threaded |
| scripts/06b_diff_isoforms.slurm → 15_phase6_diff_isoforms.R | 6 | edgeR DGE/DTE, DRIMSeq DTU, switches, LRP2 comparison (negated), Ybx1, dot plots, rule-8 checks |
| scripts/06c_rmats_summary.slurm → 16_phase6_rmats_summary.R | 6 | rMATS-long import, direction checks, DRIMSeq agreement, event classes |
| scripts/07a_sr_setup.slurm | 7 | expanded GTF (final 171,224 isoforms) → transcript FASTA, Salmon decoy-aware index, STAR index |
| scripts/07b_sr_quant_align.slurm | 7 | 6-task array: Salmon (ISR), STAR to GRCm39, U2af1 Q157R genotype |
| scripts/07c_sr_rmats_prep.slurm | 7 | 6-task array: rMATS-turbo `--task prep` per BAM |
| scripts/07c_sr_rmats.slurm | 7 | rMATS-turbo `--task post` (b1 = Q157R, b2 = WT → ΔPSI = Q157R − WT), expanded GTF, `--novelSS` |
| scripts/07d_sr_compare.slurm → 17_phase7_compare.R | 7 | genotype, novel-isoform detection, gene DE and isoform-usage replication, A3SS +1 G test |
| scripts/08a_summary.slurm → 18_phase8_candidates.R | 8 | isoform-switch candidates with evidence tiers (both methods, replicate separation, read support, IsoQuant, short reads, NMD) |
| scripts/07e_leafcutter.slurm | 7b | LeafCutter 0.2.9: regtools junctions (short reads RF; long reads XS from read flag) → clustering → `leafcutter_ds` (WT baseline → ΔPSI = Q157R − WT), both datasets |
| scripts/07f_leafcutter_summary.slurm → 19_phase7b_leafcutter.R | 7b | coordinate + rule-8 checks, significant clusters, +1 G / −3 C tests on shared-donor acceptor pairs, short vs long, LRP2, rMATS-turbo, tier A, Cd34 |

## Key results
- Genotypes confirmed from the reads: U2af1 Q157R (chr17:31,867,169 T>C) at 40–42% allele fraction in A258/A309/X504, and 0–0.2% in WT.
- Mapping ≥99.94% in all samples; HiFi mismatch rate 0.21–0.22%.
- ESPRESSO (unfiltered): 253,272 isoforms (126,317 GENCODE, 126,955 novel); 81–85% of reads assigned; novel isoforms carry 6.2–6.8% of reads in every sample.
- The "mutant-exclusive" U2af1 isoform is GENCODE U2af1-213. It is present only in Q157R (UP, not down), and every read using its junction carries the mutant allele: the Q157R base creates a 5′ splice site in U2af1 itself (cis effect).
- Final transcriptome (Phase 3): 171,224 isoforms = 126,317 GENCODE + 44,907 new (rules filter; ≥2 samples, ≥5 reads). New isoforms are 4.2–4.7% of expression in every sample; 36% of new coding isoforms are predicted NMD vs 14% of GENCODE. Primary counts = all reads (R2); the paper's FSM-read filter removes 15–20% of reads.
- NMD (Phase 4): 36.0% of novel coding isoforms are predicted NMD vs 14.0% of GENCODE (OR 3.44). The NMD share of expression is the same in Q157R and WT (1.88% vs 1.91%), and no gene shows a significant NMD-share change (5,907 tested; n = 3 vs 3).
- **U2AF1 Q157 splice signature (Phase 5):** new 3′ splice sites gained in Q157R have G at +1 in 82.5% of cases vs 44.9% for WT-associated sites (Fisher OR 5.70, p = 2.2e-6). The S34F-type −3 C/T control shows no difference (p = 1).
- **Differential isoforms (Phase 6):** rMATS-long finds 866 isoform-usage changes in 576 genes, DRIMSeq 130. The 122 found by both have the same direction in every case (r = 0.90). Most switches have no gene-level change. Top switch: Cd34 (last-exon alternative 3′SS, Cd34-201 +0.22). Gene level: 97 genes DE (89 up, mostly ribosomal/OXPHOS). Results replicate the earlier LRP2 run once its effect sizes are negated (r = 0.95; U2af1-213 +8.86 vs +8.80). No Ybx1 shift.
- A310 (WT) has about 2× Mpo/Elane, which suggests more promyelocyte-like cells in that sort. Treat granule-gene differences with caution.
- **Short-read cross-check (Phase 7; separate cohort, LSK cells, 3 WT vs 3 Q157R):** genotypes and U2af1-213 (up in Q157R)
  replicate. Isoform-usage switches replicate (83/121 same direction, binomial p = 5e-5; Cd34 A3SS ΔPSI −0.16, FDR ≈ 0, same
  acceptor shift). Gene-level changes, including the ribosomal/OXPHOS increase, do not (r = 0.18; 0 DE genes). Q157R-favoured
  A3SS acceptors have +1 G more often (61% vs 56%, p = 0.013). 61% of novel isoforms are detected (GENCODE 64%).
- **LeafCutter, annotation-free (Phase 7b):** 295 significant clusters in long reads, 168 in short reads. The U2af1-213 cis junction
  is up in both (ΔPSI +0.12 / +0.10). Q157R-favoured acceptors (shared donor, significant clusters) have +1 G in 70.7% vs 19.8%
  (OR 9.7, p = 5e-36, long reads) and 69.1% vs 38.9% (OR 3.5, p = 4e-5, short reads); −3 C/T control null in both. Cd34 A3SS
  replicates a fourth time (ΔPSI −0.24 long, −0.17 short, upstream acceptor).
- **Follow-up candidates (Phase 8):** 918 isoforms in 611 genes are significant in rMATS-long or DRIMSeq. **36 tier A genes** are
  significant in both methods, separate all replicates in long reads, and replicate in direction in short reads, with IsoQuant
  not contradicting. Top: Cd34 (A3SS, also by short-read rMATS-turbo), Cdc14a, Mtm1, Erbin, Crlf3, Atn1, Thada, Tcf19, Tmpo.
  None of the tier A/B isoforms is predicted NMD; 31 of 36 tier A genes have no gene-level change.

## Main figures
| Figure | Shows |
|---|---|
| results/01_align/phase1_read_length_distribution, phase1_depth_mapping | read length, depth, mapping QC |
| results/03_sqanti/phase3_novel_categories | structural categories of new isoforms |
| results/03_sqanti/phase3_expression_known_vs_novel | expression, known vs new |
| results/05_splice/phase5_3ss_logos | U2AF1 Q157 3′ splice-site signature (+1 G) |
| results/06_diff/phase6_dotplot_U2af1 | U2af1-213, the cis Q157R isoform, every replicate |
| results/06_diff/phase6_rmatslong_Cd34_structure / _abundance | top switch, Cd34 A3SS |
| results/08_summary/phase8_lr_vs_sr_delta | long- vs short-read Δ for all candidates |
| results/08_summary/phase8_evidence_top_genes | evidence grid, top 30 candidates |
| results/08_summary/phase8_tierA_proportions | tier A switches, every replicate, both datasets |
| results/07_leafcutter/compare/phase7b_short_vs_long_dpsi | LeafCutter ΔPSI, long vs short reads |

The candidate list with evidence columns is `results/08_summary/phase8_candidates_genes.tsv` (one row per gene) and
`phase8_candidates_isoforms.tsv`; see `results/reports/phase8_summary.md`.
