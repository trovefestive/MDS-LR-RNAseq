# Phase 7 — Short-read cross-check report (2026-10-04)

**Sign convention:** every effect (logFC, ΔPSI, Δ isoform proportion) is **Q157R − WT** (WT = reference level).
Outputs: `results/07_shortread/compare/` (tracked); Salmon, STAR and rMATS-turbo outputs stay local.

## Dataset
- A separate in-house Illumina dataset: paired-end 150 bp, dUTP-stranded, 3 WT (SR_WT1–3) and 3 U2af1 Q157R knock-in
  (SR_KI1–3) mice. See `samples_shortread.tsv`.
- **Cell type differs:** the short reads are from **LSK** (Lin⁻ Sca1⁺ Kit⁺) cells, while the long reads are from **LK**
  (Lin⁻ Kit⁺) cells. These are different animals and a different, more primitive, population, so this checks
  replication of mutation effects, not matched validation.
- 56–94 M read pairs per sample.

## Jobs
| Step | Script | SLURM ID | Outcome |
|---|---|---|---|
| Salmon env, expanded GTF/FASTA (171,224 isoforms), Salmon decoy index, STAR index | `07a_sr_setup.slurm` | 20828827 | COMPLETED (1 h 41 min) |
| Salmon + STAR + genotype, per sample | `07b_sr_quant_align.slurm` | 20828828 (+ 20841172 for 2 Salmon reruns) | COMPLETED |
| rMATS-turbo prep, per BAM | `07c_sr_rmats_prep.slurm` | 20841173 | COMPLETED (3–11 min each) |
| rMATS-turbo post | `07c_sr_rmats.slurm` | 20841174 | COMPLETED |
| Comparison | `07d_sr_compare.slurm` → `17_phase7_compare.R` | 20842456 | COMPLETED |

Problems fixed on the way:
- **20828835:** the single all-BAM rMATS-turbo run segfaulted after 33 s with no message. Reduced tests (full GTF,
  8 threads, unplaced contigs) all ran cleanly, so the run was split into per-BAM `--task prep` jobs plus one
  `--task post` (rMATS's documented route). That completed without errors.
- **Library type:** all 6 libraries are dUTP-stranded (86–92% of strand-informative Salmon fragments are ISR), but
  Salmon's `-l A` auto-detection called SR_WT2 and SR_KI3 "IU". The type is now fixed in `config.sh` (Salmon `ISR`,
  rMATS `fr-firststrand`), and those two samples were requantified. `07c` checks the ISR fraction (≥ 0.8) of every sample.
- **Salmon duplicate collapsing:** Salmon merges isoforms with identical spliced sequence. 170 of 171,224 (51 novel)
  were merged into another isoform. They get 0-count rows (gene sums unchanged), are excluded from the detection rate,
  and both members of a pair are excluded from isoform-usage comparisons.
- **20841175, 20842359:** 07d stopped on the duplicate IDs above, then on two R bugs (detection % computed from the
  summed count; a `tibble()` column shadowing a data frame). Fixed; no result was used from those runs.

**Rule-8 checks:**
- rMATS-turbo: `--b1` = Q157R, `--b2` = WT; `IncLevelDifference` = mean PSI b1 − mean PSI b2 (checked on Cd34:
  0.096 − 0.256 = −0.16).
- Positive control U2af1-213 (ENSMUST00000468653.1), Salmon reads: Q157R 205 / 246 / 116 (TPM 10.3 / 13.9 / 5.6) vs
  WT 34 / 0 / 6 (TPM 1.5 / 0 / 0.4) → up in Q157R, as in the long reads.
- Short-read gene DE uses `~ genotype` with WT as the reference level (asserted in the script).

## 1) Genotype
| Sample | Expected | Depth | Alt C fraction | Call |
|---|---|---|---|---|
| SR_WT1 | WT | 1,105 | 0.001 | WT |
| SR_WT2 | WT | 324 | 0.000 | WT |
| SR_WT3 | WT | 801 | 0.002 | WT |
| SR_KI1 | Q157R | 830 | 0.394 | Q157R |
| SR_KI2 | Q157R | 668 | 0.431 | Q157R |
| SR_KI3 | Q157R | 1,475 | 0.411 | Q157R |

U2af1 chr17:31,867,169 T>C. All six match; the knock-ins are heterozygous.

## 2) Mapping
- STAR (GRCm39, expanded-GTF junctions): 67–73% uniquely mapped, 3–4% multi-mapped.
- Salmon (expanded transcriptome, genome decoys): 35–50% mapped. 16–23 M fragments per sample hit the genome decoy
  (introns/intergenic), consistent with pre-mRNA-rich libraries rather than an index problem.

## 3) Are the long-read isoforms seen in short reads?
Detected = ≥ 1 read in ≥ 2 samples and ≥ 5 reads in total (Salmon).

| Origin | Isoforms | Detected | % |
|---|---|---|---|
| GENCODE | 126,198 | 80,270 | 63.6 |
| Novel | 44,856 | 27,527 | 61.4 |

Novel isoforms are detected as often as GENCODE ones. Caveat: Salmon distributes ambiguous reads by EM, so detection
is compatible-read support, not proof of the full novel structure.

## 4) Gene-level expression: does not replicate
- 10,854 genes shared; long- vs short-read logFC r = **0.18**.
- **0 genes** at FDR < 0.05 in short reads (minimum FDR 0.88).
- Of the 97 long-read DE genes, 30 have the same sign in short reads and none is significant.
- Of the 34 ribosomal-protein genes up in long reads, 9 go the same way (median short-read logFC −0.15).
- The ribosomal/OXPHOS increase of Phase 6 is therefore **not replicated**. Possible reasons: a library/length effect
  in the long reads (the Phase 6 caveat), or the different cell population (LSK vs LK). This analysis cannot tell them apart.

## 5) Isoform-usage switches: replicate
Long-read high-confidence switches = rMATS-long ∩ DRIMSeq (122 isoforms, Phase 6). Short-read Δ = Salmon isoform
share Q157R − WT within the same gene grouping.

| High-confidence switches | Testable (≥ 10 gene reads) | Same direction | r (Δ) | Binomial p |
|---|---|---|---|---|
| 122 | 121 | 83 (69%) | 0.44 | 5.3e-5 |

| Gene | Isoform | Long-read Δ | Short-read WT → Q157R | Short-read Δ |
|---|---|---|---|---|
| Cd34 | Cd34-201 (ENSMUST00000016638) | +0.22 | 0.56 → 0.64 | +0.08 |
| Cd34 | Cd34-202 (ENSMUST00000110815) | −0.23 | 0.11 → 0.00 | −0.11 |
| Tmpo | ENSMUST00000072239 | −0.18 | 0.32 → 0.21 | −0.11 |
| Atrx | novel ESPRESSO:chrX:1412:82 | +0.16 | 0.10 → 0.10 | −0.01 (not replicated) |

## 6) rMATS-turbo (short reads, expanded GTF, `--novelSS`)
| Event | Tested (JCEC) | FDR < 0.05, \|ΔPSI\| ≥ 0.1 | …higher inclusion in Q157R |
|---|---|---|---|
| SE | 29,704 | 2,852 | 1,122 |
| RI | 19,341 | 2,671 | 1,025 |
| A3SS | 21,197 | 1,624 | 762 |
| A5SS | 13,782 | 1,144 | 519 |
| MXE | 5,523 | 595 | 290 |

rMATS-turbo is known to be liberal with n = 3 per group and about 60 M read pairs per sample; read the counts as an
upper bound and rely on the agreement tests.

- **Cd34, independent replication of the top long-read switch:** the same last-exon A3SS (acceptors chr1:194,642,088
  and 194,642,244). PSI of the upstream acceptor is 0.094 / 0.122 / 0.071 in Q157R vs 0.262 / 0.278 / 0.227 in WT,
  ΔPSI **−0.16**, FDR ≈ 0, with no overlap between replicates. Q157R shifts usage to the downstream acceptor (Cd34-201),
  as in the long reads, now in a second cell population (LSK) and a second cohort.
- **U2AF1 Q157 +1 G signature:** in the 1,624 significant A3SS events, the acceptor favoured in Q157R has G at +1 in
  **61.3%** of events vs **55.9%** for the disfavoured acceptor (Fisher OR **1.25**, p = **0.013**). 99% of both
  acceptor sets have AG at −2/−1 (coordinate check).
  - Same direction as Phase 5 but much weaker (long reads: 82.5% vs 44.9%, OR 5.70). The tests differ: here the two
    acceptors within each event are compared; Phase 5 compared novel junctions associated with Q157R vs with WT.

## Conclusions
- Genotypes and the U2af1-213 positive control replicate.
- **Splicing effects replicate across cohorts and cell types** (isoform-usage direction 69%, p = 5e-5; Cd34 A3SS;
  a weaker +1 G preference), while **gene-level expression changes do not**. That fits a mutation-intrinsic
  splicing effect plus population-specific (or technical) expression differences.
- Novel long-read isoforms are as detectable in short reads as GENCODE isoforms (61% vs 64%).

## Caveats
- Different animals and different cells (LSK vs LK): this is replication, not matched validation.
- n = 3 vs 3 in both datasets; rMATS-turbo counts are liberal.
- Salmon detection of novel isoforms relies on EM read assignment.

## Tables (`results/07_shortread/compare/`)
- `phase7_genotype_check.tsv`, `phase7_salmon_mapping.tsv`, `phase7_isoform_detection.tsv`
- `phase7_gene_DE_lr_vs_sr.tsv`, `phase7_gene_DE_summary.tsv`
- `phase7_usage_highconf_lr_vs_sr.tsv`, `phase7_usage_summary.tsv`
- `phase7_rmats_turbo_events.tsv`, `phase7_rmats_A3SS_significant.tsv`, `phase7_rmats_A3SS_plus1G.tsv`
