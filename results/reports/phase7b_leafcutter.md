# Phase 7b — LeafCutter (annotation-free intron usage), short and long reads (2026-10-05)

**Sign convention:** ΔPSI is **Q157R − WT**. The groups file lists WT first, which `leafcutter_ds.R` codes as the
baseline (0), so `deltapsi = perturbed − baseline` = Q157R − WT. The summary script checks the effect-size column order
(`WT`, `Q157R`) and the U2af1-213 positive control before anything else.
Outputs: `results/07_leafcutter/compare/` (tracked); junction files, counts and raw LeafCutter tables stay local.

## Why LeafCutter
Every splicing result so far depends on a transcript annotation: rMATS-long and DRIMSeq test the long-read isoforms;
rMATS-turbo tests events built from the expanded GTF. LeafCutter clusters introns directly from split reads, so it
is independent of ESPRESSO, SQANTI3 and GENCODE. It is run on both datasets: short reads (LSK cohort, STAR BAMs from
Phase 7) and long reads (LK, Phase 1 minimap2 BAMs).

## Jobs
| Step | Script | SLURM ID | Outcome |
|---|---|---|---|
| env + R package; regtools → clustering → `leafcutter_ds` for both datasets | `07e_leafcutter.slurm` | 20900963 (short reads), 20902429 (long-read rerun) | COMPLETED |
| summary, rule-8 checks, 3′SS tests, comparisons | `07f_leafcutter_summary.slurm` → `19_phase7b_leafcutter.R` | 20904579 | COMPLETED |

Settings: regtools `-a 8 -m 50 -M 500000`; clusters ≥ 50 reads (short) / ≥ 30 (long); `leafcutter_ds -i 3 -g 3`
(n = 3 per group: below the n = 4 the authors calibrated p-values for, so FDR is read with caution).

Problems fixed on the way:
- **Environment:** LeafCutter 0.2.9 (R package, github 2c9907e) needs `TailRank`, whose current CRAN dependencies
  require R ≥ 4.4 while rstan 2.32 pins R 4.3. Used CRAN-archived oompaBase 3.2.9 / oompaData 3.1.4. The package's
  2018 `Makevars` lacks the TBB/StanHeaders flags current rstan needs: patched (repo clone only) with the current
  rstantools template and regenerated Stan exports. Recorded in `software_versions.txt`.
- **Long-read strands (20900963):** regtools' `RF`/`FR` modes assume paired reads. On the unpaired FLNC reads every
  plus-strand junction was labelled "−" and every minus-strand one "?", giving **0% GT-AG** at the summary step's
  coordinate check (flipping restored 98%). Fixed by adding an `XS` tag from the read flag (FLNC reads are
  transcript-oriented and aligned with `-uf`, so read strand = transcript strand) and running regtools `-s XS`.
  After the fix: 97.6% (short) / 98.2% (long) of sampled introns are GT…AG on the reported strand.
- Two summary-script bugs (a column shadowing a data frame inside `tibble()`; LRP2 intron IDs with comma-separated
  multi-intron entries) fixed; no result was used from those runs.

**Rule-8 positive control (U2af1-213):** the intron whose 5′SS is created by the Q157R base
(chr17:31,867,070–31,867,168, − strand) is in a tested cluster in both datasets and is **up in Q157R**:

| Dataset | PSI WT | PSI Q157R | ΔPSI | adj p |
|---|---|---|---|---|
| Short reads | 0.001 | 0.100 | **+0.098** | 0.015 |
| Long reads | 0.001 | 0.120 | **+0.120** | 0.001 |

## 1) Differential intron usage
| Dataset | Clusters tested | Clusters adj p < 0.05 | Genes | Introns in those clusters with \|ΔPSI\| ≥ 0.10 | …up in Q157R |
|---|---|---|---|---|---|
| Short reads (LSK) | 7,270 | **168** | 183 | 214 | 108 |
| Long reads (LK) | 13,543 | **295** | 298 | 127 | 79 |

- The earlier LRP2 long-read LeafCutter run found 0 significant clusters (11,388 tested). This run's long-read input
  differs (new alignments, per-read XS strands, ≥ 30-read clusters); across 18,458 shared introns its ΔPSI
  (negated to Q157R − WT) correlates with ours at **r = 0.79**.

## 2) U2AF1 Q157 signature, independent of any annotation
Within significant clusters, every donor with two alternative acceptors was taken; the acceptor with the larger
ΔPSI is "favoured in Q157R", the one with the smaller is "disfavoured". Fisher's exact test on the base after the AG
(+1) and the base before it (−3, the S34F signature, negative control):

| Dataset | Acceptor pairs | +1 G, favoured | +1 G, disfavoured | OR | p |
|---|---|---|---|---|---|
| Long reads | 356 | **70.7%** | **19.8%** | **9.7** | **4.9 × 10⁻³⁶** |
| Short reads | 132 | **69.1%** | **38.9%** | **3.5** | **4.4 × 10⁻⁵** |

| Dataset | −3 C, favoured | −3 C, disfavoured | OR | p |
|---|---|---|---|---|
| Long reads | 61.8% | 63.0% | 0.95 | 0.75 |
| Short reads | 62.0% | 67.2% | 0.80 | 0.42 |

- 99–100% of both acceptor sets carry AG at −2/−1 (coordinate check).
- The Phase 5 result (82.5% vs 44.9%, OR 5.7, on ESPRESSO novel junctions) therefore does not depend on the
  assembler or the annotation: the same preference appears from raw split reads, in both cohorts and both cell types.
  The short-read effect here (OR 3.5) is stronger than the rMATS-turbo within-event test gave in Phase 7 (OR 1.25),
  because LeafCutter restricts to acceptors that actually share a donor in a significant cluster.

## 3) Short vs long reads, and vs the other methods
- 15,045 introns are tested in both datasets: ΔPSI r = 0.21 overall; among the 244 introns significant in either,
  **r = 0.52 and 75% same direction**; 19 are significant in both.
- Genes with a significant LeafCutter cluster that also have a significant rMATS-turbo event (same short-read BAMs):
  117 of 183 (short reads), 202 of 298 (long reads).
- **Phase 8 tier A genes:** 16 of 36 have a significant LeafCutter cluster in long reads, 2 of 36 in short reads
  (LeafCutter's cluster definition and its n = 3 power differ from the isoform-level tests, so this is a floor).
- **Cd34, the same last-exon A3SS** (donor chr1:194,641,539; acceptors 194,642,088 vs 194,642,244):

| Dataset | PSI upstream acceptor, WT → Q157R | ΔPSI | adj p |
|---|---|---|---|
| Long reads | 0.42 → 0.18 | **−0.24** (+0.24 downstream) | 2.5 × 10⁻⁴ |
| Short reads | 0.27 → 0.11 | **−0.17** (+0.17 downstream) | 0.022 |

  Fourth method (after rMATS-long, DRIMSeq, rMATS-turbo) and second cohort showing the shift to the downstream
  acceptor.

## Conclusions
- An annotation-free method reproduces the two main splicing results: the U2af1-213 cis junction is up, and
  Q157R-favoured 3′ splice sites prefer +1 G (OR 9.7 long reads, 3.5 short reads) with a null −3 control.
- Cd34 replicates again. Long- and short-read ΔPSI agree in direction for 75% of significant introns.

## Caveats
- n = 3 per group is below LeafCutter's calibrated range (n ≥ 4); treat adj p as approximate.
- Long reads were given to a tool designed for short reads; clustering uses split alignments only, so read length
  does not matter, but the per-read coverage model does.
- Different animals and cells (LSK vs LK) in the short-read cohort.

## Outputs (`results/07_leafcutter/compare/`)
- `phase7b_leafcutter_summary.tsv`, `phase7b_{shortread,longread}_significant_introns.tsv`
- `phase7b_3ss_tests.tsv`, `phase7b_alt_acceptor_pairs.tsv`, `phase7b_U2af1_213_control.tsv`
- `phase7b_short_vs_long.tsv`, `phase7b_short_vs_long_dpsi.{pdf,png}`, `phase7b_longread_vs_LRP2.tsv`
- `phase7b_vs_rmats_turbo_genes.tsv`, `phase7b_tierA_leafcutter.tsv`, `phase7b_Cd34_introns.tsv`
