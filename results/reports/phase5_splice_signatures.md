# Phase 5 — Genotype-associated novel isoforms and U2AF1 3′ splice-site signatures (2026-10-02)

**Sign convention:** every difference is **Q157R − WT** (WT = reference level).

## Job
`scripts/05a_splice_signatures.slurm` → `scripts/14_phase5_splice_signatures.R`, job 20775376 (final, on the
duplicate-corrected Phase 3 set).
- Earlier attempts 20765910/20766004/20766111 stopped on script bugs: an NA in a rule-8 check, a missing count join, and
  a clusterProfiler/dplyr `select` clash. They gave the same motif result.
- Inputs: the Phase 3 final set (R2 counts), the SQANTI3 junction table, and BSgenome mm39.
- Outputs are in `results/05_splice/`.

## A) Which isoforms go with which genotype
- **Specific:** ≥1 read in ≥2 of the 3 mice of one genotype, and 0 reads in every mouse of the other.
- **Enriched:** the isoform's share of its gene shifts (edgeR `diffSpliceDGE`, isoforms as "exons", `~genotype`). Genes
  need ≥10 reads in ≥3 samples. Calls need FDR < 0.05 and |Δ share| > 0.10.
  - The paper used a Wilcoxon test across many samples, which cannot reach significance with 3 vs 3.

| | specific only | enriched only | both |
|---|---|---|---|
| **Novel**, Q157R | 854 | 3 | 0 |
| **Novel**, WT | 801 | 0 | 2 |
| GENCODE, Q157R | 4,360 | 6 | 7 |
| GENCODE, WT | 4,486 | 1 | 4 |

- With n = 3, formal usage shifts are rare (23 isoforms in total). Almost all associations are **detection patterns**
  ("specific"). Because WT and Q157R have nearly equal numbers, most of these are likely sampling noise at low counts.
  **The motif test below is what tells signal from noise**: noise would not produce a genotype-specific 3′SS sequence.
- **Rule-8 checks passed.** All 4,760 Q157R-specific isoforms with a usage test have Δ share > 0 and logFC > 0. Every
  enriched call has matching signs.
- U2af1-213 (the cis isoform created by Q157R, Phase 2) is called Q157R-specific (usage FDR 1e-11, Δ share +0.099).

## B) Junction classes in genotype-associated novel isoforms (unique junctions)
Classes follow the paper's code: **A3SS** = known donor + novel acceptor; A5SS = novel donor + known acceptor; NE = both
sites novel (novel exon); NJ = known sites in a new pairing (e.g. exon skipping). Site novelty comes from SQANTI3.

| Group | A3SS | A5SS | NE | NJ | Intron-retention isoforms |
|---|---|---|---|---|---|
| Q157R | 127 | 119 | 45 | 133 | 365 of 857 |
| WT | 108 | 146 | 64 | 75 | 372 of 803 |

- Q157R-associated isoforms have more NJ junctions (133 vs 75), which mostly means exon skipping or inclusion. The other
  classes are similar between groups.
- Mutually exclusive exons are not separated from NJ here. Event-level classification is part of Phase 6 (rMATS-long).

## C) U2AF1 3′SS signatures at novel A3SS acceptors (paper method)
Bases −3/−2/−1/+1 around each unique novel A3SS acceptor, strand-aware, from BSgenome mm39. **100% are AG** (the
check stops the script below 95%).

| Acceptor set | n | +1 G / (G+A) | −3 C / (C+T) |
|---|---|---|---|
| Q157R-associated | 127 | **82.5%** (66/80) | 73.4% |
| WT-associated | 108 | 44.9% (31/69) | 72.6% |
| All novel A3SS | 3,726 | 62.6% | 72.6% |
| Annotated 3′SS | 166,124 | 64.1% | 68.4% |

| Test (Fisher's exact) | Odds ratio | p |
|---|---|---|
| **Q157 +1 G vs A: Q157R vs WT** | **5.70** | **2.2 × 10⁻⁶** |
| Q157 +1 G vs A: Q157R vs annotated | 2.64 | 4.1 × 10⁻⁴ |
| Q157 +1 G vs A: WT vs annotated | 0.46 | 1.5 × 10⁻³ |
| **S34 −3 C vs T (negative control): Q157R vs WT** | 1.04 | **1.0** |
| S34 −3 C vs T: Q157R vs annotated | 1.28 | 0.30 |

- **This is the expected U2AF1 Q157 signature.** New 3′ splice sites gained in Q157R cells have G after the AG far more
  often. Sites associated with WT, i.e. used relatively less in the mutant, are depleted of G (enriched for A) at +1.
  This mirrors the published preference of Q157-mutant U2AF1 (Q157P in the paper; G at +1).
- The S34F-type −3 C/T position does not differ between genotypes, as expected for a Q157 mutant.
- All novel A3SS acceptors together have slightly more C at −3 than annotated ones (OR 1.22, p = 4 × 10⁻⁷), the same in
  both genotypes. This is a general property of the novel set, not genotype-related.
- Caveats: the groups are defined by detection pattern and enrichment with 3 mice per genotype, and each test counts
  about 100 acceptors. The effect is large and fits the known biology, but should be confirmed on independent data
  (Phase 7 short reads).
- Logos (−20…+3): `phase5_3ss_logos.{pdf,png}`. The +1 position is G-dominated in the Q157R set and A-dominated in the WT
  set.

## D) Over-representation analysis (GO BP; universe = all testable genes)
- Q157R-associated novel isoforms come from 717 genes; WT-associated from 693.
- **No GO term reaches FDR < 0.05.**
- The top Q157R terms are chromosome segregation (FDR 0.10) and **regulation of RNA splicing** (FDR 0.14; 23 genes). The
  splicing term fits the known autoregulation of splicing factors, but it is suggestive only.
- The top WT terms have FDR ≥ 0.25.
- Full table: `phase5_ora_go_bp.tsv`.

## Files (`results/05_splice/`)
`phase5_genotype_associated_isoforms.tsv` (per isoform: specific/enriched call, shares, Δ, FDR, per-sample counts),
`phase5_event_types.tsv`, `phase5_intron_retention.tsv`, `phase5_3ss_position_counts.tsv`, `phase5_3ss_fisher_tests.tsv`,
`phase5_ora_go_bp.tsv`, `phase5_3ss_logos.{pdf,png}`.
