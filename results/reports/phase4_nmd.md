# Phase 4 — ORF / NMD characterization report (2026-10-02)

**Sign convention:** every difference is **Q157R − WT** (WT = reference level).

## Job
`scripts/04a_nmd.slurm` → `scripts/13_phase4_nmd.R`, job 20775375 (final, on the duplicate-corrected Phase 3 set;
earlier runs 20757910/20757955 gave the same conclusions). Input: Phase 3 final set (171,224 isoforms, R2 matrix, TMM CPM) with SQANTI3/TD2 ORF and `predicted_NMD`
calls. Outputs are in `results/04_nmd/`.

**Rule-8 checks (both passed):**
- All **754 mutant-exclusive isoforms** (0 reads in all WT, >0 in all Q157R, ≥5 reads) get logFC > 0 under the
  `~genotype` design. The U2af1-213 positive control has logFC = **+8.86**.
- Every significant NMD-share call must have matching logFC and Δproportion signs. There were no significant calls to
  check.

## 1) NMD is enriched among novel coding isoforms
| Origin | Coding isoforms | Predicted NMD | % NMD |
|---|---|---|---|
| GENCODE | 108,187 | 15,185 | 14.0 |
| Novel | 41,900 | 15,072 | **36.0** |

- **Fisher's exact test: odds ratio 3.44 (novel vs GENCODE), p < 1e-300.**
- By novel category: NIC **42.1%** (31,006 coding), NNC 18.3% (10,124), fusion 26.8% (456).
  - NIC isoforms are new combinations of known splice sites (for example, exon skipping or inclusion). They often
    shift the reading frame or add a premature stop.
- This matches the paper's central observation that novel isoforms are NMD-enriched. Here it holds in normal and
  mutant mouse LK cells alike.
- `predicted_NMD` is SQANTI3's sequence-based rule (a stop codon >50 nt upstream of the last junction). It is a
  prediction, not a measurement of decay.

## 2) No genotype difference in the expression share from NMD isoforms
| Metric (% of TMM CPM) | WT mean | Q157R mean | Q157R − WT | Welch p |
|---|---|---|---|---|
| All NMD isoforms | 1.91 | 1.88 | −0.04 | 0.77 |
| Novel NMD isoforms | 1.03 | 1.03 | +0.002 | 0.98 |
| Novel share of NMD expression | 53.8% | 54.9% | +1.1 | 0.29 |

- Per sample, NMD isoforms are 1.70–2.00% of expression, and novel ones carry about 54% of it in both genotypes.
- With n = 3 per group this is descriptive. There is no sign of a global NMD-isoform increase in Q157R LK cells. This
  differs from the paper's tumor-vs-normal comparison, which is a different contrast (whole AML/MDS vs sorted normal
  cells).

## 3) Per-gene NMD-isoform share: no significant changes
- Method:
  - For each gene with both NMD and non-NMD coding isoforms (≥10 reads in ≥3 samples), the reads were split into an
    NMD bin and a non-NMD bin.
  - The bins were tested with edgeR `diffSpliceDGE` (`~genotype`).
  - A gene-level edgeR QL test was used to flag genes that are down overall.
- **5,907 genes tested; 0 pass FDR < 0.05 and |Δ share| > 0.10.** The lowest FDR is 0.27.
- The p-values are close to uniform: 9 genes at P < 0.001 (6 expected by chance), 66 at P < 0.01 (59 expected). With
  3 vs 3 there is **no detectable gene-specific NMD shift**.
- 93 genes have |Δ share| > 0.10 (56 up, 37 down), but none are statistically supported.
- Top nominal hits are mostly small changes in ribosomal/translation genes (e.g. Rpl23a −0.07, Rps14 −0.03, Rps19
  −0.02; Bicd1 +0.09, Dbf4 +0.11). Treat these as hypotheses only.
- Full table: `phase4_nmd_share_by_gene.tsv`, with per-sample shares, Δ, FDR and gene-level logFC.

## Caveats
- With n = 3 vs 3, effects below roughly 10–15 percentage points in a gene's NMD share are not detectable.
- Bulk LK cells: the A310 composition difference (higher Mpo/Elane, Phase 2) adds noise between samples.
- Most isoform-specific effects of U2AF1 Q157R are expected at the level of individual splice sites and events rather
  than as a global NMD shift. Those are tested in Phase 5 (3′SS motifs) and Phase 6 (differential isoform usage).

## Figures (`results/04_nmd/`, PDF + PNG)
`phase4_nmd_known_vs_novel`, `phase4_nmd_share_per_sample`, `phase4_nmd_share_volcano`. There is no top-gene panel
because no gene was significant.
