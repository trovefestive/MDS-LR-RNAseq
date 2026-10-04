# Phase 6 — Differential isoform analysis report (2026-10-04)

**Sign convention:** every effect (logFC, Δ isoform proportion) is **Q157R − WT** (WT = reference level).
Input: the Phase 3 final set (171,224 isoforms; R2 counts, duplicate-collapsed). Outputs are in `results/06_diff/`.

## Jobs
| Step | Script | SLURM ID | Outcome |
|---|---|---|---|
| rMATS-long | `06a_rmats_long.slurm` | 20775877 | COMPLETED (32 min, single-threaded) |
| edgeR DGE/DTE, DRIMSeq DTU, switches, LRP2, Ybx1, dot plots | `06b` → `15_phase6_diff_isoforms.R` | 20775882 | COMPLETED |
| rMATS-long import, direction checks, DRIMSeq agreement, event classes | `06c` → `16_phase6_rmats_summary.R` | 20828597 | COMPLETED |

Problems fixed on the way:
- **20775799:** my GTF filter matched `$9` without a tab separator, so the rMATS-long GTF input was empty. Cancelled and
  fixed before any result was used.
- **20775824:** rMATS-long's `cleanup_threads()` waits only 1 s per worker (`join(1)`) and raises "thread(s) were not
  joined" when one exits slowly. Rerun with `--num-threads 1`.

**Rule-8 checks (all passed):**
- edgeR DTE: every mutant-exclusive isoform that was tested has logFC > 0. U2af1-213 logFC = **+8.86**.
- rMATS-long: group 1 = Q157R and group 2 = WT. The source computes `delta = group_1_avg − group_2_avg`. Δ matches Q157R − WT
  recomputed from the per-sample proportions for **all 115,313 isoforms**. U2af1-213 Δ = **+0.099** (adj p = 1e-52).
- DRIMSeq: Δ is computed directly as mean(Q157R) − mean(WT) proportion.

## 1) Gene and transcript expression (edgeR QL, `~genotype`)
- **Gene level (DGE):** 97 of 10,916 genes at FDR < 0.05, **89 up** and 8 down in Q157R.
  - Top up genes are short, highly expressed ribosomal and OXPHOS transcripts: Rps26 +2.2, Rpl38 +2.6, Rps21 +2.5,
    Rpl37a +2.1, Atp5me +2.5, Romo1 +2.6, Cox6b1 +1.9 (log2FC).
  - Down: e.g. Msmo1 −3.0, Mterf3 −2.3, Cdc14a −1.9.
  - **Caveat:** a coordinated shift in short ribosomal-protein transcripts can also come from library or length-related
    differences between samples. Treat it as an observation to confirm (Phase 7 short reads), not as an established
    biological effect.
- **Transcript level (DTE):** 71 of 28,271 isoforms at FDR < 0.05, 66 up.

## 2) Differential isoform usage (DTU)
| Method | Significant isoforms (FDR < 0.05, \|Δ\| ≥ 0.10) | Genes |
|---|---|---|
| **rMATS-long** (paper's tool; its own summary, which adds a read-count filter) | **866** | **576** |
| DRIMSeq (Dirichlet-multinomial) | 130 | 120 genes at gene-level FDR < 0.05 |

- **Agreement:** across 22,979 isoforms tested by both, Δ correlates at **r = 0.90**.
  - **122 of DRIMSeq's 131** significant isoforms are also significant in rMATS-long, all with the same sign.
  - rMATS-long calls about 6× more. With n = 3, the **122 isoforms found by both form the high-confidence set**.
- **Switches without a gene-level change (YBX1-type):** 109 of the 130 DRIMSeq switches, and 577 of 606 genes with an
  rMATS-long switch (by the looser adj p < 0.05, |Δ| ≥ 0.10 filter), show no significant gene-level change.
  Usage shifts are therefore largely separate from expression changes.
- Top usage switches (rMATS-long, all without gene-level DE):

| Gene | Isoform | WT share | Q157R share | Δ |
|---|---|---|---|---|
| Cd34 | Cd34-201 (ENSMUST00000016638) | 0.56 | 0.78 | **+0.22** |
| Cd34 | Cd34-202 (ENSMUST00000110815) | 0.41 | 0.17 | **−0.23** |
| Set | ENSMUST00000067996 / ENSMUST00000102866 | 0.29 / 0.68 | 0.46 / 0.52 | +0.17 / −0.16 |
| Tmpo | ENSMUST00000072239 | 0.59 | 0.41 | −0.18 |
| Atrx | novel ESPRESSO:chrX:1412:85 / :82 | 0.14 / 0.01 | 0.00 / 0.17 | −0.14 / +0.16 |
| Usp1 | novel ESPRESSO:chr4:1127:16 | 0.27 | 0.42 | +0.15 |
| Dock10 | novel ESPRESSO:chr1:908:45 | 0.03 | 0.26 | +0.23 |

## 3) What kind of isoform differences? (rMATS-long classification, top 100 significant genes)
| Class | Events |
|---|---|
| Complex (several changes at once) | 40 |
| Exon skipping (SE) | 35 |
| Alternative 3′SS (A3SS) | 13 |
| Intron retention (RI) | 10 |
| Alternative first exon (AFE) | 9 |
| Alternative ends | 4 |
| Alternative last exon (ALE) | 2 |
| Alternative 5′SS (A5SS) | 1 |

- Exon skipping and A3SS dominate the simple events. That fits a 3′-splice-site factor (U2AF1) and the Phase 5 result.
  There were no mutually exclusive exons.
- **Cd34 (top switch):** Cd34-201 and Cd34-202 differ only by an **A3SS in the last exon**: two acceptors 156 nt apart
  (exon starts chr1:194,642,088 and 194,642,244). Q157R shifts usage to the downstream acceptor (Cd34-201).
  - Both acceptors have **G at +1** (TAG|GAG vs CAG|GGT), so the Q157 +1 G preference does not explain this switch on its
    own. They differ at −3 (T vs C).
  - Cd34 is the classic LK/HSPC marker. An isoform shift with unchanged total Cd34 is worth checking at the protein level,
    but its mechanism is not established here.

## 4) Comparison with the earlier LRP2 run (effects negated to Q157R − WT)
LRP2 outputs named `A_vs_B` use A as the reference. The old DTU `delta_proportion` also correlates negatively with its own
(prop_Q157R − prop_WT), which independently confirms the convention.

| Level | Shared features | r (Δ or logFC) | Significant here | …also significant in LRP2 | Same sign |
|---|---|---|---|---|---|
| Gene (DGE) | 9,900 | 0.95 | 94 | 92 | 100% |
| Isoform (DTE) | 19,700 | 0.95 | 60 | 60 | 100% |
| Isoform usage (DTU) | 13,061 | 0.76 | 102 | 77 | 100% |

- Old `U2af1::sedfdca6…` = **U2af1-213** (exact intron chain): logFC **+8.86** here vs **+8.80** in LRP2 after negation.
- Isoforms are matched by exact intron chain (gffcompare `=`); genes by Ensembl ID without version.
- The two pipelines agree closely once the direction label is read correctly. Their sign disagreement was entirely the
  naming convention.

## 5) Paper-reported event: YBX1 exon skipping
- Ybx1 is dominated by its canonical isoform (ENSMUST00000079644: **92% in both genotypes**). No other isoform exceeds 5%,
  so DRIMSeq's filter (≥2 isoforms at ≥5%) leaves the gene out.
- The largest Ybx1 isoform change is under 0.5 percentage points. The novel exon-skipping-type isoforms are at 0.2–2.3% in
  both genotypes.
- **No YBX1 exon-skipping shift in Q157R mouse LK cells.** In the paper this event was strongest in SRSF2-mutant samples
  and weaker in U2AF1-mutant ones, so the null is plausible rather than contradictory.

## Caveats
- n = 3 vs 3. rMATS-long is more liberal than DRIMSeq, so prefer the 122 isoforms found by both.
- A310 (WT) has a higher promyelocyte/granulocyte signature (Phase 2). Usage or expression differences in granule genes
  may partly reflect sorting composition.
- No orthogonal short-read confirmation yet (Phase 7).

## Figures and tables (`results/06_diff/`)
- `phase6_DGE_edgeR_Q157R_vs_WT.tsv`, `phase6_DTE_edgeR_Q157R_vs_WT.tsv`, `phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv`
- `phase6_isoform_switches.tsv`, `phase6_rmatslong_significant.tsv`, `phase6_rmatslong_vs_DRIMSeq.tsv`
- `phase6_rmatslong_isoform_differences.tsv`, `phase6_rmatslong_event_counts.tsv`, `phase6_vs_LRP2.tsv`,
  `phase6_Ybx1_isoforms.tsv`
- Per-isoform CPM dot plots (every replicate): U2af1, Ybx1, Atrx, Cd34, Rpl37, Atn1, Mrpl33, Lsm7
  (`phase6_dotplot_*.{pdf,png}`)
- rMATS-long Fig 6-style structure and abundance plots for the top 100 genes: `rmats_long/results_by_gene/<gene>/`
  (local; large)
