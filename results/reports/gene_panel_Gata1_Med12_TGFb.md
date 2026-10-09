# Gene panel: Gata1, Med12, Pik3r5, Gdf11, Acvr2b, Smad2, Smad3 (2026-10-09)

**Sign convention:** every logFC, Δ proportion and ΔPSI is **Q157R − WT**.
Script: `scripts/20_gene_panel.R` (reads existing Phase 3/6/7/7b outputs; refits short-read gene DE over all genes).
Tables: `results/09_gene_panel/panel_{summary,DGE,isoforms,rmats_turbo,leafcutter}.tsv`.
Long reads = LK cells (n = 3 vs 3); short reads = separate LSK cohort (n = 3 vs 3).

## Bottom line
**None of the seven genes changes significantly at gene level in either dataset, and none has a robust splicing
change.** Two genes are too lowly expressed in LK cells to assess (Gdf11, Acvr2b). The only nominally significant
splicing event with |ΔPSI| ≥ 0.1 (Gata1 retained intron, short reads) is driven by a single replicate.

## Expression
| Gene | LR CPM WT (reps) | LR CPM Q157R (reps) | LR logFC | LR FDR | SR TPM WT | SR TPM Q157R | SR logFC | SR FDR |
|---|---|---|---|---|---|---|---|---|
| Gata1 | 210 (239/268/122) | 177 (144/160/227) | −0.26 | 0.91 | 8.4 | 13.4 | +0.69 | 0.88 |
| Med12 | 26 (23/32/24) | 36 (35/39/33) | +0.46 | 0.62 | 33.2 | 36.6 | +0.20 | 0.88 |
| Smad2 | 40 (41/32/46) | 37 (30/35/47) | −0.09 | 0.96 | 46.9 | 50.9 | +0.07 | 0.96 |
| Smad3 | 11 (9/13/10) | 14 (21/14/9) | +0.41 | 0.84 | 37.0 | 39.0 | +0.14 | 0.91 |
| Pik3r5 | 6.3 (4.7/3.2/10.9) | 5.3 (9.0/3.9/2.9) | −0.32 | 0.93 | 9.8 | 12.3 | +0.35 | 0.95 |
| Acvr2b | 0.2 (0.4/0.0/0.3) | 0.5 (0.0/0.5/1.1) | not tested (too low) | — | 6.4 | 6.1 | −0.03 | 0.99 |
| Gdf11 | 0.0 (0/0/0) | 0.1 (0/0.2/0) | not tested (too low) | — | 1.7 | 1.8 | +0.11 | 0.98 |

- Gata1 is high in LK cells (erythroid/megakaryocyte progenitors are in the LK gate) and much lower in LSK (CPM and TPM are not directly comparable across platforms, but the rank difference is large). The two
  datasets disagree in direction (−0.26 vs +0.69) and neither is close to significant. In short reads one Q157R
  replicate (KI3, 19 TPM vs 7–13 in the others) drives the mean.
- Med12 is the most consistent: higher in Q157R in both datasets (+0.46 / +0.20), every long-read Q157R replicate above
  every WT replicate, but FDR 0.62 / 0.88. With n = 3 this is a trend at most.
- Gdf11 is essentially absent from LK long reads and very low in LSK; Acvr2b is barely detected in long reads.
  Absence of evidence here is not evidence of no change.

## Splicing
**Long reads (rMATS-long, DRIMSeq, LeafCutter):** no isoform-usage change in any gene (all adj p ≥ 0.22, all
|Δ| < 0.02). Gata1 is 95% canonical isoform (ENSMUST00000033502) in both genotypes. Smad3 and Pik3r5 have one
dominant isoform (≥ 95%) and are not testable by rMATS-long. LeafCutter long-read clusters (Smad2, Med12, Gata1): all
adj p ≥ 0.58, |ΔPSI| ≤ 0.012.

**Short reads (rMATS-turbo, LeafCutter):**

| Gene | Event | ΔPSI | FDR | Comment |
|---|---|---|---|---|
| Gata1 | RI chrX:7,828,163–7,829,480 | +0.36 | 0.016 | PSI Q157R 0.07 / 0.00 / **1.00**, WT 0 / 0 / 0; ≤ 31 reads per sample. Driven by KI3 alone; not supported by long reads (LeafCutter adj p 0.58) |
| Gata1 | A3SS chrX:7,828,410–7,829,071 | +0.33 | 0.12 | same sample, same reads |
| Smad2 | SE chr18:76,419,931–76,420,020 | +0.06 | 1e-9 | below the 0.1 threshold; PSI 0.94–1.00 vs 0.90–0.94. LeafCutter short reads: the flanking intron, ΔPSI +0.06, adj p 0.061. Long reads: no change |
| Med12 | A3SS chrX:100,337,399–476 | +0.07 | 0.001 | below threshold; LeafCutter n.s. in both datasets |
| Acvr2b | A3SS chr9:119,257,033–142 | +0.06 | 2e-4 | below threshold; driven by 14 skipping reads in one WT sample |

- Smad2 is the only gene with a small, method-consistent signal in short reads (slightly more inclusion of a 90-nt
  exon in Q157R; rMATS-turbo and LeafCutter point the same way), but it is small (≈ 6 percentage points), absent in
  long reads, and below the project's |ΔPSI| ≥ 0.1 threshold.
- None of these genes is among the 918 Phase 8 candidate isoforms.

## Caveats
- n = 3 vs 3 per dataset; low-expressed genes (Gdf11, Acvr2b, Pik3r5) have little power.
- LK vs LSK: Gata1 and the TGF-β/activin receptor components are expressed at very different levels in the two
  populations, so the two datasets are not expected to agree on expression.
- Salmon isoform shares for genes whose isoforms share most of their sequence (e.g. Smad2 canonical 0.38 → 0.16 in
  short reads) are EM estimates and are not confirmed by any event-level test; they are not interpreted here.
