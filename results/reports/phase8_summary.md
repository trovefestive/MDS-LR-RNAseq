# Phase 8 — Summary and follow-up candidates (2026-10-04)

**Sign convention:** every Δ (isoform proportion) and logFC is **Q157R − WT** (WT = reference level).
Script: `08a_summary.slurm` → `18_phase8_candidates.R` (job 20843265). Outputs: `results/08_summary/`.

## Candidate list
**Candidates** = isoforms significant in rMATS-long (adj p < 0.05, |Δ| ≥ 0.10; 909) or DRIMSeq (FDR < 0.05, |Δ| > 0.10;
131): **918 isoforms in 611 genes**, 122 significant in both.
- The rMATS-long count here (909) is the plain adj p / |Δ| filter. The 866 in the Phase 6 report is rMATS-long's own summary, which adds a read-count filter.

### Evidence per isoform
| Evidence | Definition |
|---|---|
| Both methods | rMATS-long and DRIMSeq both significant, same sign |
| LR replicates separate | all 3 Q157R long-read proportions above (or below) all 3 WT |
| LR ≥ 10 reads | ≥ 10 long reads per sample on average in the genotype where the isoform is higher |
| IsoQuant agrees | IsoQuant built the identical intron chain (gffcompare `=`) and its Δ has the same sign; n/a without such a model |
| SR same sign | short reads (separate cohort, **LSK** cells, Salmon): same sign, \|Δ\| ≥ 0.05, gene ≥ 10 reads |
| SR replicates separate | all 3 short-read Q157R proportions above (or below) all 3 WT |
| Predicted NMD, gene-level DE | reported, not scored |

### Tier definitions
| Tier | Rule | Isoforms | Genes (best isoform per gene) |
|---|---|---|---|
| **A** | both methods + LR replicates separate + SR same sign + IsoQuant not contradicting | 41 | **36** |
| **B** | both methods + (SR same sign or IsoQuant agrees) | 48 | 35 |
| C | everything else | 829 | 540 |

- **IsoQuant:** 484 candidates have an identical IsoQuant model; for 411 of 475 with a usable Δ (87%), the IsoQuant Δ
  has the same sign.
- **Short reads:** 907 candidates are testable. Across all of them the agreement is weak (r = 0.18, same sign 57%), and much
  stronger in the high-confidence set (Phase 7: 69% same sign, r = 0.44). So rMATS-long-only calls should be treated
  as weak.
- **NMD:** none of the tier A or B isoforms is predicted NMD. 39 tier C isoforms are, 17 of them up in Q157R.
- **Gene-level change:** 5 tier A genes also change at gene level (Cdc14a, Msmo1, Atp5if1, Ndufa7, Gng5). The other 31
  are switches without a gene-level change.

## Tier A genes (`phase8_candidates_genes.tsv`)
Score = number of the 6 scored evidence types met. Δ = long-read Δ of the listed isoform; SR Δ = short-read Δ.

| Score | Gene | Isoform | Long-read Δ | Short-read Δ | Event with partner isoform (rMATS-long) |
|---|---|---|---|---|---|
| 6 | Cdc14a | ENSMUST00000090464 | −0.42 | −0.20 | complex (gene also down) |
| 6 | Mtm1 | ENSMUST00000061970 | −0.19 | −0.37 | |
| 6 | Erbin | ENSMUST00000188997 | +0.14 | +0.22 | |
| 6 | Crlf3 | ENSMUST00000061283 | −0.11 | −0.19 | exon skipping |
| 5 | Atn1 | ENSMUST00000507053 | +0.74 | +0.08 | A3SS |
| 5 | Gm16758 | novel ESPRESSO:chr17:533:12 | +0.59 | +0.12 | |
| 5 | Thada | ENSMUST00000047524 | +0.31 | +0.09 | A3SS |
| 5 | Ppp3r1 | ENSMUST00000102880 | −0.29 | −0.10 | |
| 5 | Rbl1 | ENSMUST00000029170 | −0.28 | −0.06 | complex |
| 5 | Cd34 | Cd34-202 ENSMUST00000110815 | −0.23 | −0.11 | **A3SS, confirmed by rMATS-turbo (ΔPSI −0.16)** |
| 5 | Tcf19 | ENSMUST00000160885 | +0.22 | +0.20 | A3SS |
| 5 | Tmpo | ENSMUST00000072239 | −0.18 | −0.11 | complex |
| 5 | Mpzl1, Emp3, Efcab14, Cdc37l1, Idh3g, Msmo1, Snhg8, 2410006H16Rik | | | | |
| 4–5 | OXPHOS / mito-ribosome: Atp5if1, Atp5mg, Cox16, Ndufv3, Ndufa7, Ndufa13 (novel), Ndufb8, Mrpl23, Mrpl33, Rpl22l1 | | | | exon skipping, alternative last exon, complex |
| 4 | Rbm34, Surf2, Gng5, Vps29, Pex2, Atp11c | | | | |

- 17 of the 36 tier A genes are among the top-100 genes whose isoform pairs rMATS-long classified. Their events are
  exon skipping (Crlf3, Idh3g, Ndufa13, Ndufv3, Pex2, Vps29), A3SS (Atn1, Cd34, Tcf19, Thada), alternative last exon
  (Cdc37l1, Ndufb8) and complex (Atp5if1, Cdc14a, Mrpl33, Rbl1, Tmpo). Tier A is real splicing, not just alternative ends.
- 34 of the 36 are GENCODE (full-splice match); two are novel (Gm16758, Ndufa13).

## Suggested wet-lab follow-up
1. **Cd34 A3SS** (last exon; acceptors 156 nt apart): the strongest single event. It is significant in rMATS-long,
   DRIMSeq and short-read rMATS-turbo, in two cohorts and two cell populations. RT-PCR across the last-exon acceptor
   region (two products 156 bp apart); then protein-level check of the two Cd34 C-termini.
2. **Score-6 switches** (Cdc14a, Mtm1, Erbin, Crlf3): every scored evidence type met; isoform-specific RT-qPCR.
3. **A3SS events** (Atn1, Tcf19, Thada): the event class expected for U2AF1. Check whether the Q157R-favoured acceptor
   has +1 G (Phase 5 signature).
4. **U2af1-213** (cis, created by the Q157R base): a built-in positive control for any assay.

## Caveats
- n = 3 vs 3 in both datasets; the short-read cohort is a different population (LSK vs LK) and different animals.
- Salmon isoform proportions can differ in absolute level from long-read proportions when isoforms share most of their
  sequence (e.g. Crlf3: long reads 0.91 → 0.80, short reads 0.18 → 0.00). Direction is the comparable quantity.
- Several tier A genes are short OXPHOS / ribosomal genes. Their usage shifts replicate in short reads, but their
  long-read gene-level increases did not (Phase 7), so read them as isoform-level observations.
- A310 (WT) has a higher promyelocyte signature (Phase 2).

## Outputs (`results/08_summary/`)
- `phase8_candidates_isoforms.tsv`: all 918 candidate isoforms with every evidence column
- `phase8_candidates_genes.tsv`: one row per gene (best isoform + strongest opposite-direction partner)
- `phase8_tier_counts.tsv`
- `phase8_lr_vs_sr_delta.{pdf,png}`: long- vs short-read Δ for all candidates, tier A labelled
- `phase8_evidence_top_genes.{pdf,png}`: evidence grid for the top 30 tier A/B genes
- `phase8_tierA_proportions.{pdf,png}`: every replicate, long and short reads, for the top 12 tier A switches
