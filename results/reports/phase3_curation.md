# Phase 3 — SQANTI3 curation, filtering and requantification report (2026-10-02)

## Jobs (final runs)
| Step | SLURM ID | Outcome |
|---|---|---|
| 03a SQANTI3 QC (ORF-ID patch) | 20719927 | COMPLETED, 16 min |
| 03b filters + keep list (P3_FILTER=rules) | 20719930 | COMPLETED; keep list identical to the earlier run |
| 03c requant ×21 chromosomes (R2 → FSM filter → fsm) | 20712722_[1-21] | COMPLETED (chr11 3 h 21 min). Reused, because the keep list is unchanged |
| 03b → 03c → 03d rerun (duplicate-row fix, see below) | 20766249 / 20766250 / 20775374 | COMPLETED; final numbers below are from this run |
| 03b2 ML rerun without ORF/NMD features (comparison) | 20699349 | COMPLETED |

Failed or superseded runs and their fixes:
- **20695215:** TD2 imports `pkg_resources`, which setuptools ≥81 removed. Fixed by pinning setuptools 80.9.0.
- **20696734:** SQANTI3 crashed on a zero-length CDS interval when a CDS ends on an exon boundary. Fixed with local
  patch 1.
- **20697921:** cancelled. The default ML filter kept only 337 new isoforms.
- **ORF-ID bug:** SQANTI3's TD2 header regex `[^\s:]+` cannot take the colons in ESPRESSO IDs, so every new isoform came
  out `non_coding` with NMD = NA. Fixed with local patch 2. QC, filter and summary were rerun; the rules filter does
  not use ORF features, so the keep list did not change.
- Both patches are recorded in `results/metadata/software_versions.txt`, and the originals are kept as `*.orig`.

## SQANTI3 QC (253,272 ESPRESSO isoforms)
- Inputs: GENCODE vM39, refTSS CAGE peaks, PolyASite lifted to mm39, polyA motifs, per-sample read counts (`-fl`),
  and `--include_ORF` (TD2).
- Categories: 126,229 FSM, 83,288 NIC, 39,265 NNC, 1,409 fusion, 917 antisense, 815 ISM, 585 intergenic, 408 genic
  intron, 356 genic.
- ISM is rare because ESPRESSO's read correction keeps truncated reads from becoming separate isoforms. The old LRP2
  set had about 44,600 ISM.
- Coding: 108,187 GENCODE and 114,762 new isoforms. Predicted NMD: 15,185 GENCODE and 40,135 new.

## Filter choice — deviation from the paper (rules filter instead of ML)
| Filter | New isoforms kept | Notes |
|---|---|---|
| ML, default (paper setting) | 337 / 126,955 (0.3%) | top features `predicted_NMD`, `coding`; also hit by the ORF-ID bug (all new isoforms "non_coding") |
| ML, without ORF/NMD and per-sample count features | 26,792 (21%) | top features `bite`, `length`, CAGE distance, A-content after the TTS, polyA distance |
| **Rules (default JSON) — primary** | **86,118 (68%)** | removes intra-priming (>59% A after the TTS), RT-switching, and non-canonical junctions without coverage |

**Why rules:**
- The default ML run learned "coding and not NMD = real". Its positive training set is reference FSMs, which are
  mostly coding. That is circular for an NMD analysis (Phase 4), and the ORF-ID bug made it worse.
- With ORF features removed, ML and rules agree on 22,703 kept and 36,748 dropped new isoforms. The 63,415 new isoforms
  that **only ML drops** have independent support as good as the agreed-kept set (76% vs 66% TSS in a CAGE peak; 86%
  vs 86% TTS at a PolyASite site). They are just shorter (median 2,461 vs 3,066 bp), so ML is learning "length like a
  reference transcript", not "artifact".
- The **agreed-dropped** isoforms are clearly worse at the 3′ end (32% at a PolyASite site), which is the intra-priming
  signature. The rules filter targets exactly that.
- ML also flagged GENCODE transcripts with *more* CAGE support (63%) than the ones it kept (38%).
- GENCODE transcripts the filter flags as artifacts (10,266 under rules) are **restored**, as in the paper.
- The setting is `P3_FILTER` in `config.sh`. Switching to `ml_noORF` means rerunning 03b–03d.

## Requantification and the paper's FSM-read filter
- Keep list: 212,435 isoforms (86,118 new + 126,317 GENCODE).
- `ESPRESSO_Q --read_ratio_cutoff 2` on the filtered GTF gives the R2 matrix. Keeping only reads that are FSM to a
  retained isoform (paper `remove_nonfsm_data.pl`, fixed copy in `scripts/paper_derived/`) and requantifying gives the
  fsm matrix.

| Sample | Genotype | R2 reads | FSM-only reads | % lost |
|---|---|---|---|---|
| V335 | WT | 8,309,834 | 6,821,068 | 17.9 |
| V334 | WT | 4,545,628 | 3,629,603 | **20.2** |
| A310 | WT | 7,639,978 | 6,487,300 | 15.1 |
| X504 | Q157R | 5,924,664 | 4,999,063 | 15.6 |
| A258 | Q157R | 8,023,090 | 6,580,202 | 18.0 |
| A309 | Q157R | 6,610,173 | 5,325,809 | 19.4 |

- With PacBio, read loss is 15–20%, against about 47% for the paper's ONT data.
- **By the pre-set rule (fsm primary only if every sample loses <20%), R2 is primary.** V334 just exceeds the cutoff.
  The fsm matrix is kept as a sensitivity check (`results/03_sqanti/quant/fsm_abundance.esp`). There is no genotype
  pattern in the loss (WT 15.1–20.2%, Q157R 15.6–19.4%).

## Final isoform set (primary R2 matrix; `results/03_sqanti/final/`)
- Rule: all GENCODE transcripts, plus new isoforms with ≥1 read in ≥2 samples and ≥5 reads in total.
- **171,224 isoforms: 126,317 GENCODE + 44,907 new.**
  - 43,347 new isoforms (96.5%) are also seen in ≥2 of the 3 mice of one genotype.
- New isoforms by category: 32,736 NIC, 11,076 NNC, 490 fusion, 217 ISM, 188 antisense, 119 intergenic, 55 genic,
  15 FSM, 11 genic intron.
- New isoforms carry **4.2–4.7% of TMM-normalized expression in every sample** (WT 4.21–4.57%, Q157R 4.26–4.65%), so
  there is no global genotype difference.
- 2,926 new isoforms account for more than 10% of their gene's expression (mean CPM).
- **Coding and NMD:** 41,900 of 44,907 new isoforms are coding, and 15,072 are predicted NMD.
  - Among coding isoforms, **36.0% of new ones are predicted NMD vs 14.0% of GENCODE**, the same direction as the
    paper. The formal test is in Phase 4.
- U2af1-213 (the cis-created Q157R isoform) is in the final set: 0/0/0 reads in WT, and 33/64/79 in X504/A258/A309.
- CPM is edgeR TMM, as in the paper. Files: `final_counts.tsv.gz`, `final_cpm_tmm.tsv.gz`,
  `final_isoform_annotation.tsv.gz`, `final_isoform_ids.txt`. Column order is samples.tsv order: V335, V334, A310,
  X504, A258, A309.

## Cross-checks
- **IsoQuant (gffcompare, exact intron chain).** IsoQuant reports fewer isoforms per gene: 86k models, about half of
  them single-exon.
  - From IsoQuant's side, our final set recovers **93.2% of IsoQuant's known multi-exon models (26,413/28,337)** and
    **56.5% of its new multi-exon models (9,658/17,097)**.
  - From our side, only 6.4% of our new isoforms match an IsoQuant model exactly. Many IsoQuant "new" models match
    transcripts that are GENCODE vM39 in our set. Most of our unmatched new isoforms share junctions with an IsoQuant
    model (class `j`) but differ in the full chain.
- **Old LRP2 filtered set:** **81.2% of our new isoforms** (36,450 of 44,907) have an exact intron-chain match in the April LRP2
  transcriptome, so new structures replicate across pipelines. (49% of our GENCODE set matches, because LRP2 kept only
  expressed models from an older annotation.)

## Correction: duplicated rows from ESPRESSO_Q (2026-10-02)
- In the requant step (`ESPRESSO_Q --read_ratio_cutoff 2`), **390 multi-exon isoforms that have no GENCODE gene were
  written twice**, with their reads split across the two records. All of them had `gene_id "NA"` in the ESPRESSO GTF.
- Giving them their ESPRESSO cluster ID as gene ID (03b) did not stop it.
- Summing the two copies reproduces the pre-filter (Phase 2) counts for 336/390 isoforms, against 16/390 for taking the
  maximum. So 03d now **sums duplicate rows** in both matrices and drops exact-duplicate GTF lines.
- `12_phase3_summary.R` stops if any duplicate ID remains.
- Effect on the final set: 171,213 → 171,224 isoforms (novel 44,896 → 44,907). Every other conclusion is unchanged.

## Caveats
- With n = 3 vs 3 and V334 at about 55% of the depth of the others, low-count new isoforms near the ≥5-read threshold
  are the least reliable. The within-genotype flag is provided for that reason.
- There is no matched short-read data, so 5′/3′ support comes from CAGE and PolyASite atlases, not from these samples.
- SQANTI3 6.0.2 carries two local patches (see above). Neither changes how isoforms are classified structurally.

## Figures (`results/03_sqanti/`, PDF + PNG)
`phase3_novel_categories`, `phase3_length_known_vs_novel`, `phase3_expression_known_vs_novel`, `phase3_fsm_read_loss`.

## Ready for Phase 4
ORF/NMD characterization on the final set (known vs novel NMD, per-sample NMD expression share in Q157R vs WT).
