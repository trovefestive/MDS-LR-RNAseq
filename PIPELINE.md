# PIPELINE — how the long-read RNA-seq analysis works, step by step

A handover guide for running, understanding and extending this analysis. It explains *what* each step does and
*why*, which files go in and out, and the traps we hit. Read it together with `README.md` (results) and the phase
reports in `results/reports/`. Nothing here needs the original analyst: every step is a script in `scripts/`.

---

## 1. The question and the data

- **Biology:** U2AF1 is the small subunit of the U2 auxiliary factor that recognises the 3′ splice site (3′SS). The
  Q157R mutation (codon 157, second zinc finger) changes which 3′ splice sites are preferred: the base right after
  the AG (+1) matters. S34F, the other common U2AF1 mutation, affects the base at −3 instead. So the main prediction
  is: **in Q157R cells, gained 3′ splice sites should prefer +1 G; the −3 position should not change** (that is the
  negative control).
- **Design:** 3 wild-type vs 3 U2af1 Q157R knock-in mice, Lin⁻ Kit⁺ (LK) bone-marrow progenitors, PacBio Kinnex
  full-length cDNA (FLNC reads, already demultiplexed and polyA-selected by the core). About 50 M reads, ~2 kb each.
- **Why long reads:** one read = one full transcript, so isoforms are observed directly instead of inferred from
  short fragments. The cost is depth (≈8 M reads per sample vs 60 M pairs for short reads) and n = 3.
- **Reference workflow:** Miller et al. 2026 (bioRxiv) for AML/MDS: ESPRESSO → SQANTI3 → rMATS-long. We follow it
  and note every adaptation for PacBio data.
- **Second dataset:** an older in-house Illumina RNA-seq set (3 WT vs 3 Q157R, **LSK** cells, different mice) is used
  in Phase 7 to test whether findings replicate. Different cells and animals, so it is replication, not validation.

## 2. How the project is organised

```
config.sh              every path and parameter; sourced by every script (edit here, not in scripts)
config.local.sh        untracked: your SLURM account/partition, short-read FASTQ path
samples.tsv            long-read samples; samples_shortread.tsv: short-read samples
scripts/NN_*.slurm     one SLURM job per step, numbered by phase; 1N_*.R / .py: the analysis behind each step
scripts/submit.sh      the only way to submit:  scripts/submit.sh scripts/02a_espresso_by_chr.slurm
data/ref/              genome, GTF, indexes, CAGE/polyA references (built by Phase 0)
raw_files/             the delivered FLNC BAM/FASTQ, read-only; data/flnc/ holds symlinks
results/<step>/        outputs per phase; results/reports/phaseN_*.md: what happened and what it means
results/logs/          SLURM logs; results/metadata/software_versions.txt: every tool version and local patch
```

Conventions every script follows:
- `set -euo pipefail`; `# === section ===` headers; `echo "[$(date)] ..."` progress lines.
- **Skip-if-exists:** each step checks for its output before running, so a failed job can be resubmitted and only
  the missing part reruns. Delete the output if you want a step to recompute.
- `source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"` first; conda envs via `activate_env <name>`.
- **Direction rule (rule 8):** every effect size is **Q157R − WT**. Each statistical script asserts this on a known
  mutant-only isoform (U2af1-213) and stops if the sign is wrong. Read `README.md` → "Comparison" and the Phase 6
  report for why: an earlier pipeline (LRP2) names contrasts `A_vs_B` with A as reference, so its "Q157R_vs_WT"
  values are WT relative to Q157R. The numbers were right; only the label was easy to misread.
- Plots use `scripts/00_plot_theme.R` (WT blue `#3B6FB6`, Q157R red `#C8553D`); saved as PDF and PNG.

Running a step:
```bash
scripts/submit.sh scripts/01a_qc_align.slurm          # submit from the project root
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # afterwards
less results/logs/01a_qc_align_<jobid>_1.log
```
Dependencies: `scripts/submit.sh --dependency=afterok:<id> scripts/02b_espresso_merge.slurm`.

## 3. Phase 0 — setup (`00a` envs, `00b` references, `00c` indexes)

| Env | Tools | Used by |
|---|---|---|
| lr-align | minimap2, samtools, seqkit, gffread, bedtools, liftOver | Phases 1, 7 |
| lr-espresso | ESPRESSO 1.6.0 | Phase 2–3 |
| lr-isoquant | IsoQuant 4.0.0 | Phase 2 |
| lr-sqanti3 | SQANTI3 6.0.2 (bioconda) | Phase 3 |
| lr-rmatslong | rMATS-long 2.1.0 | Phase 6 |
| lr-r | R + edgeR, DRIMSeq, clusterProfiler, BSgenome mm39 | all summaries |
| lr-salmon, rmatsenv-2.7 | Salmon; rMATS-turbo 4.1.2 + STAR | Phase 7 |
| lr-leafcutter | LeafCutter 0.2.9 (R), regtools | Phase 7b |

References: GRCm39 primary assembly + **GENCODE vM39** (keep one release throughout; vM39 is much larger than older
releases, so "novel" counts are not comparable with older annotations). SQANTI3 support files: refTSS v4.1 CAGE
peaks (mm39) and PolyASite 2.0 lifted mm10 → mm39, plus the SQANTI3 polyA-motif list. `00c` builds the minimap2
index and a junction BED from the GTF (`paftools.js gff2bed`), which guides spliced alignment.

Traps: SQANTI3 needed `setuptools<81` and two one-line source patches (zero-length CDS; `:` in ORF IDs), all
recorded in `software_versions.txt`. The login node kills conda solves; build envs inside a job.

## 4. Phase 1 — QC, alignment, genotype (`01a` array of 6, `01b` summary → `10_phase1_qc_summary.R`)

- **In:** `data/flnc/<sample>.flnc.fastq.gz`. **Out:** `results/01_align/bam/<sample>.sorted.bam`, QC tables/plots.
- `seqkit stats` for read counts and lengths. Alignment:
  `minimap2 -ax splice:hq -uf --secondary=no --junc-bed <gencode.bed>`. `splice:hq` = HiFi preset (low error);
  `-uf` = reads are transcript-oriented, so only the forward strand is a candidate transcript strand (true for FLNC,
  not for ONT cDNA); `--junc-bed` nudges junctions onto annotated splice sites.
- **Genotype check:** `samtools mpileup` at chr17:31,867,169 (T>C on the + strand = Q157R on the − strand gene).
  Expect ~40% C in the knock-ins (heterozygous) and ~0 in WT. The summary stops if any sample disagrees. Do this in
  every dataset you ever get; it caught nothing here, but it is the cheapest sample-swap check that exists.
- Expected numbers: ≥99.9% mapped, 0.2% mismatch rate, 5.7–10.3 M reads per sample.

## 5. Phase 2 — transcript assembly (`02a` ESPRESSO per chromosome, `02b` merge, `02c` IsoQuant; `11_phase2_summary.py`)

- **ESPRESSO** builds isoforms from the alignments and quantifies them, using the GTF as prior knowledge. It has
  three steps: `S` (collect splice junctions per sample), `C` (correct junctions per sample), `Q` (assemble and count
  across samples). Following the paper, each chromosome is a separate run (21-task array: chr1–19, X, Y), then the
  per-chromosome abundance tables and GTFs are concatenated.
- **Why per chromosome:** memory and wall time; chromosomes are independent for this purpose.
- Input BAMs are filtered with `samtools view -F 0x904` (no secondary/supplementary/unmapped).
- **Out:** `results/02_espresso/merged_updated.gtf` + `merged_N2_R0_abundance.esp` (isoform × sample counts; novel
  isoforms are named `ESPRESSO:chr:cluster:n`).
- **IsoQuant** runs on the same BAMs (`--data_type pacbio_ccs`, `--fl_data`, stranded) as an independent assembler.
  Isoforms reconstructed by both are higher confidence; the overlap is measured in Phase 3 with gffcompare (`=` class
  code = identical intron chain).
- Traps: the cluster has strict memory overcommit, so ESPRESSO_C's forks fail at high thread counts (`-T 8` with a
  retry loop). IsoQuant's `--check_canonical` crashed and was dropped. ESPRESSO_Q can emit the same isoform twice
  for genes without an annotated ID; the Phase 3 merge sums those rows.
- Sanity result: the "mutant-only" U2af1 isoform is GENCODE **U2af1-213**: every read across its junction carries
  the mutant base, i.e. Q157R creates a 5′ splice site in U2af1 itself. It is the positive control from here on.

## 6. Phase 3 — curation (`03a` SQANTI3 QC, `03b` filter, `03c` requantify, `03d` merge → `12_phase3_summary.R`)

- **SQANTI3 QC** classifies every ESPRESSO isoform against GENCODE: FSM (matches a known isoform), ISM (truncated
  match), NIC (new combination of known junctions), NNC (at least one new splice site), fusion, antisense,
  intergenic, genic. It also adds support evidence: CAGE peak at the 5′ end, polyA site/motif at the 3′ end,
  RT-switching risk, non-canonical junctions, ORF and NMD prediction (via TD2).
- **Filter:** SQANTI3 has an ML filter (random forest trained on FSM vs likely artefacts) and a rules filter
  (explicit thresholds). We ran both. The ML default kept only 337 novel isoforms and favoured coding/NMD features,
  so the **rules filter is primary** (`P3_FILTER=rules` in config.sh; `03b2` is the ML-without-ORF comparison).
  Every GENCODE transcript the filter rejected is restored: the reference is treated as truth, as in the paper.
- **Requantify:** ESPRESSO_Q again on the filtered GTF with `--read_ratio_cutoff 2`, so reads are assigned only to
  retained isoforms (output "R2"). The paper then keeps only reads that are FSM to a retained isoform
  (`scripts/paper_derived/remove_nonfsm_data.pl`) and requantifies once more ("fsm" matrix). That cost them ~47% of
  ONT reads; here it costs 15–20%, so R2 stays primary and fsm is a sensitivity check.
- **Novel-isoform rule:** keep every GENCODE isoform, plus novel isoforms seen in ≥ 2 samples with ≥ 5 reads in total.
  Also flag support in ≥ 2 replicates of the same genotype. Final set: 171,224 isoforms (44,907 novel).
- **Out:** `results/03_sqanti/final/final_counts.tsv.gz` (isoform × sample), `final_cpm_tmm.tsv.gz` (edgeR TMM
  CPM), `final_isoform_annotation.tsv.gz` (category, support, ORF/NMD, per-genotype detection), `final_isoform_ids.txt`.
  These three files are the input to everything downstream.

## 7. Phase 4 — ORF / NMD (`04a` → `13_phase4_nmd.R`)

- NMD prediction comes from SQANTI3 (premature stop > 50 nt upstream of the last exon junction).
- Tests: novel vs GENCODE NMD fraction (Fisher); per-sample NMD share of expression by genotype (Welch t-test, not
  Wilcoxon: with 3 vs 3 Wilcoxon's smallest two-sided p is 0.10, so it can never reach significance); per-gene NMD
  share (edgeR `diffSpliceDGE`).
- Result pattern to expect: novel isoforms are more NMD-prone (36% vs 14%), but the genotype does not change the
  NMD share. This is the first script with a rule-8 assertion (754 mutant-only isoforms, all logFC > 0).

## 8. Phase 5 — mutation-associated isoforms and the 3′SS signature (`05a` → `14_phase5_splice_signatures.R`)

- "Q157R-associated" novel isoforms = detected in ≥ 2 Q157R replicates and no WT (specific), or with a significant
  usage shift (edgeR `diffSpliceDGE`, enriched). Same for WT.
- For each associated isoform, find its novel 3′ splice sites (acceptors absent from GENCODE), take −20…+3 around
  the AG from `BSgenome.Mmusculus.UCSC.mm39`, and run Fisher's exact tests: **+1 G vs A** (Q157 signature) and
  **−3 C vs T** (S34F signature, negative control). Compare Q157R-associated vs WT-associated sites and vs 5,000
  random annotated acceptors. Sequence logos per group; event classes (A3SS, A5SS, SE, RI, MXE); GO ORA.
- Expected: strong +1 G enrichment for Q157R (82.5% vs 44.9%, OR 5.7), nothing at −3. If a future dataset were S34F,
  the two tests would swap roles.

## 9. Phase 6 — differential isoform usage (`06a` rMATS-long; `06b` → `15_phase6_diff_isoforms.R`; `06c` → `16_phase6_rmats_summary.R`)

- **rMATS-long** (paper's tool): per-isoform proportion within gene, Q157R vs WT, from the ESPRESSO counts.
  Group 1 must be Q157R so its Δ = group 1 − group 2 = Q157R − WT. Run with `--num-threads 1` (a thread-join race
  crashes multithreaded runs). It also classifies the isoform pair differences (SE, A3SS, RI, AFE, complex…) and
  draws structure/abundance plots for the top genes.
- **DRIMSeq** (Dirichlet-multinomial DTU) on the same counts, as a stricter second method; **edgeR** at gene level
  (DGE) and isoform level (DTE). Isoforms significant in both rMATS-long and DRIMSeq form the high-confidence set.
- Thresholds: FDR < 0.05 and |Δ proportion| ≥ 0.10; pre-filter ≥ 10 counts in ≥ 3 samples.
- "YBX1-type" findings = isoform switches without a gene-level change; most switches here are of that kind.
- Also compares with the earlier LRP2 outputs after negating their effect sizes (see rule 8); they agree (r ≈ 0.95).
- Caveat to carry forward: gene-level DGE here was dominated by short ribosomal/OXPHOS transcripts; Phase 7 did not
  replicate it, so treat such a signal as possibly technical until confirmed.

## 10. Phase 7 — short-read replication (`07a` setup, `07b` per sample, `07c` rMATS-turbo prep + post, `07d` → `17_phase7_compare.R`)

- **Expanded GTF** = the 171,224 final isoforms (GENCODE + novel). `gffread` makes the transcript FASTA; Salmon
  index with the genome as decoy; STAR index with the expanded GTF junctions.
- Per sample: Salmon (`-l ISR`, stranded dUTP; `--gcBias --seqBias`), STAR (sorted BAM with XS strand attribute),
  and the same genotype pileup as Phase 1.
- **rMATS-turbo** on the STAR BAMs with the expanded GTF and `--novelSS`; `--b1` Q157R, `--b2` WT so
  `IncLevelDifference` = Q157R − WT. The all-BAM run segfaulted; it runs as per-BAM `--task prep` then one
  `--task post` (rMATS's documented route for large runs; use it from the start).
- Comparison: genotype confirmation; % of long-read isoforms detected (≥ 1 read in ≥ 2 samples, ≥ 5 total);
  gene-level logFC agreement; direction agreement of the high-confidence usage switches (binomial test);
  rMATS-turbo A3SS +1 G test within events; Cd34 specifically.
- Traps: Salmon's automatic library-type detection labelled 2 of 6 stranded libraries as unstranded; always fix the
  library type from what you know about the prep. Salmon merges isoforms with identical sequence
  (`salmon_index/duplicate_clusters.tsv`): those cannot be compared individually.

## 11. Phase 7b — LeafCutter, annotation-free (`07e`, `07f` → `19_phase7b_leafcutter.R`)

- Everything above depends on a transcript model. LeafCutter instead clusters introns straight from split reads
  (regtools `junctions extract` → `leafcutter_cluster_regtools.py` → `leafcutter_ds.R`), so it is an independent
  check of the splicing results. Run on both the short-read BAMs and the long-read BAMs.
- Groups file lists WT first: LeafCutter codes the first group as baseline, so `deltapsi` = Q157R − WT. The summary
  checks the column order and the U2af1-213 junction before anything else, and verifies intron coordinates by
  requiring GT…AG at > 90% of sampled introns.
- Signature test, annotation-free: inside significant clusters, take acceptor pairs that share a donor; the acceptor
  with the larger ΔPSI is "favoured in Q157R". +1 G at favoured vs disfavoured acceptors (and −3 C/T control).
- **Trap (important for any long-read use of regtools):** its `RF`/`FR` strand modes assume paired reads. On unpaired
  long reads every strand was mislabelled (0% GT-AG). Fix: add an `XS:A:+/-` tag from the read flag (FLNC reads are
  transcript-oriented, aligned with `-uf`) and run `-s XS`. The build also needed R-4.3-compatible `oompaBase`/
  `oompaData` from the CRAN archive and an updated `Makevars` for rstan 2.32 (both in `07e`).
- With n = 3 per group LeafCutter is below its calibrated range (n ≥ 4); adjusted p-values are approximate.

## 12. Phase 8 — candidate list (`08a` → `18_phase8_candidates.R`)

Candidates = isoforms significant in rMATS-long or DRIMSeq. Each is scored on six kinds of evidence (both methods;
long-read replicates fully separated; ≥ 10 reads; IsoQuant built the same intron chain with the same sign;
short-read same direction; short-read replicates separated) and placed in tier A/B/C (rules in
`results/reports/phase8_summary.md`). Outputs: `results/08_summary/phase8_candidates_{genes,isoforms}.tsv` and
three figures. NMD prediction and gene-level DE are reported but not scored.

## 13. What the analysis found (short version)

1. Q157R creates a splice site in U2af1 itself (U2af1-213): the built-in positive control.
2. 3′ splice sites favoured in Q157R prefer +1 G (OR 5.7 from ESPRESSO junctions; 9.7 / 3.5 from LeafCutter on
   long / short reads); the −3 control is null everywhere.
3. Isoform-usage switches replicate in the independent LSK short-read cohort; gene-level expression changes do not.
4. Cd34 last-exon alternative 3′SS: significant by four methods in two cohorts; top wet-lab candidate.
5. 36 tier A candidate genes; none predicted NMD; most without gene-level change.

## 14. If you want to rerun or extend

- **Rerun one step:** delete its output directory (or the specific file the skip-if-exists check looks for), then
  `scripts/submit.sh` it again. Downstream steps must be rerun too.
- **New long-read samples:** add rows to `samples.tsv`, create the `data/flnc/` symlinks, rerun from Phase 1.
  ESPRESSO Q is cross-sample, so Phase 2 onwards must be recomputed.
- **New short-read cohort:** edit `samples_shortread.tsv`, set `SR_FASTQ_DIR` in `config.local.sh`, rerun Phase 7.
  Check strandedness (`lib_format_counts.json`) and set `SR_SALMON_LIBTYPE` / `SR_RMATS_LIBTYPE` in `config.sh`.
- **Another mutation (e.g. S34F):** Phase 5 and 7b already compute both +1 and −3 tests; swap which is the
  expected signature and which the control in the reports.
- **Different genome/annotation:** change `GENCODE_REL` and the reference URLs in `config.sh`, rerun Phase 0; the
  genotype coordinate (chr17:31,867,169, GRCm39) must be lifted too.
- **Before pushing anything:** the repository is public. No names, computing IDs, allocation names, lab names or
  absolute paths in tracked files; `.gitignore` is a whitelist, so new result folders must be added explicitly.

## 15. Reading order for a newcomer

1. `README.md` (samples, direction rule, key results, script table).
2. This file.
3. `results/reports/phase1_qc_alignment.md` → `phase3_curation.md` → `phase5_splice_signatures.md` →
   `phase6_differential_isoforms.md` → `phase7_shortread.md` → `phase7b_leafcutter.md` → `phase8_summary.md`.
4. The scripts for one phase, side by side with its report: start with `01a_qc_align.slurm` and
   `10_phase1_qc_summary.R`, the simplest pair.
5. Miller et al. 2026 (bioRxiv 10.64898/2026.05.20.726635) and its code (zenodo 10.5281/zenodo.20314532), to see what we
   adapted and why.
