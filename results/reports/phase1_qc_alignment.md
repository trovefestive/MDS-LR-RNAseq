# Phase 1 — QC and alignment report (2026-09-30)

## Jobs
| Job | SLURM ID | State | Wall time | Max RSS |
|---|---|---|---|---|
| 01a_qc_align, task 1 V335 | 20677389_1 | COMPLETED | 53m | 72.5 GB |
| 01a_qc_align, task 2 V334 | 20677389_2 | COMPLETED | 27m | 45.4 GB |
| 01a_qc_align, task 3 A310 | 20677389_3 | COMPLETED | 45m | 60.9 GB |
| 01a_qc_align, task 4 X504 | 20677389_4 | COMPLETED | 38m | 58.3 GB |
| 01a_qc_align, task 5 A258 | 20677389_5 | COMPLETED | 47m | 61.7 GB |
| 01a_qc_align, task 6 A309 | 20677389_6 | COMPLETED | 41m | 69.8 GB |
| 01b_qc_summary | 20677390, rerun 20679301 | COMPLETED | <1m | – |

Resources: 16 CPU / 128 GB / 8 h per task. Peak use was 72.5 GB, so 128 GB left plenty of headroom.

Command: `minimap2 -ax splice:hq -uf --secondary=no --junc-bed gencode.vM39.junctions.bed` (minimap2 2.31, prebuilt
splice:hq index), then `samtools sort` and `samtools index`. Each BAM carries `@RG ID/SM = sample`.
Outputs: `results/01_align/bam/<sample>.sorted.bam`, 27 GB in total.

## QC table (`results/01_align/phase1_qc_table.tsv`)
| Sample | Genotype | FLNC reads | Median len | N50 | Mapped | % mapped | Supplementary | Error rate |
|---|---|---|---|---|---|---|---|---|
| V335 | WT | 10,265,655 | 1,929 | 2,229 | 10,264,442 | 99.99 | 20,041 | 0.21% |
| V334 | WT | 5,679,950 | 1,967 | 2,260 | 5,679,352 | 99.99 | 12,167 | 0.22% |
| A310 | WT | 9,145,749 | 2,016 | 2,265 | 9,144,619 | 99.99 | 18,597 | 0.22% |
| X504 | Q157R | 7,164,196 | 2,139 | 2,483 | 7,160,012 | 99.94 | 15,259 | 0.22% |
| A258 | Q157R | 9,831,483 | 1,984 | 2,262 | 9,827,798 | 99.96 | 20,407 | 0.21% |
| A309 | Q157R | 8,297,963 | 1,885 | 2,197 | 8,295,755 | 99.97 | 20,878 | 0.21% |

- Mapping is ≥99.94% in every sample, and the mismatch rate is 0.21–0.22%, as expected for HiFi reads.
- Supplementary alignments, a proxy for chimeric or fusion reads, are 0.2% of reads per sample.
- Depth: V334 has 5.7 M reads, about 55–60% of the others. Expect lower power for V334-driven calls; library size
  normalization handles the rest.
- Read length: X504 is the longest (median 2,139 bp, N50 2,483), and A309 is the shortest (median 1,885 bp). There is
  no consistent genotype difference in medians: WT 1,929–2,016 bp vs Q157R 1,885–2,139 bp. RIN is not available to
  test whether the length variation tracks RNA quality.

## U2af1 Q157R genotype check (`results/01_align/phase1_genotype_check.tsv`)
Position: chr17:31,867,169 (GRCm39), the middle base of codon 157 in ENSMUST00000014684 (minus strand). The reference
codon is CAG (Gln), and Q157R is CGG, which appears as T>C on the genome's + strand. Filters: MAPQ ≥ 20, all base
qualities.

| Sample | Expected | Depth | T (ref) | C (Q157R) | C fraction | Call |
|---|---|---|---|---|---|---|
| V335 | WT | 877 | 877 | 0 | 0.000 | WT ✓ |
| V334 | WT | 336 | 336 | 0 | 0.000 | WT ✓ |
| A310 | WT | 523 | 522 | 1 | 0.002 | WT ✓ |
| X504 | Q157R | 501 | 298 | 203 | 0.405 | Q157R ✓ |
| A258 | Q157R | 674 | 390 | 284 | 0.421 | Q157R ✓ |
| A309 | Q157R | 730 | 423 | 307 | 0.421 | Q157R ✓ |

**All 6 genotypes match samples.tsv.**
- The mutants carry the allele at 40–42% of U2af1 transcripts, consistent with a heterozygous knock-in. The slight
  shortfall from 50% could come from allele-specific expression or reduced stability of the mutant transcript. We
  can't tell which from these data.
- The single C read in A310 (1 of 523) is at sequencing-error level.

## Issue found and fixed
The read-length histogram step sorted its header line together with the data rows. In A310, A258 and A309, the header
ended up below the <100 bp bin, and the R summary silently dropped those three samples from the length plot.
- Fixes: the header is now written outside `sort` in `01a_qc_align.slurm`. The six existing files were repaired in
  place, and their totals now equal the FLNC counts exactly. `10_phase1_qc_summary.R` now stops if any sample is
  missing or its histogram total differs from the read count. The summary was rerun as job 20679301.
- Alignments, the QC table and the genotype results were not affected.

## Figures (`results/01_align/`, PDF + PNG)
- `phase1_read_length_distribution`: FLNC length distribution per sample, colored by genotype.
- `phase1_depth_mapping`: reads and % mapped per sample.
- `phase1_u2af1_Q157R_genotype`: fraction of C reads at codon 157.

## Read-length plot (corrected, job 20679301)
- All six samples now appear. The WT and Q157R distributions overlap, with a shared mode at about 1.7 kb and shoulders
  at about 1.1 kb and 2.2 kb that are present in every sample. These most likely come from highly expressed LK
  transcripts.
- X504 has a heavier tail from 3 to 5 kb, which explains its higher median and N50.
- A310 has a single-bin spike at 2,200–2,299 bp: 6.9% of its reads, vs 4.9–5.4% in the other samples. That fits one
  transcript being highly abundant in that mouse rather than a library-wide problem. Check which gene it is in
  Phase 2, once reads are assigned to isoforms. This is a flag, not a failure.

## Ready for Phase 2
The six sorted, indexed BAMs are ready for ESPRESSO, which needs sorted BAMs with the NM tag (present), and for
IsoQuant.
