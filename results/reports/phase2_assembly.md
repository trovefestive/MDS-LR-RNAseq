# Phase 2 — Transcript assembly and quantification report (2026-09-30)

## Jobs
| Job | SLURM ID | Outcome |
|---|---|---|
| 02a ESPRESSO per chromosome, attempt 1 | 20679625 | FAILED: bioconda `ESPRESSO_*.pl` have no shebang (exit 2 before S). Fixed by calling them via `perl`. |
| 02a attempt 2 | 20679746 | chrY COMPLETED. chr3/chr15 FAILED inside ESPRESSO_C; the rest were cancelled. |
| 02a attempt 3 | 20680397 | chr14 FAILED with the same error; the rest were cancelled. |
| **02a final** | **20680599_[1-20]** + chrY from 20679746 | **All 21 COMPLETED**, 0 retries used |
| 02b merge | 20680600 | COMPLETED (2.5 min) |
| 02c IsoQuant, attempt 1 | 20679627 | FAILED after read assignment: IsoQuant 4.0.0 bug in the `--check_canonical` code path |
| **02c IsoQuant final** | **20683816** | **COMPLETED** (2 h 02 min), run without `--check_canonical` |

**ESPRESSO_C failure: root cause and fix.**
- The error was `Failed to run sort ... Exit code is -1`. Perl's `system()` returning −1 means the fork itself failed; `sort` never ran.
- The compute nodes use strict memory overcommit (`vm.overcommit_memory=2`). Every ESPRESSO_C worker thread forks the
  whole multi-GB Perl process to launch `sort`, and those forks were refused at random when several tasks shared a node.
- My first guess, an oversized `--sort_buffer_size 8G`, was wrong; the same error recurred at the default 2G.
- Final settings: `ESPRESSO_C -T 8` (as in the paper's scripts) and up to 3 retries per sample. No retry was needed in
  the final run.

**Wall time (final array).** chr20–chrY took 40–60 min, most chromosomes 1.5–2.5 h, and chr11 took 7.5 h. chr11 was the
slowest at every step (35–120 min per sample in C, 65 min in Q). Memory stayed within 128 GB.

## Method (mirrors the paper's `assembly/assemble_transcripts.sh`)
1. Per chromosome (chr1–19, X, Y), keep primary alignments only (`samtools view -F 0x900`) from each sample BAM.
   chrM and unplaced scaffolds are excluded, as in the paper.
2. `ESPRESSO_S` runs on all 6 samples together, using the GENCODE vM39 GTF for that chromosome.
3. `ESPRESSO_C` runs once per sample.
4. `ESPRESSO_Q` uses the default `-N 2 -R 0` and writes `*_N2_R0_abundance.esp`, `*_updated.gtf` and the per-read
   compatible-isoform table (needed for the Phase 3 FSM-read filter).
5. Per-chromosome outputs are merged into `results/02_espresso/merged/` (header once, rows concatenated, GTF sorted).
- IsoQuant 4.0.0 runs on the full BAMs: `--data_type pacbio_ccs --fl_data --stranded forward --polya_trimmed all
  --complete_genedb --sqanti_output`.

## Results (`scripts/11_phase2_summary.py`; tables in `results/02_espresso/`)
| Sample | Genotype | FLNC reads | Assigned to isoforms | % assigned | % reads on novel isoforms |
|---|---|---|---|---|---|
| V335 | WT | 10,265,655 | 8,443,063 | 82.2 | 6.5 |
| V334 | WT | 5,679,950 | 4,619,127 | 81.3 | 6.6 |
| A310 | WT | 9,145,749 | 7,759,683 | 84.8 | 6.2 |
| X504 | Q157R | 7,164,196 | 6,028,925 | 84.2 | 6.8 |
| A258 | Q157R | 9,831,483 | 8,156,113 | 83.0 | 6.5 |
| A309 | Q157R | 8,297,963 | 6,714,330 | 80.9 | 6.3 |

- **ESPRESSO, before any filtering:** 253,272 isoforms. 126,317 are GENCODE vM39 transcripts with reads, and 126,955 are
  novel ESPRESSO models.
  - Novel isoforms carry only 6.2–6.8% of reads, with no difference between genotypes.
  - 62,848 novel isoforms already pass "≥2 samples, ≥5 reads". This is a preview only: in the paper the rule is applied
    after SQANTI3 filtering and requantification (Phase 3).
- About 16–19% of reads are not assigned to any isoform. That covers chrM and scaffold reads (excluded), MAPQ 0 reads,
  and reads ESPRESSO could not match to any isoform. Phase 3 breaks this down.
- **IsoQuant:** 88,225 transcript models (28,357 GENCODE, 59,868 novel). It is more conservative than ESPRESSO
  before filtering. Overlap is assessed after the SQANTI3 curation in Phase 3.

## Finding 1 — the "mutant-exclusive U2af1 isoform" arises from the Q157R allele itself (cis), and it is UP in Q157R
The isoform from the original figure (old LRP2 `U2af1::sedfdca6…`, reported there as logFC −8.80, "down in Q157R") has
**the same intron chain as the GENCODE vM39 transcript U2af1-213 (ENSMUST00000468653, TAGENE, protein-coding)**. Only
the transcript ends differ.

| Transcript | A258 | A309 | X504 | A310 | V334 | V335 |
|---|---|---|---|---|---|---|
| U2af1-213 (reads) | 64 | 79 | 33 | 0 | 0 | 0 |
| U2af1-201 canonical (reads) | 372 | 482 | 260 | 344 | 213 | 652 |

- The isoform is **present only in Q157R**, so it is up in the mutants. The old "down" call was the sign inversion.
- U2af1-213 differs from the canonical transcript at a single 5′ splice site. Exon 3 ends at **chr17:31,867,169**, the
  exact base mutated by Q157R, instead of at 31,867,157. The result is a predicted in-frame 12-nt (4-aa) shorter exon,
  to be confirmed with the ORF calls in Phase 4.
- In transcript orientation the WT sequence there is `CA|GTATG`. The Q157R base (A→G at the exon's last position)
  makes it `CG|GTATG`, giving a G at −1 next to a GT donor. That is a much stronger 5′ splice site.
- **The reads confirm an allele-specific (cis) effect**
  (`results/02_espresso/phase2_u2af1_junction_alleles.tsv`):
  - All 224 reads using the U2af1-213 junction in the mutants carry the mutant C (85, 87 and 52 reads). None carry T.
  - Reads using the canonical junction carry both alleles (C 26–34%).
  - WT samples show this junction once in 1,680 reads (1 read in A310, carrying T).
  - So about **26–32% of mutant-allele U2af1 transcripts** use the new donor.
- **Interpretation:** this is a direct cis consequence of the knock-in base change on U2af1's own pre-mRNA. It is not
  evidence of trans splicing dysregulation by mutant U2AF1. **Treat this isoform as a genotype marker, not as a
  downstream target.** It also gives a clean positive control for the rule-8 sign check in Phase 6: it must come out
  with logFC > 0.
- The reads using U2af1-213 still contain the mutated base, so they are already counted in Phase 1's 40–42% mutant
  allele fraction. The shortfall from 50% therefore isn't explained by this isoform. Its cause (allele-specific
  expression or mutant-transcript decay) can't be resolved from these data.

## Finding 2 — the A310 read-length spike is Mpo
The 2.2–2.3 kb spike in A310 (Phase 1) is **Mpo-201**: 35,387 CPM in A310 versus a mean of 15,207 in the others, about
3.5% vs 1.5% of reads. A310 also has higher Elane (2,733 vs 1,213 CPM), with Alas1 and other granule/secretory genes
modestly up (`phase2_A310_overrepresented.tsv`).
- This points to a larger promyelocyte/granulocyte-committed fraction in the A310 LK sort, a **cell-composition
  difference**, not a library artifact.
- A310 is a WT replicate, so this composition effect could look like a "down in Q157R" signal for neutrophil-granule
  genes. Phase 6 should check sample-level clustering and report composition-sensitive genes with that caveat.

## Outputs
- `results/02_espresso/merged/merged_N2_R0_abundance.esp` (19.6 MB) and `merged_N2_R0_updated.gtf` (395 MB). These are
  local; a gzipped copy of the abundance table is in the repo.
- `results/02_espresso/work/<chr>/<chr>_compatible_isoform.tsv`: per-read isoform compatibility, needed in Phase 3.
- `results/02_isoquant/LR/`: IsoQuant transcript models, counts, TPM and the SQANTI-like table (local).
- `results/02_espresso/phase2_*.tsv`: summary tables from `scripts/11_phase2_summary.py`.

## Ready for Phase 3
SQANTI3 QC with the ML filter (plus the rules filter for comparison), restoring GENCODE transcripts, ESPRESSO requant
with `--read_ratio_cutoff 2`, the FSM-read filter with per-sample read loss, and then the novel-isoform rule.
