# Phase 0 — Setup report (2026-09-30)

## Jobs
| Job | SLURM ID | State | Wall time |
|---|---|---|---|
| 00a_setup_envs | 20676343 | COMPLETED (0:0) | 18m51s |
| 00b_download_refs | 20676344 | COMPLETED (0:0) | 10m01s |
| 00c_build_indexes | 20676345 | COMPLETED (0:0) | 1m33s |

## Environments (full list: `results/metadata/software_versions.txt`)
| Env | Key tools |
|---|---|
| lr-align | minimap2 2.31, samtools 1.24, seqkit 2.14.0, pbmm2 26.2.99, NanoPlot 1.48.0, bedtools 2.31.1, gffread 0.12.9, liftOver 482 |
| lr-espresso | ESPRESSO 1.6.0 |
| lr-isoquant | IsoQuant 4.0.0 |
| lr-sqanti3 | SQANTI3 6.0.2 (bioconda; R 4.5.2) |
| lr-rmatslong | rMATS-long 2.1.0 (R 4.5.3, DRIMSeq 1.38, edgeR 4.8.2) |
| lr-r | R 4.4.3, edgeR 4.4.0, limma 3.62.1, DESeq2 1.46.0, DRIMSeq 1.34.0, DEXSeq 1.52.0, IsoformSwitchAnalyzeR 2.6.0, clusterProfiler 4.14.0, BSgenome.Mmusculus.UCSC.mm39 1.4.3, org.Mm.eg.db 3.20.0 |

Command names in these versions, which differ from older docs and the paper's scripts:
- ESPRESSO: `ESPRESSO_S.pl`, `ESPRESSO_C.pl` and `ESPRESSO_Q.pl`, all on PATH.
- IsoQuant 4.0: `isoquant`. There is no `isoquant.py`.
- SQANTI3: `sqanti3_qc.py`, `sqanti3_filter.py` and `sqanti3`.
- rMATS-long: `rmats-long`. The paper calls `python rmats_long.py` from inside the repo.

SQANTI3 comes from bioconda, not the GitHub clone plus env yml that the plan named. It is the same release (v6.0.2), and
it was simpler to install.

## References (`data/ref/`, md5 in `results/metadata/reference_provenance.txt`)
| File | Content |
|---|---|
| GRCm39.primary_assembly.genome.fa (+ .fai, chrom.sizes) | 61 contigs |
| gencode.vM39.primary_assembly.annotation.gtf | 78,348 genes; 481,956 transcripts |
| gencode.vM39.junctions.bed | BED12 for minimap2 `--junc-bed` (one line per transcript) |
| GRCm39.primary_assembly.splice.mmi | minimap2 index, splice:hq preset, 11 GB |
| refTSS_v4.1_mouse_coordinate.mm39.noheader.bed | 172,324 CAGE peaks (header line removed; the original is kept) |
| atlas.clusters.2.0.GRCm38.96.lifted_mm39.bed | PolyASite 2.0: 300,758 of 301,006 clusters lifted mm10 → mm39 (248 unlifted, 0.08%) |
| mouse_and_human.polyA_motif.txt | 16 standard polyA hexamers (SQANTI3 v6 no longer bundles this list) |

Sanity check: U2af1 is at chr17:31,866,036–31,878,236 (−) in vM39, which matches the old LRP2 GRCm39 coordinates.

**Note: GENCODE vM39 is much larger than older mouse releases.** It has 481,956 transcripts: 180k protein_coding,
152k lncRNA, 89k NMD and 22k retained_intron. This reflects the large long-read-based additions in recent
Ensembl/GENCODE releases.
- Consequence: many isoforms that older annotations would call "novel" will now be FSM or ISM.
- So our known/novel split is **not directly comparable** to the paper (Ensembl 95) or to the old LRP2 run (older
  vM release). Each comparison will state the annotation used.
- Downstream steps must also handle ~2.7× more reference transcripts. That matters for ESPRESSO memory and for the
  rMATS-long `--gencode-gtf`.

## Paper code (data/paper_code/chrisamiller-longread-8ff254a) — parameters to reuse
Read from `assembly/assemble_transcripts.sh` and the `espresso_step*.sh` scripts:
1. **Split BAMs by chromosome.** They keep primary alignments only (`samtools view -F 0x900`) on the main chromosomes
   (1–22/X/Y in human; chrM and scaffolds are dropped). ESPRESSO then runs **per chromosome** with all samples
   together.
2. **ESPRESSO steps and resources:**
   - S: `-T 6`, 96 GB.
   - C: per chromosome × sample, `-T 8`, 72 GB.
   - Q: per chromosome, 1 thread, 300 GB. Human scale with 71 samples; ours is 6 samples, so much smaller.
3. **SQANTI3 v5.2.1 QC** uses short-read SJ.out.tab and BAMs. We have no matched short reads, so we use CAGE, polyA
   sites and motifs instead.
4. **SQANTI3 `filter ml`** with default settings. Their filter is the ML filter, not the rules filter.
5. **Restore reference transcripts:** any Ensembl transcript flagged "Artifact" goes back into the filtered GTF.
6. **Requant:** `ESPRESSO_Q --read_ratio_cutoff 2` on the filtered GTF, so only GTF transcripts are quantified.
7. **FSM-read filter:** reads not classified as FSM to a retained transcript (`remove_nonfsm_data.pl`) are removed,
   then Q runs again. This is where the paper's ~47% read loss comes from.
8. **rMATS-long:** `--delta-proportion 0.1 --adj-pvalue 0.1`. Note that adjusted p is 0.1 here, not 0.05. The groups
   are g1/g2 from the config.
9. **CPM:** edgeR TMM CPM from the final abundance table (`espresso_to_cpm.R`).
10. The "novel isoform in ≥2 samples, ≥5 reads" rule is not in the deposited assembly code. We apply it after the final
    requant, as the plan says.

### Implications for our plan
- **SQANTI3 filter:** the plan currently says the rules filter. The paper used the ML filter. Proposal: run ML as the
  primary filter to match the paper, and rules as a comparison, then report how much the two agree.
- **FSM-read filter:** the paper's step removes 35–65% of ONT reads. Our FLNC reads are already full-length with a
  required polyA tail. Proposal: run it and report the read loss per sample, but make it the primary quantification
  only if the loss is modest (<20%). Otherwise keep the step-6 requant as primary and use the FSM-only table as a
  sensitivity check.
- **rMATS-long direction (rule 8):** before Phase 6, confirm which group is subtracted from which in rMATS-long's
  delta proportion. Test it on a known mutant-only isoform, and set group-1/group-2 so Δ = Q157R − WT.
- **Chromosome split:** follow the paper, using chr1–19, chrX and chrY with primary alignments only. We have only 6
  samples, so it is optional for memory, but it parallelizes S/C/Q across 21 jobs.

## Status
Phase 0 is complete. Nothing is aligned yet.
