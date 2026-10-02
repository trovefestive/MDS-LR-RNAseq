# =============================================================================
# 12_phase3_summary.R — Phase 3: primary matrix, novel-isoform rule, CPM, summaries
# Inputs : results/03_sqanti/quant/{R2,fsm}_abundance.esp, R2_updated.gtf
#          results/03_sqanti/primary_filter_classification.txt (link to the P3_FILTER table, set by 03b)
#          results/03_sqanti/compare/*.tmap (gffcompare vs IsoQuant / old LRP2)
# Outputs: results/03_sqanti/final/*  and  results/03_sqanti/phase3_*.tsv / plots
# Rules (ANALYSIS_PLAN Phase 3):
#   - primary = fsm matrix only if every sample loses <20% of reads to the FSM filter; else R2
#   - final set = all GENCODE transcripts + novel isoforms with >=1 read in >=2 samples and >=5 reads
#   - flag novel isoforms with >=1 read in >=2 of 3 replicates of one genotype
#   - CPM = edgeR TMM (as in the paper's espresso_to_cpm.R)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
d  <- file.path(proj, "results/03_sqanti"); q <- file.path(d, "quant")
fo <- file.path(d, "final"); dir.create(fo, showWarnings = FALSE)
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))

read_esp <- function(f) {
  read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) |>
    mutate(across(all_of(samples$sample_id), as.numeric))
}

# === 1) read loss from the FSM filter -> primary matrix ===
message("reading R2 and fsm matrices")
r2  <- read_esp(file.path(q, "R2_abundance.esp"))
fsm <- read_esp(file.path(q, "fsm_abundance.esp"))
loss <- tibble(sample_id = samples$sample_id,
               reads_R2  = colSums(r2[samples$sample_id]),
               reads_fsm = colSums(fsm[samples$sample_id])) |>
  left_join(select(samples, sample_id, genotype), by = "sample_id") |>
  mutate(pct_lost_fsm = 100 * (1 - reads_fsm / reads_R2))
print(loss)
primary_name <- if (all(loss$pct_lost_fsm < 20)) "fsm" else "R2"
message("PRIMARY MATRIX: ", primary_name,
        " (max FSM read loss ", sprintf("%.1f", max(loss$pct_lost_fsm)), "%; rule: fsm only if all < 20%)")
write_tsv(loss |> mutate(primary = primary_name), file.path(d, "phase3_fsm_read_loss.tsv"))
m <- if (primary_name == "fsm") fsm else r2

# === 2) annotation from SQANTI3 (ML filter table) ===
cls <- read_tsv(file.path(d, "primary_filter_classification.txt"), show_col_types = FALSE,   # P3_FILTER (03b)
                guess_max = 1e5) |>
  select(isoform, chrom, strand, length, exons, structural_category, subcategory, associated_gene,
         associated_transcript, coding, predicted_NMD, any_of(c("within_CAGE_peak", "polyA_motif_found",
         "perc_A_downstream_TTS", "RTS_stage", "all_canonical")), filter_result)
ann <- m |> select(transcript_ID, transcript_name, gene_ID) |>
  left_join(cls, by = c("transcript_ID" = "isoform")) |>
  mutate(origin = if_else(str_starts(transcript_ID, "ENSMUST"), "GENCODE", "novel"))

# === 3) novel-isoform rule ===
cnt <- as.matrix(m[samples$sample_id]); rownames(cnt) <- m$transcript_ID
n_samp <- rowSums(cnt >= 1); tot <- rowSums(cnt)
wt <- samples$sample_id[samples$genotype == "WT"]; mu <- samples$sample_id[samples$genotype == "Q157R"]
ann <- ann |> mutate(n_samples_ge1 = n_samp, total_reads = tot,
                     n_WT_ge1 = rowSums(cnt[, wt] >= 1), n_Q157R_ge1 = rowSums(cnt[, mu] >= 1),
                     pass_novel_rule = origin == "GENCODE" | (n_samples_ge1 >= 2 & total_reads >= 5),
                     within_genotype_support = n_WT_ge1 >= 2 | n_Q157R_ge1 >= 2)
keep <- ann$pass_novel_rule
message("final set: ", sum(keep), " isoforms (GENCODE ", sum(keep & ann$origin == "GENCODE"),
        ", novel ", sum(keep & ann$origin == "novel"), "); novel with within-genotype support: ",
        sum(keep & ann$origin == "novel" & ann$within_genotype_support))

# === 4) TMM CPM on the final set ===
y <- DGEList(counts = cnt[keep, ]); y <- calcNormFactors(y)
cpm_m <- cpm(y, normalized.lib.sizes = TRUE)
fin <- ann[keep, ]
write_tsv(bind_cols(m[keep, 1:3], as_tibble(round(cnt[keep, ], 2))), file.path(fo, "final_counts.tsv.gz"))
write_tsv(bind_cols(m[keep, 1:3], as_tibble(round(cpm_m, 4))), file.path(fo, "final_cpm_tmm.tsv.gz"))
write_tsv(fin, file.path(fo, "final_isoform_annotation.tsv.gz"))
writeLines(fin$transcript_ID, file.path(fo, "final_isoform_ids.txt"))

# === 5) summaries ===
fin <- fin |> mutate(mean_cpm = rowMeans(cpm_m))
cat_tab <- fin |> count(origin, structural_category, name = "n") |> arrange(origin, desc(n))
write_tsv(cat_tab, file.path(d, "phase3_structural_categories.tsv")); print(cat_tab, n = 40)
nov_frac <- tibble(sample_id = colnames(cpm_m),
                   pct_cpm_novel = 100 * colSums(cpm_m[fin$origin == "novel", , drop = FALSE]) / colSums(cpm_m)) |>
  left_join(select(samples, sample_id, genotype), by = "sample_id")
write_tsv(nov_frac, file.path(d, "phase3_novel_expression_fraction.tsv")); print(nov_frac)
gene_share <- fin |> mutate(gene = coalesce(associated_gene, gene_ID)) |> group_by(gene) |>
  mutate(share = mean_cpm / sum(mean_cpm)) |> ungroup()
n_major_novel <- sum(gene_share$origin == "novel" & gene_share$share > 0.10, na.rm = TRUE)
summ <- tibble(metric = c("primary_matrix", "isoforms_final", "gencode_final", "novel_final",
                          "novel_within_genotype_support", "novel_gt10pct_of_gene", "novel_coding",
                          "novel_predicted_NMD"),
               value = c(primary_name, sum(keep), sum(fin$origin == "GENCODE"), sum(fin$origin == "novel"),
                         sum(fin$origin == "novel" & fin$within_genotype_support), n_major_novel,
                         sum(fin$origin == "novel" & fin$coding == "coding", na.rm = TRUE),
                         sum(fin$origin == "novel" & fin$predicted_NMD %in% c(TRUE, "TRUE"), na.rm = TRUE)))
write_tsv(summ, file.path(d, "phase3_summary.tsv")); print(summ)

# === 6) overlap with IsoQuant and old LRP2 (gffcompare intron-chain match "=") ===
cmp_dir <- file.path(d, "compare")
for (nm in c("vs_isoquant", "vs_lrp2")) {
  tm <- list.files(cmp_dir, pattern = paste0("^", nm, ".*\\.tmap$"), full.names = TRUE)
  if (!length(tm)) next
  t <- read_tsv(tm[1], show_col_types = FALSE) |> select(qry_id, class_code)
  o <- fin |> left_join(t, by = c("transcript_ID" = "qry_id")) |>
    mutate(match = class_code == "=") |> group_by(origin) |>
    summarise(n = n(), n_intron_chain_match = sum(match, na.rm = TRUE), pct = 100 * mean(match, na.rm = TRUE))
  write_tsv(o, file.path(d, paste0("phase3_overlap_", nm, ".tsv"))); message(nm); print(o)
}

# === 6b) IsoQuant-centric recovery: share of IsoQuant multi-exon models whose intron chain is in our set ===
iq_gtf <- file.path(proj, "results/02_isoquant/LR/LR.transcript_models.gtf")
rmap <- list.files(cmp_dir, pattern = "^vs_isoquant.*\\.refmap$", full.names = TRUE)
if (file.exists(iq_gtf) && length(rmap)) {
  iq <- read_tsv(iq_gtf, comment = "#", col_names = FALSE, show_col_types = FALSE,
                 col_types = cols(.default = "c")) |>
    filter(X3 == "exon") |> mutate(tid = str_match(X9, 'transcript_id "([^"]+)"')[, 2]) |>
    count(tid, name = "n_exons") |> filter(n_exons > 1) |>
    mutate(iq_class = if_else(str_starts(tid, "ENSMUST"), "IsoQuant known", "IsoQuant novel"))
  rm_eq <- read_tsv(rmap[1], show_col_types = FALSE) |> filter(class_code == "=") |> pull(ref_id) |> unique()
  iq_rec <- iq |> group_by(iq_class) |>
    summarise(n_multi_exon = n(), recovered_in_final_set = sum(tid %in% rm_eq),
              pct = 100 * recovered_in_final_set / n_multi_exon)
  write_tsv(iq_rec, file.path(d, "phase3_isoquant_recovery.tsv")); message("IsoQuant-centric recovery"); print(iq_rec)
}

# === 7) plots ===
p1 <- ggplot(cat_tab |> filter(origin == "novel"),
             aes(reorder(structural_category, n), n)) + geom_col(fill = "grey40") + coord_flip() +
  labs(x = NULL, y = "novel isoforms (final set)", title = "Novel isoforms by SQANTI3 category") + theme_lr()
save_plot(p1, file.path(d, "phase3_novel_categories"), 6, 4)
p2 <- ggplot(fin, aes(length, colour = origin)) + geom_density() + scale_x_log10() +
  labs(x = "transcript length (bp)", y = "density", title = "Length: GENCODE vs novel") + theme_lr()
save_plot(p2, file.path(d, "phase3_length_known_vs_novel"), 6, 4)
p3 <- ggplot(fin |> filter(mean_cpm > 0), aes(mean_cpm, colour = origin)) + geom_density() + scale_x_log10() +
  labs(x = "mean TMM CPM", y = "density", title = "Expression: GENCODE vs novel") + theme_lr()
save_plot(p3, file.path(d, "phase3_expression_known_vs_novel"), 6, 4)
p4 <- ggplot(loss, aes(sample_id, pct_lost_fsm, fill = genotype)) + geom_col(width = 0.7) +
  geom_hline(yintercept = 20, linetype = 2) + scale_fill_manual(values = GENO_COLORS) +
  facet_grid(~genotype, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "% reads removed by FSM filter", title = "Paper FSM-read filter: read loss (dashed = 20% rule)") +
  theme_lr()
save_plot(p4, file.path(d, "phase3_fsm_read_loss"), 7, 4)
message("phase 3 summary done")
