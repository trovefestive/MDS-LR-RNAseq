# =============================================================================
# 20_gene_panel.R — expression and splicing status of a gene panel in long reads (LK) and short reads (LSK)
#   usage: Rscript scripts/20_gene_panel.R Gata1 Med12 Pik3r5 Gdf11 Acvr2b Smad2 Smad3   (run via submit or lr-r env)
# Every effect is Q157R − WT (rule 8). Reads existing outputs only; nothing is recomputed except per-gene CPM.
#   DGE   : long reads edgeR (Phase 6); short reads edgeR on Salmon gene sums (Phase 7)
#   Usage : rMATS-long + DRIMSeq (long reads); Salmon isoform shares (short reads)
#   Events: rMATS-turbo JCEC (short reads); LeafCutter clusters (both, Phase 7b)
# Out    : results/09_gene_panel/panel_{summary,DGE,isoforms,rmats_turbo,leafcutter}.tsv
# =============================================================================
suppressPackageStartupMessages(library(tidyverse))
proj <- Sys.getenv("PROJ"); if (proj == "") proj <- getwd()
genes <- commandArgs(trailingOnly = TRUE)
if (!length(genes)) genes <- c("Gata1", "Med12", "Pik3r5", "Gdf11", "Acvr2b", "Smad2", "Smad3")
out <- file.path(proj, "results/09_gene_panel"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
f <- \(...) file.path(proj, ...)
gn <- read_tsv(f("results/04_nmd/gene_names.tsv"), show_col_types = FALSE) |> filter(gene_name %in% genes)
miss <- setdiff(genes, gn$gene_name); if (length(miss)) message("not in annotation: ", paste(miss, collapse = ", "))
ids <- set_names(gn$gene_id, gn$gene_name)
samples <- read_tsv(f("samples.tsv"), show_col_types = FALSE)
wt <- samples$sample_id[samples$genotype == "WT"]; mu <- samples$sample_id[samples$genotype == "Q157R"]
ss <- read_tsv(f("samples_shortread.tsv"), show_col_types = FALSE)
swt <- ss$sample_id[ss$genotype == "WT"]; smu <- ss$sample_id[ss$genotype == "Q157R"]

# === 1) gene-level expression ===
ann <- read_tsv(f("results/03_sqanti/final/final_isoform_annotation.tsv.gz"), show_col_types = FALSE, guess_max = 2e5) |>
  select(transcript_ID, gene = associated_gene, origin, structural_category, predicted_NMD)
cpm <- read_tsv(f("results/03_sqanti/final/final_cpm_tmm.tsv.gz"), show_col_types = FALSE)
lr_cpm <- cpm |> inner_join(select(ann, transcript_ID, gene), by = "transcript_ID") |> filter(gene %in% ids) |>
  group_by(gene) |> summarise(across(all_of(c(wt, mu)), sum), .groups = "drop") |>
  mutate(LR_CPM_WT = rowMeans(across(all_of(wt))), LR_CPM_Q157R = rowMeans(across(all_of(mu))),
         LR_CPM_WT_reps = pmap_chr(across(all_of(wt)), \(...) paste(sprintf("%.1f", c(...)), collapse = "/")),
         LR_CPM_Q157R_reps = pmap_chr(across(all_of(mu)), \(...) paste(sprintf("%.1f", c(...)), collapse = "/"))) |>
  select(gene, LR_CPM_WT, LR_CPM_Q157R, LR_CPM_WT_reps, LR_CPM_Q157R_reps)
lr_dge <- read_tsv(f("results/06_diff/phase6_DGE_edgeR_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  filter(gene %in% ids) |> transmute(gene, LR_logFC = logFC_Q157R_vs_WT, LR_FDR = FDR)
# short-read gene DE: same model as Phase 7 (edgeR QL on Salmon gene sums, ~ genotype, WT reference), all genes,
# so genes not tested in long reads (e.g. Acvr2b, Gdf11) are included too
suppressPackageStartupMessages(library(edgeR))
sq_all <- map(set_names(ss$sample_id), \(s) read_tsv(f("results/07_shortread/salmon", s, "quant.sf"), show_col_types = FALSE))
cnt <- sapply(sq_all, \(x) x$NumReads); rownames(cnt) <- sq_all[[1]]$Name
g_of <- set_names(ann$gene, ann$transcript_ID); keep <- rownames(cnt) %in% names(g_of)
gg <- g_of[rownames(cnt)[keep]]; okg <- str_starts(gg, "ENSMUSG") & !str_detect(gg, "_")
g_sr <- rowsum(round(cnt[keep, ][okg, ]), gg[okg])
grp <- factor(ss$genotype[match(colnames(g_sr), ss$sample_id)], levels = c("WT", "Q157R"))
design <- model.matrix(~ grp); stopifnot(colnames(design)[2] == "grpQ157R")
y <- DGEList(g_sr); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
y <- estimateDisp(y, design); fit <- glmQLFit(y, design)
sr_dge <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table |> rownames_to_column("gene") |> filter(gene %in% ids) |>
  transmute(gene, SR_logFC = logFC, SR_FDR = FDR)
message("short-read DE genes tested: ", nrow(y), "; panel genes passing filterByExpr: ", nrow(sr_dge))
# short-read TPM per gene from Salmon (also for genes edgeR filtered out)
sq <- map_dfr(ss$sample_id, \(s) read_tsv(f("results/07_shortread/salmon", s, "quant.sf"), show_col_types = FALSE) |>
                transmute(transcript_ID = Name, TPM, NumReads, sample = s))
sr_tx <- sq |> inner_join(select(ann, transcript_ID, gene), by = "transcript_ID") |> filter(gene %in% ids)
sr_tpm <- sr_tx |> group_by(gene, sample) |> summarise(TPM = sum(TPM), reads = sum(NumReads), .groups = "drop") |>
  group_by(gene) |> summarise(SR_TPM_WT = mean(TPM[sample %in% swt]), SR_TPM_Q157R = mean(TPM[sample %in% smu]),
                              SR_reads_mean = mean(reads), .groups = "drop")
dge <- tibble(gene_name = names(ids), gene = unname(ids)) |> left_join(lr_cpm, by = "gene") |> left_join(lr_dge, by = "gene") |>
  left_join(sr_tpm, by = "gene") |> left_join(sr_dge, by = "gene")
write_tsv(dge, file.path(out, "panel_DGE.tsv"))

# === 2) isoform usage: rMATS-long + DRIMSeq (long reads), Salmon shares (short reads) ===
rl <- read_tsv(f("results/06_diff/rmats_long/differential_transcripts.tsv"), show_col_types = FALSE) |>
  filter(gene_id %in% ids) |>
  transmute(gene = gene_id, transcript_ID = feature_id, LR_prop_WT = group_2_average_proportion,
            LR_prop_Q157R = group_1_average_proportion, LR_delta = delta_isoform_proportion, LR_adjp = adj_pvalue)
dtu <- read_tsv(f("results/06_diff/phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  filter(gene %in% ids) |> select(transcript_ID, DRIM_delta = dprop_Q157R_minus_WT, DRIM_FDR = DTU_tx_FDR, DRIM_gene_FDR = DTU_gene_FDR)
srp <- sr_tx |> group_by(gene, sample) |> mutate(share = NumReads / sum(NumReads)) |> ungroup() |>
  group_by(gene, transcript_ID) |> summarise(SR_share_WT = mean(share[sample %in% swt]), SR_share_Q157R = mean(share[sample %in% smu]),
                                             .groups = "drop") |> mutate(SR_delta = SR_share_Q157R - SR_share_WT)
iso <- full_join(rl, select(srp, -gene), by = "transcript_ID") |> mutate(gene = coalesce(gene, srp$gene[match(transcript_ID, srp$transcript_ID)])) |>
  left_join(dtu, by = "transcript_ID") |> left_join(select(ann, transcript_ID, origin, structural_category, predicted_NMD), by = "transcript_ID") |>
  left_join(tibble(gene = unname(ids), gene_name = names(ids)), by = "gene") |>
  filter(coalesce(LR_prop_WT, 0) + coalesce(LR_prop_Q157R, 0) + coalesce(SR_share_WT, 0) + coalesce(SR_share_Q157R, 0) >= 0.05) |>
  relocate(gene_name) |> arrange(gene_name, desc(coalesce(LR_prop_WT, SR_share_WT)))
write_tsv(iso, file.path(out, "panel_isoforms.tsv"))

# === 3) rMATS-turbo events (short reads) ===
ev_cols <- list(SE = c("exonStart_0base", "exonEnd"), A3SS = c("longExonStart_0base", "longExonEnd"),
                A5SS = c("longExonStart_0base", "longExonEnd"), RI = c("riExonStart_0base", "riExonEnd"),
                MXE = c("1stExonStart_0base", "2ndExonEnd"))
rt <- imap_dfr(ev_cols, \(cc, e) {
  x <- read_tsv(f("results/07_shortread/rmats/out", paste0(e, ".MATS.JCEC.txt")), show_col_types = FALSE, name_repair = "unique_quiet") |>
    mutate(GeneID = str_remove_all(GeneID, '"')) |> filter(GeneID %in% ids)
  tibble(event = e, gene = x$GeneID, chr = x$chr, strand = x$strand, region = paste0(x[[cc[1]]] + 1, "-", x[[cc[2]]]),
         IJC_Q157R = x$IJC_SAMPLE_1, SJC_Q157R = x$SJC_SAMPLE_1, IJC_WT = x$IJC_SAMPLE_2, SJC_WT = x$SJC_SAMPLE_2,
         PSI_Q157R = x$IncLevel1, PSI_WT = x$IncLevel2, dPSI_Q157R_minus_WT = x$IncLevelDifference, FDR = x$FDR)
})
rt <- rt |> left_join(tibble(gene = unname(ids), gene_name = names(ids)), by = "gene") |> relocate(gene_name) |> arrange(gene_name, FDR)
write_tsv(rt, file.path(out, "panel_rmats_turbo.tsv"))

# === 4) LeafCutter clusters (both datasets) ===
lc <- map_dfr(c("shortread", "longread"), \(set) {
  d <- f("results/07_leafcutter", set)
  cs <- read_tsv(file.path(d, paste0(set, "_ds_cluster_significance.txt")), show_col_types = FALSE)
  es <- read_tsv(file.path(d, paste0(set, "_ds_effect_sizes.txt")), show_col_types = FALSE) |>
    mutate(cluster = str_replace(intron, "^([^:]+):[^:]+:[^:]+:", "\\1:"))
  cs |> filter(map_lgl(str_split(coalesce(genes, ""), ","), \(g) any(g %in% names(ids)))) |>
    left_join(es |> group_by(cluster) |> summarise(n_introns = n(), max_abs_dPSI = max(abs(deltapsi)),
                                                    top_intron = intron[which.max(abs(deltapsi))],
                                                    top_dPSI_Q157R_minus_WT = deltapsi[which.max(abs(deltapsi))], .groups = "drop"),
              by = "cluster") |>
    transmute(set, gene_name = genes, cluster, status, p.adjust, n_introns, max_abs_dPSI, top_intron, top_dPSI_Q157R_minus_WT)
})
write_tsv(lc, file.path(out, "panel_leafcutter.tsv"))

# === 5) one-line summary per gene ===
sig_rt <- rt |> filter(FDR < 0.05, abs(dPSI_Q157R_minus_WT) >= 0.1)
summ <- dge |> transmute(gene_name, LR_CPM_WT, LR_CPM_Q157R, LR_logFC, LR_FDR, SR_TPM_WT, SR_TPM_Q157R, SR_logFC, SR_FDR) |>
  left_join(iso |> group_by(gene_name) |> summarise(LR_isoforms_tested = sum(!is.na(LR_adjp)),
              LR_usage_sig = sum(LR_adjp < 0.05 & abs(LR_delta) >= 0.1, na.rm = TRUE),
              DRIM_gene_FDR = suppressWarnings(min(DRIM_gene_FDR, na.rm = TRUE)), .groups = "drop"), by = "gene_name") |>
  left_join(rt |> group_by(gene_name) |> summarise(SR_events_tested = n(), .groups = "drop"), by = "gene_name") |>
  left_join(sig_rt |> group_by(gene_name) |> summarise(SR_events_sig = n(),
              SR_sig_events = paste0(event, " ", region, " dPSI ", sprintf("%+.2f", dPSI_Q157R_minus_WT), collapse = "; "), .groups = "drop"),
            by = "gene_name") |>
  left_join(lc |> separate_rows(gene_name, sep = ",") |> group_by(gene_name) |>
              summarise(LC_SR_min_padj = suppressWarnings(min(p.adjust[set == "shortread"], na.rm = TRUE)),
                        LC_LR_min_padj = suppressWarnings(min(p.adjust[set == "longread"], na.rm = TRUE)), .groups = "drop"), by = "gene_name") |>
  mutate(across(where(is.numeric), \(v) if_else(is.infinite(v), NA_real_, v)), across(c(SR_events_sig), \(v) coalesce(v, 0L)))
write_tsv(summ, file.path(out, "panel_summary.tsv"))
options(width = 250); print(as.data.frame(summ))
message("sign convention: every logFC / delta / dPSI is Q157R - WT")
