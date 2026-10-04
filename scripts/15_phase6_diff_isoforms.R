# =============================================================================
# 15_phase6_diff_isoforms.R — Phase 6: differential isoform analysis cross-checks (Q157R − WT; WT = reference)
# Inputs : results/03_sqanti/final/{final_counts,final_isoform_annotation}.tsv.gz (R2 primary)
#          results/06_diff/rmats_long/ (rMATS-long, group 1 = Q157R, group 2 = WT)
#          results/04_nmd/gene_names.tsv; old LRP2 outputs (local only, comparison)
# Steps  : 1) edgeR gene-level (DGE) and transcript-level (DTE) QL tests
#          2) DRIMSeq DTU (isoform proportions within gene)
#          3) rMATS-long import + direction check
#          4) isoform switches without gene-level change; paper events (Ybx1); LRP2 comparison (negated)
#          5) per-isoform CPM dot plots
# Rule 8 : every effect = Q157R − WT. Assertions: mutant-exclusive isoforms -> DTE logFC > 0; rMATS-long delta for
#          U2af1-213 > 0; DRIMSeq Δprop sign = sign(mean prop Q157R − mean prop WT) by construction.
# Thresholds: FDR < 0.05; |Δprop| > 0.10 for usage; prefilter ≥10 reads in ≥3 samples.
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
fin <- file.path(proj, "results/03_sqanti/final"); out <- file.path(proj, "results/06_diff")
dir.create(out, showWarnings = FALSE)
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))
S <- samples$sample_id; wt <- S[samples$genotype == "WT"]; mu <- S[samples$genotype == "Q157R"]
design <- model.matrix(~ genotype, data = samples); stopifnot(colnames(design)[2] == "genotypeQ157R")
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)
strip_v <- function(x) sub("\\.\\d+$", "", x)

cnt <- read_tsv(file.path(fin, "final_counts.tsv.gz"), show_col_types = FALSE)
ann <- read_tsv(file.path(fin, "final_isoform_annotation.tsv.gz"), show_col_types = FALSE, guess_max = 2e5) |>
  mutate(gene = associated_gene, gene_ok = str_starts(gene, "ENSMUSG") & !str_detect(gene, "_")) |>
  left_join(select(gn, gene_id, gene_name), by = c("gene" = "gene_id"))
stopifnot(identical(cnt$transcript_ID, ann$transcript_ID), !anyDuplicated(cnt$transcript_ID))
M <- as.matrix(cnt[S]); rownames(M) <- cnt$transcript_ID

# === 1a) DGE: isoforms summed per gene ===
message("== 1) edgeR DGE / DTE")
qlf <- function(counts) {
  y <- DGEList(counts); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  y <- estimateDisp(y, design); f <- glmQLFit(y, design)
  list(y = y, tab = topTags(glmQLFTest(f, coef = 2), n = Inf)$table |> rownames_to_column("id"))
}
dge <- qlf(rowsum(M[ann$gene_ok, ], ann$gene[ann$gene_ok]))
dge_tab <- dge$tab |> rename(gene = id, logFC_Q157R_vs_WT = logFC) |> left_join(select(gn, gene_id, gene_name), by = c("gene" = "gene_id"))
# === 1b) DTE: isoforms ===
dte <- qlf(M)
dte_tab <- dte$tab |> rename(transcript_ID = id, logFC_Q157R_vs_WT = logFC) |>
  left_join(select(ann, transcript_ID, gene, gene_name, origin, structural_category, predicted_NMD), by = "transcript_ID")
# rule 8: mutant-exclusive isoforms that were tested must have logFC > 0
excl <- names(which(rowSums(M[, wt] > 0) == 0 & rowSums(M[, mu] > 0) == length(mu)))
chk <- dte_tab |> filter(transcript_ID %in% excl)
if (any(chk$logFC_Q157R_vs_WT <= 0)) stop("RULE 8 VIOLATION in DTE")
message("rule 8 DTE OK: ", nrow(chk), " mutant-exclusive isoforms tested, all logFC > 0; U2af1-213 logFC = ",
        round(dte_tab$logFC_Q157R_vs_WT[dte_tab$transcript_ID == "ENSMUST00000468653.1"], 2))
write_tsv(dge_tab, file.path(out, "phase6_DGE_edgeR_Q157R_vs_WT.tsv"))
write_tsv(dte_tab, file.path(out, "phase6_DTE_edgeR_Q157R_vs_WT.tsv"))
message("DGE: ", sum(dge_tab$FDR < 0.05), "/", nrow(dge_tab), " genes FDR<0.05 (up ", sum(dge_tab$FDR < 0.05 & dge_tab$logFC_Q157R_vs_WT > 0),
        "); DTE: ", sum(dte_tab$FDR < 0.05), "/", nrow(dte_tab), " isoforms FDR<0.05 (up ", sum(dte_tab$FDR < 0.05 & dte_tab$logFC_Q157R_vs_WT > 0), ")")

# === 2) DRIMSeq DTU ===
message("== 2) DRIMSeq DTU")
suppressPackageStartupMessages(library(DRIMSeq))
dd <- tibble(gene_id = ann$gene, feature_id = ann$transcript_ID) |> bind_cols(as_tibble(M)) |>
  filter(ann$gene_ok) |> as.data.frame()
ds <- dmDSdata(counts = dd, samples = data.frame(sample_id = S, group = samples$genotype))
ds <- dmFilter(ds, min_samps_gene_expr = 3, min_gene_expr = 10, min_samps_feature_expr = 2, min_feature_expr = 2,
               min_samps_feature_prop = 2, min_feature_prop = 0.05)
message("DRIMSeq genes after filter: ", length(ds))
set.seed(1)
ds <- dmPrecision(ds, design = design, BPPARAM = BiocParallel::MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))))
ds <- dmFit(ds, design = design, verbose = 0)
ds <- dmTest(ds, coef = "genotypeQ157R", verbose = 0)
dg <- results(ds) |> as_tibble() |> transmute(gene = gene_id, DTU_gene_p = pvalue, DTU_gene_FDR = adj_pvalue)
dfe <- results(ds, level = "feature") |> as_tibble() |> transmute(transcript_ID = feature_id, gene = gene_id, DTU_tx_p = pvalue, DTU_tx_FDR = adj_pvalue)
pr <- proportions(ds) |> as_tibble()
dtu_tab <- pr |> transmute(transcript_ID = feature_id, gene = gene_id,
                           prop_WT = rowMeans(across(all_of(wt))), prop_Q157R = rowMeans(across(all_of(mu))),
                           dprop_Q157R_minus_WT = prop_Q157R - prop_WT) |>
  left_join(dfe, by = c("transcript_ID", "gene")) |> left_join(dg, by = "gene") |>
  left_join(select(ann, transcript_ID, gene_name, origin, structural_category, predicted_NMD), by = "transcript_ID")
write_tsv(dtu_tab, file.path(out, "phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv"))
message("DTU: ", n_distinct(dtu_tab$gene[dtu_tab$DTU_gene_FDR < 0.05]), " genes FDR<0.05; ",
        sum(dtu_tab$DTU_tx_FDR < 0.05 & abs(dtu_tab$dprop_Q157R_minus_WT) > 0.10, na.rm = TRUE), " isoforms FDR<0.05 & |Δprop|>0.10")

# === 4a) isoform switches without gene-level change (YBX1-type) ===
message("== 4) switches, LRP2 comparison, paper events")
sw <- dtu_tab |> filter(DTU_tx_FDR < 0.05, abs(dprop_Q157R_minus_WT) > 0.10) |>
  left_join(select(dge_tab, gene, gene_logFC = logFC_Q157R_vs_WT, gene_FDR = FDR), by = "gene") |>
  mutate(gene_level_change = !is.na(gene_FDR) & gene_FDR < 0.05,
         switch_class = if_else(gene_level_change, "usage shift + gene DE", "usage shift, no gene-level change"))
write_tsv(sw, file.path(out, "phase6_isoform_switches.tsv")); print(count(sw, switch_class))

# === 4b) old LRP2 results (local, read-only); LRP2 A_vs_B uses A as reference -> negate to get Q157R − WT ===
lrp <- file.path(proj, "2026_04/260327_MDS_lrp2/S4_MULTISAMPLE_ANALYSIS/M2_DIFFERENTIAL_EXPRESSION")
cmp <- list()
if (dir.exists(lrp)) {
  od <- read_tsv(file.path(lrp, "differential_gene_expression/Q157R_vs_WT_DGE_edgeR_results.txt"), show_col_types = FALSE) |>
    transmute(g = strip_v(reference_gene_id), old_logFC_Q157R_vs_WT = -logFC, old_FDR = FDR)   # negated
  m1 <- dge_tab |> mutate(g = strip_v(gene)) |> inner_join(od, by = "g")
  cmp$DGE <- tibble(level = "gene (DGE)", n_shared = nrow(m1), r_logFC = cor(m1$logFC_Q157R_vs_WT, m1$old_logFC_Q157R_vs_WT),
                    sig_new = sum(m1$FDR < 0.05), sig_old = sum(m1$old_FDR < 0.05), sig_both = sum(m1$FDR < 0.05 & m1$old_FDR < 0.05),
                    sign_agree_both = mean(sign(m1$logFC_Q157R_vs_WT[m1$FDR < 0.05 & m1$old_FDR < 0.05]) ==
                                           sign(m1$old_logFC_Q157R_vs_WT[m1$FDR < 0.05 & m1$old_FDR < 0.05])))
  tm <- list.files(file.path(proj, "results/03_sqanti/compare"), "^vs_lrp2.*\\.tmap$", full.names = TRUE)
  map <- read_tsv(tm[1], show_col_types = FALSE) |> filter(class_code == "=") |> distinct(qry_id, ref_id)
  ot <- read_tsv(file.path(lrp, "differential_transcript_expression/Q157R_vs_WT_DTE_edgeR_results.txt"), show_col_types = FALSE) |>
    transmute(ref_id = isoform_id, old_logFC_Q157R_vs_WT = -logFC, old_FDR = FDR)                # negated
  m2 <- dte_tab |> inner_join(map, by = c("transcript_ID" = "qry_id")) |> inner_join(ot, by = "ref_id")
  cmp$DTE <- tibble(level = "isoform (DTE)", n_shared = nrow(m2), r_logFC = cor(m2$logFC_Q157R_vs_WT, m2$old_logFC_Q157R_vs_WT),
                    sig_new = sum(m2$FDR < 0.05), sig_old = sum(m2$old_FDR < 0.05), sig_both = sum(m2$FDR < 0.05 & m2$old_FDR < 0.05),
                    sign_agree_both = mean(sign(m2$logFC_Q157R_vs_WT[m2$FDR < 0.05 & m2$old_FDR < 0.05]) ==
                                           sign(m2$old_logFC_Q157R_vs_WT[m2$FDR < 0.05 & m2$old_FDR < 0.05])))
  ou <- read_tsv(file.path(lrp, "differential_transcript_usage/Q157R_vs_WT_DTU_transcript_DRIMSeq_summary.txt"), show_col_types = FALSE)
  # verify the old delta direction from its own group proportions before negating
  old_sign <- sign(cor(ou$delta_proportion, ou$group_prop_Q157R - ou$group_prop_WT, use = "complete.obs"))
  message("old LRP2 DTU delta_proportion vs (prop_Q157R - prop_WT): correlation sign = ", old_sign, " (-1 => negate)")
  ou <- ou |> transmute(ref_id = isoform_id, old_dprop_Q157R_minus_WT = old_sign * delta_proportion, old_FDR = adj_pvalue_transcript)
  m3 <- dtu_tab |> inner_join(map, by = c("transcript_ID" = "qry_id")) |> inner_join(ou, by = "ref_id")
  cmp$DTU <- tibble(level = "isoform usage (DTU)", n_shared = nrow(m3), r_logFC = cor(m3$dprop_Q157R_minus_WT, m3$old_dprop_Q157R_minus_WT, use = "complete.obs"),
                    sig_new = sum(m3$DTU_tx_FDR < 0.05, na.rm = TRUE), sig_old = sum(m3$old_FDR < 0.05, na.rm = TRUE),
                    sig_both = sum(m3$DTU_tx_FDR < 0.05 & m3$old_FDR < 0.05, na.rm = TRUE),
                    sign_agree_both = mean(with(filter(m3, DTU_tx_FDR < 0.05, old_FDR < 0.05),
                                                sign(dprop_Q157R_minus_WT) == sign(old_dprop_Q157R_minus_WT))))
  cmp_tab <- bind_rows(cmp); write_tsv(cmp_tab, file.path(out, "phase6_vs_LRP2.tsv")); print(cmp_tab, width = Inf)
  u <- m2 |> filter(str_detect(ref_id, "sedfdca6"))
  message("old U2af1::sedfdca6 <-> ", paste(u$transcript_ID, collapse = ","), "; new logFC ", round(u$logFC_Q157R_vs_WT, 2),
          " vs old (negated) ", round(u$old_logFC_Q157R_vs_WT, 2))
}

# === 4c) paper event: YBX1 exon-skipping novel isoform (SRSF2/U2AF1-mutant enriched in human AML/MDS) ===
ybx <- ann |> filter(gene_name == "Ybx1") |> select(transcript_ID, origin, structural_category, subcategory, predicted_NMD) |>
  left_join(select(dtu_tab, transcript_ID, prop_WT, prop_Q157R, dprop_Q157R_minus_WT, DTU_tx_FDR), by = "transcript_ID") |>
  left_join(select(dte_tab, transcript_ID, logFC_Q157R_vs_WT, FDR), by = "transcript_ID") |>
  bind_cols(as_tibble(M[match(ann$transcript_ID[ann$gene_name %in% "Ybx1"], rownames(M)), , drop = FALSE]))
write_tsv(ybx, file.path(out, "phase6_Ybx1_isoforms.tsv")); print(ybx, width = Inf)

# === 5) per-isoform CPM dot plots (every replicate shown) for key genes ===
message("== 5) dot plots")
y_all <- calcNormFactors(DGEList(M)); cpm_m <- cpm(y_all, normalized.lib.sizes = TRUE)
top_sw <- sw |> filter(!is.na(gene_name)) |> arrange(DTU_tx_FDR) |> distinct(gene_name) |> slice_head(n = 6) |> pull(gene_name)
key_genes <- unique(c("U2af1", "Ybx1", top_sw))
for (g in key_genes) {
  ids <- ann |> filter(gene_name == g) |> pull(transcript_ID)
  if (!length(ids)) next
  keep <- ids[rowMeans(cpm_m[ids, , drop = FALSE]) >= 1]                       # isoforms with mean CPM >= 1
  if (!length(keep)) next
  keep <- head(keep[order(-rowMeans(cpm_m[keep, , drop = FALSE]))], 8)
  pd <- as_tibble(cpm_m[keep, , drop = FALSE], rownames = "transcript_ID") |>
    pivot_longer(all_of(S), names_to = "sample_id", values_to = "CPM") |>
    left_join(select(samples, sample_id, genotype), by = "sample_id") |>
    left_join(select(ann, transcript_ID, origin), by = "transcript_ID") |>
    left_join(select(dtu_tab, transcript_ID, dprop_Q157R_minus_WT, DTU_tx_FDR), by = "transcript_ID") |>
    mutate(label = sprintf("%s (%s)\nΔprop %s, FDR %s", transcript_ID, origin,
                           ifelse(is.na(dprop_Q157R_minus_WT), "NA", sprintf("%+.2f", dprop_Q157R_minus_WT)),
                           ifelse(is.na(DTU_tx_FDR), "NA", formatC(DTU_tx_FDR, format = "g", digits = 2))))
  p <- ggplot(pd, aes(genotype, CPM, colour = genotype)) + geom_point(size = 2.3, position = position_jitter(width = 0.08, seed = 1)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = "black", linewidth = 0.25) +
    facet_wrap(~label, scales = "free_y", ncol = 4) + scale_colour_manual(values = GENO_COLORS) +
    labs(x = NULL, y = "TMM CPM", title = paste0(g, ": isoform expression per sample (Δprop = Q157R − WT)")) + theme_lr(8)
  save_plot(p, file.path(out, paste0("phase6_dotplot_", g)), 10, 2.6 * ceiling(length(keep) / 4) + 1)
}
message("phase 6 (edgeR/DRIMSeq part) done")
