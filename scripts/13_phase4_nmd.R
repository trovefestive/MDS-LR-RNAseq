# =============================================================================
# 13_phase4_nmd.R — Phase 4: ORF / NMD characterization (Q157R − WT; WT = reference, rule 8)
# Inputs : results/03_sqanti/final/final_counts.tsv.gz, final_isoform_annotation.tsv.gz (R2 primary)
#          results/04_nmd/gene_names.tsv
# Outputs: results/04_nmd/*.tsv + plots
#   1) NMD fraction known vs novel among coding isoforms (Fisher's exact), also by SQANTI3 category
#   2) per-sample % of expression (TMM CPM) from NMD isoforms / novel NMD isoforms; Welch t-test (n=3 vs 3,
#      descriptive emphasis; no Wilcoxon: min two-sided p = 0.10)
#   3) genes with higher NMD-isoform share in Q157R: per gene, reads binned NMD vs non-NMD (coding);
#      edgeR diffSpliceDGE (bin usage, ~genotype); gene-level edgeR QL to exclude genes significantly down.
#      Δprop = mean(Q157R) − mean(WT). Thresholds: FDR < 0.05, |Δprop| > 0.10, prefilter ≥10 reads in ≥3 samples.
# Rule-8 assertions: (a) isoforms with 0 reads in all WT and >0 in all Q157R get logFC > 0 under the same
#   design; (b) significant bins have sign(logFC) == sign(Δprop).
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
fin <- file.path(proj, "results/03_sqanti/final"); out <- file.path(proj, "results/04_nmd")
dir.create(out, showWarnings = FALSE)
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))            # WT first = reference
S <- samples$sample_id; wt <- S[samples$genotype == "WT"]; mu <- S[samples$genotype == "Q157R"]
gn <- read_tsv(file.path(out, "gene_names.tsv"), show_col_types = FALSE)

cnt <- read_tsv(file.path(fin, "final_counts.tsv.gz"), show_col_types = FALSE)
ann <- read_tsv(file.path(fin, "final_isoform_annotation.tsv.gz"), show_col_types = FALSE, guess_max = 2e5) |>
  mutate(nmd = predicted_NMD %in% c(TRUE, "TRUE"), is_coding = coding == "coding",
         gene = associated_gene,
         gene_ok = str_starts(gene, "ENSMUSG") & !str_detect(gene, "_"))
stopifnot(identical(cnt$transcript_ID, ann$transcript_ID))
M <- as.matrix(cnt[S]); rownames(M) <- cnt$transcript_ID
y_all <- calcNormFactors(DGEList(M))
cpm_m <- cpm(y_all, normalized.lib.sizes = TRUE)

# === 1) NMD fraction, known vs novel (coding isoforms) ===
message("== 1) NMD known vs novel")
cod <- ann |> filter(is_coding)
tab <- table(origin = cod$origin, NMD = cod$nmd); print(tab)
ft <- fisher.test(tab)
nmd_kn <- cod |> group_by(origin) |> summarise(n_coding = n(), n_NMD = sum(nmd), pct_NMD = 100 * mean(nmd)) |>
  mutate(odds_ratio_novel_vs_GENCODE = (tab["novel", "TRUE"] / tab["novel", "FALSE"]) /
                                       (tab["GENCODE", "TRUE"] / tab["GENCODE", "FALSE"]),
         fisher_OR = unname(ft$estimate), fisher_p = ft$p.value)   # both = odds(NMD|novel)/odds(NMD|GENCODE)
write_tsv(nmd_kn, file.path(out, "phase4_nmd_known_vs_novel.tsv")); print(nmd_kn)
by_cat <- cod |> filter(origin == "novel") |> group_by(structural_category) |>
  summarise(n_coding = n(), pct_NMD = 100 * mean(nmd)) |> arrange(desc(n_coding))
write_tsv(by_cat, file.path(out, "phase4_nmd_novel_by_category.tsv")); print(by_cat)

# === 2) per-sample share of expression from NMD isoforms ===
message("== 2) per-sample NMD expression share (TMM CPM)")
share <- tibble(sample_id = S,
                pct_cpm_NMD = 100 * colSums(cpm_m[ann$nmd, ]) / colSums(cpm_m),
                pct_cpm_novel_NMD = 100 * colSums(cpm_m[ann$nmd & ann$origin == "novel", ]) / colSums(cpm_m),
                pct_NMD_cpm_that_is_novel = 100 * colSums(cpm_m[ann$nmd & ann$origin == "novel", ]) /
                  colSums(cpm_m[ann$nmd, ])) |>
  left_join(select(samples, sample_id, genotype), by = "sample_id")
write_tsv(share, file.path(out, "phase4_nmd_share_per_sample.tsv")); print(share)
tt <- map_dfr(c("pct_cpm_NMD", "pct_cpm_novel_NMD", "pct_NMD_cpm_that_is_novel"), \(v) {
  x <- share[[v]]; g <- share$genotype; t <- t.test(x[g == "Q157R"], x[g == "WT"])   # Welch; Q157R − WT
  tibble(metric = v, mean_WT = mean(x[g == "WT"]), mean_Q157R = mean(x[g == "Q157R"]),
         diff_Q157R_minus_WT = mean(x[g == "Q157R"]) - mean(x[g == "WT"]), welch_p = t$p.value)
})
write_tsv(tt, file.path(out, "phase4_nmd_share_ttest.tsv")); print(tt)

# === 3a) gene-level DE (all isoforms summed per gene) — used to exclude genes significantly DOWN ===
message("== 3a) gene-level edgeR (Q157R − WT)")
design <- model.matrix(~ genotype, data = samples); stopifnot(colnames(design)[2] == "genotypeQ157R")
gsum <- rowsum(M[ann$gene_ok, ], ann$gene[ann$gene_ok])
yg <- DGEList(gsum); yg <- yg[filterByExpr(yg, design), , keep.lib.sizes = FALSE]; yg <- calcNormFactors(yg)
yg <- estimateDisp(yg, design); fg <- glmQLFit(yg, design)
gde <- topTags(glmQLFTest(fg, coef = 2), n = Inf)$table |> rownames_to_column("gene") |>
  transmute(gene, gene_logFC_Q157R_vs_WT = logFC, gene_FDR = FDR)

# === rule-8 assertion (a): mutant-exclusive isoforms must get logFC > 0 under this design ===
excl <- rowSums(M[, wt] > 0) == 0 & rowSums(M[, mu] > 0) == length(mu) & rowSums(M) >= 5
message("mutant-exclusive isoforms (0 in all WT, >0 in all Q157R, >=5 reads): ", sum(excl))
ye <- DGEList(M[excl, , drop = FALSE], lib.size = y_all$samples$lib.size, norm.factors = y_all$samples$norm.factors)
le <- glmLRT(glmFit(ye, design, dispersion = yg$common.dispersion), coef = 2)
if (!all(le$table$logFC > 0)) stop("RULE 8 VIOLATION: mutant-exclusive isoform with logFC <= 0")
message("rule 8 (a) OK; U2af1-213 logFC = ",
        round(le$table[rownames(le$table) == "ENSMUST00000468653.1", "logFC"], 2))

# === 3b) NMD-bin usage per gene: diffSpliceDGE on {NMD, non-NMD coding} bins ===
message("== 3b) per-gene NMD share (diffSpliceDGE)")
bins <- ann |> filter(gene_ok, is_coding) |> mutate(bin = if_else(nmd, "NMD", "nonNMD"))
genes2 <- bins |> distinct(gene, bin) |> count(gene) |> filter(n == 2) |> pull(gene)
bins <- bins |> filter(gene %in% genes2)
bc <- rowsum(M[bins$transcript_ID, ], paste(bins$gene, bins$bin, sep = "|"))
bg <- tibble(id = rownames(bc)) |> separate(id, c("gene", "bin"), sep = "\\|", remove = FALSE)
gtot <- rowsum(bc, bg$gene)
keep_g <- rownames(gtot)[rowSums(gtot >= 10) >= 3]                         # prefilter: >=10 reads in >=3 samples
kb <- bg$gene %in% keep_g
yb <- DGEList(bc[kb, ], genes = bg[kb, ]); yb <- calcNormFactors(yb)
yb <- estimateDisp(yb, design); fb <- glmQLFit(yb, design)
ds <- diffSpliceDGE(fb, coef = 2, geneid = "gene", exonid = "bin", verbose = FALSE)
nb <- topSpliceDGE(ds, test = "exon", number = Inf) |> as_tibble() |> filter(bin == "NMD") |>
  transmute(gene, NMDbin_logFC = logFC, NMDbin_P = P.Value, NMDbin_FDR = FDR)
prop <- (bc[paste0(keep_g, "|NMD"), ] / gtot[keep_g, ]) |> as.data.frame() |> rownames_to_column("gene") |>
  mutate(gene = str_remove(gene, "\\|NMD$"))
res <- prop |> rowwise() |>
  mutate(prop_WT = mean(c_across(all_of(wt))), prop_Q157R = mean(c_across(all_of(mu))),
         dprop_Q157R_minus_WT = prop_Q157R - prop_WT) |> ungroup() |>
  left_join(nb, by = "gene") |> left_join(gde, by = "gene") |>
  left_join(gn, by = c("gene" = "gene_id")) |>
  mutate(gene_sig_down = !is.na(gene_FDR) & gene_FDR < 0.05 & gene_logFC_Q157R_vs_WT < 0,
         call = case_when(NMDbin_FDR < 0.05 & dprop_Q157R_minus_WT >  0.10 & !gene_sig_down ~ "NMD_share_UP_in_Q157R",
                          NMDbin_FDR < 0.05 & dprop_Q157R_minus_WT >  0.10 &  gene_sig_down ~ "NMD_share_up_but_gene_down",
                          NMDbin_FDR < 0.05 & dprop_Q157R_minus_WT < -0.10 ~ "NMD_share_DOWN_in_Q157R",
                          TRUE ~ "ns")) |>
  relocate(gene, gene_name, call, dprop_Q157R_minus_WT, NMDbin_FDR) |> arrange(NMDbin_P)

# === rule-8 assertion (b): significant NMD bins — logFC sign must match Δprop sign ===
chk <- res |> filter(NMDbin_FDR < 0.05, abs(dprop_Q157R_minus_WT) > 0.10)
if (any(sign(chk$NMDbin_logFC) != sign(chk$dprop_Q157R_minus_WT))) stop("RULE 8 VIOLATION: logFC/Δprop sign mismatch")
message("rule 8 (b) OK on ", nrow(chk), " significant genes")
write_tsv(res, file.path(out, "phase4_nmd_share_by_gene.tsv"))
message("genes tested: ", nrow(res)); print(count(res, call))
print(res |> filter(call != "ns") |> select(gene_name, call, prop_WT, prop_Q157R, dprop_Q157R_minus_WT, NMDbin_FDR,
                                            gene_logFC_Q157R_vs_WT) |> head(25), width = Inf)

# === 4) plots ===
p1 <- ggplot(nmd_kn, aes(origin, pct_NMD, fill = origin)) + geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(%s / %s)", pct_NMD, format(n_NMD, big.mark = ","),
                                format(n_coding, big.mark = ","))), vjust = -0.2, size = 3) +
  scale_fill_manual(values = c(GENCODE = "grey55", novel = "grey25")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "% of coding isoforms predicted NMD",
       title = sprintf("NMD: novel vs GENCODE (Fisher OR = %.2f, p = %.1e)", nmd_kn$fisher_OR[1], nmd_kn$fisher_p[1])) +
  theme_lr()
save_plot(p1, file.path(out, "phase4_nmd_known_vs_novel"), 5.5, 4)

sl <- share |> pivot_longer(c(pct_cpm_NMD, pct_cpm_novel_NMD), names_to = "metric", values_to = "pct") |>
  mutate(metric = recode(metric, pct_cpm_NMD = "all NMD isoforms", pct_cpm_novel_NMD = "novel NMD isoforms"))
p2 <- ggplot(sl, aes(genotype, pct, colour = genotype)) +
  geom_point(size = 3, position = position_jitter(width = 0.08, seed = 1)) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2.8, show.legend = FALSE, seed = 1) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "black", linewidth = 0.3) +
  scale_colour_manual(values = GENO_COLORS) + facet_wrap(~metric, scales = "free_y") +
  labs(x = NULL, y = "% of TMM CPM", title = "Expression share from NMD isoforms per sample") + theme_lr()
save_plot(p2, file.path(out, "phase4_nmd_share_per_sample"), 7, 4)

p3 <- ggplot(res |> filter(!is.na(NMDbin_FDR)), aes(dprop_Q157R_minus_WT, -log10(NMDbin_FDR), colour = call)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
  ggrepel::geom_text_repel(data = res |> filter(call != "ns") |> slice_min(NMDbin_FDR, n = 15),
                           aes(label = gene_name), size = 2.6, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = c(NMD_share_UP_in_Q157R = unname(GENO_COLORS["Q157R"]),
                                 NMD_share_DOWN_in_Q157R = unname(GENO_COLORS["WT"]),
                                 NMD_share_up_but_gene_down = "orange", ns = "grey75")) +
  labs(x = "Δ NMD-isoform share (Q157R − WT)", y = "-log10 FDR (diffSpliceDGE, NMD bin)", colour = NULL,
       title = "Per-gene NMD-isoform share, Q157R vs WT") + theme_lr()
save_plot(p3, file.path(out, "phase4_nmd_share_volcano"), 7, 5)

top <- res |> filter(call != "ns") |> slice_min(NMDbin_FDR, n = 12) |> pull(gene)
if (length(top)) {
  pt <- res |> filter(gene %in% top) |> select(gene, gene_name, all_of(S)) |>
    pivot_longer(all_of(S), names_to = "sample_id", values_to = "prop") |>
    left_join(select(samples, sample_id, genotype), by = "sample_id")
  p4 <- ggplot(pt, aes(genotype, prop, colour = genotype)) + geom_point(size = 2.5) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "black", linewidth = 0.3) +
    facet_wrap(~gene_name, scales = "free_y") + scale_colour_manual(values = GENO_COLORS) +
    labs(x = NULL, y = "NMD-isoform share of gene reads", title = "Top genes: NMD-isoform share per sample") +
    theme_lr(9)
  save_plot(p4, file.path(out, "phase4_nmd_top_genes"), 8, 6)
}
message("phase 4 done")
