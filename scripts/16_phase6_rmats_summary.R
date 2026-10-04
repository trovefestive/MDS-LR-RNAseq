# =============================================================================
# 16_phase6_rmats_summary.R — Phase 6: rMATS-long import, direction checks, comparison with DRIMSeq/edgeR
# rMATS-long run with group 1 = Q157R, group 2 = WT -> delta_isoform_proportion = Q157R − WT (rule 8).
# Inputs : results/06_diff/rmats_long/{differential_transcripts.tsv, results_by_gene/*/*isoform_differences*.tsv}
#          results/06_diff/phase6_{DTU_DRIMSeq,DGE_edgeR}_Q157R_vs_WT.tsv (15_phase6_diff_isoforms.R)
# =============================================================================
suppressPackageStartupMessages(library(tidyverse))
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
out <- file.path(proj, "results/06_diff"); rl <- file.path(out, "rmats_long")
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE)
wt <- samples$sample_id[samples$genotype == "WT"]; mu <- samples$sample_id[samples$genotype == "Q157R"]
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)
ann <- read_tsv(file.path(proj, "results/03_sqanti/final/final_isoform_annotation.tsv.gz"), show_col_types = FALSE,
                guess_max = 2e5) |> select(transcript_ID, origin, structural_category, subcategory, predicted_NMD)

# === 1) import + rule-8 checks ===
rt <- read_tsv(file.path(rl, "differential_transcripts.tsv"), show_col_types = FALSE, na = c("", "NA")) |>
  mutate(chk = rowMeans(across(all_of(paste0(mu, "_proportion")))) - rowMeans(across(all_of(paste0(wt, "_proportion")))))
bad <- rt |> filter(!is.na(delta_isoform_proportion), !is.na(chk), abs(chk) > 1e-3,
                    sign(delta_isoform_proportion) != sign(chk))
if (nrow(bad)) stop("RULE 8 VIOLATION: rMATS-long delta not equal to Q157R − WT for ", nrow(bad), " isoforms")
u <- rt |> filter(feature_id == "ENSMUST00000468653.1")
if (!nrow(u) || u$delta_isoform_proportion <= 0) stop("RULE 8 VIOLATION: U2af1-213 delta not > 0")
message("rule 8 OK: delta = Q157R − WT for all ", sum(!is.na(rt$delta_isoform_proportion)),
        " isoforms; U2af1-213 delta = ", u$delta_isoform_proportion, ", adj p = ", signif(u$adj_pvalue, 3))

# === 2) significant isoforms, annotated; usage switch with / without gene-level DE ===
dge <- read_tsv(file.path(out, "phase6_DGE_edgeR_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  transmute(gene_id = gene, gene_logFC = logFC_Q157R_vs_WT, gene_FDR = FDR)
sig <- rt |> filter(adj_pvalue < 0.05, abs(delta_isoform_proportion) >= 0.10) |>
  select(gene_id, feature_id, pvalue, adj_pvalue, Q157R_prop = group_1_average_proportion,
         WT_prop = group_2_average_proportion, delta_Q157R_minus_WT = delta_isoform_proportion) |>
  left_join(select(gn, gene_id, gene_name), by = "gene_id") |>
  left_join(ann, by = c("feature_id" = "transcript_ID")) |> left_join(dge, by = "gene_id") |>
  mutate(gene_level_change = !is.na(gene_FDR) & gene_FDR < 0.05) |> arrange(adj_pvalue)
write_tsv(sig, file.path(out, "phase6_rmatslong_significant.tsv"))
message("rMATS-long significant: ", nrow(sig), " isoforms in ", n_distinct(sig$gene_id), " genes; up in Q157R ",
        sum(sig$delta_Q157R_minus_WT > 0), "; novel ", sum(sig$origin == "novel", na.rm = TRUE),
        "; genes without gene-level DE ", n_distinct(sig$gene_id[!sig$gene_level_change]))

# === 3) agreement with DRIMSeq (same isoforms, same groups) ===
dtu <- read_tsv(file.path(out, "phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv"), show_col_types = FALSE)
m <- inner_join(select(rt, feature_id, delta_isoform_proportion, adj_pvalue),
                select(dtu, transcript_ID, dprop_Q157R_minus_WT, DTU_tx_FDR), by = c("feature_id" = "transcript_ID"))
s_r <- m$adj_pvalue < 0.05 & abs(m$delta_isoform_proportion) >= 0.10
s_d <- m$DTU_tx_FDR < 0.05 & abs(m$dprop_Q157R_minus_WT) > 0.10
agree <- tibble(n_tested_by_both = nrow(m),
                r_delta = cor(m$delta_isoform_proportion, m$dprop_Q157R_minus_WT, use = "complete.obs"),
                sig_rmatslong = sum(s_r, na.rm = TRUE), sig_drimseq = sum(s_d, na.rm = TRUE),
                sig_both = sum(s_r & s_d, na.rm = TRUE),
                sign_agree_both = mean(sign(m$delta_isoform_proportion[which(s_r & s_d)]) ==
                                       sign(m$dprop_Q157R_minus_WT[which(s_r & s_d)])))
write_tsv(agree, file.path(out, "phase6_rmatslong_vs_DRIMSeq.tsv")); print(agree, width = Inf)

# === 4) rMATS-long isoform-difference classes (top 100 significant genes) ===
fs <- list.files(file.path(rl, "results_by_gene"), "isoform_differences.*\\.tsv$", recursive = TRUE, full.names = TRUE)
ev <- map_dfr(fs, \(f) read_tsv(f, show_col_types = FALSE, col_types = cols(.default = "c")) |>
                mutate(gene_id = basename(dirname(f)))) |>
  left_join(select(gn, gene_id, gene_name), by = "gene_id")
write_tsv(ev, file.path(out, "phase6_rmatslong_isoform_differences.tsv"))
ev_tab <- ev |> count(event, name = "n_events") |> arrange(desc(n_events))
write_tsv(ev_tab, file.path(out, "phase6_rmatslong_event_counts.tsv")); print(ev_tab)
print(sig |> filter(!is.na(gene_name)) |> select(gene_name, feature_id, origin, structural_category, WT_prop, Q157R_prop,
                                                 delta_Q157R_minus_WT, adj_pvalue, gene_level_change) |> head(25), width = Inf)
message("phase 6 rMATS-long summary done")
