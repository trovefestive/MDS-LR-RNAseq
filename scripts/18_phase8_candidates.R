# =============================================================================
# 18_phase8_candidates.R — Phase 8: isoform-switch candidates for wet-lab follow-up, with evidence levels
# Every Δ is Q157R − WT (rule 8). Candidates = rMATS-long significant (adj p < 0.05, |Δ| ≥ 0.10) ∪ DRIMSeq significant.
# Evidence per isoform:
#   both_methods   rMATS-long and DRIMSeq both significant, same sign
#   lr_separated   long-read proportions do not overlap between genotypes (all 3 vs all 3)
#   lr_reads       ≥ 10 long reads per sample on average in the genotype where the isoform is higher
#   iq_agree       IsoQuant built the same intron chain (gffcompare "=") and its Δ has the same sign (NA: no model)
#   sr_replicated  short reads (separate LSK cohort, Salmon): same sign, |Δ| ≥ 0.05, gene ≥ 10 reads (NA: untestable)
#   sr_separated   short-read proportions do not overlap between genotypes
#   predicted_NMD, gene-level DE: reported, not scored
# Tier A: both_methods & lr_separated & sr_replicated & iq_agree not FALSE
# Tier B: both_methods & (sr_replicated | iq_agree)        Tier C: the rest
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(ggrepel) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
out <- file.path(proj, "results/08_summary"); dir.create(out, showWarnings = FALSE)
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE)
wt <- samples$sample_id[samples$genotype == "WT"]; mu <- samples$sample_id[samples$genotype == "Q157R"]
ss <- read_tsv(file.path(proj, "samples_shortread.tsv"), show_col_types = FALSE)
swt <- ss$sample_id[ss$genotype == "WT"]; smu <- ss$sample_id[ss$genotype == "Q157R"]
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)
ann <- read_tsv(file.path(proj, "results/03_sqanti/final/final_isoform_annotation.tsv.gz"), show_col_types = FALSE,
                guess_max = 2e5) |> select(transcript_ID, gene_ID, associated_gene, origin, structural_category, predicted_NMD)
sep <- function(hi, lo) apply(hi, 1, min) > apply(lo, 1, max)          # rows: all of hi above all of lo

# === 1) long reads: rMATS-long (all tested isoforms) + DRIMSeq ===
rl <- read_tsv(file.path(proj, "results/06_diff/rmats_long/differential_transcripts.tsv"), show_col_types = FALSE)
pm <- as.matrix(rl[paste0(mu, "_proportion")]); pw <- as.matrix(rl[paste0(wt, "_proportion")])
cm <- as.matrix(rl[paste0(mu, "_count")]);      cw <- as.matrix(rl[paste0(wt, "_count")])
rl <- rl |> transmute(transcript_ID = feature_id, gene_id, lr_delta = delta_isoform_proportion, lr_adjp = adj_pvalue,
                      lr_prop_WT = group_2_average_proportion, lr_prop_Q157R = group_1_average_proportion,
                      lr_reads_WT = rowMeans(cw), lr_reads_Q157R = rowMeans(cm),
                      lr_separated = if_else(lr_delta > 0, sep(pm, pw), sep(pw, pm)),
                      lr_reads = if_else(lr_delta > 0, lr_reads_Q157R, lr_reads_WT) >= 10)
dtu <- read_tsv(file.path(proj, "results/06_diff/phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  select(transcript_ID, drim_delta = dprop_Q157R_minus_WT, drim_FDR = DTU_tx_FDR, drim_gene_name = gene_name)
dge <- read_tsv(file.path(proj, "results/06_diff/phase6_DGE_edgeR_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  transmute(gene_id = gene, gene_logFC = logFC_Q157R_vs_WT, gene_FDR = FDR)
cand <- full_join(rl, dtu, by = "transcript_ID") |>
  mutate(sig_rl = !is.na(lr_adjp) & lr_adjp < 0.05 & abs(lr_delta) >= 0.10,
         sig_drim = !is.na(drim_FDR) & drim_FDR < 0.05 & abs(drim_delta) > 0.10) |>
  filter(sig_rl | sig_drim) |>
  mutate(delta = coalesce(lr_delta, drim_delta),
         both_methods = sig_rl & sig_drim & sign(lr_delta) == sign(drim_delta)) |>
  left_join(ann, by = "transcript_ID") |> mutate(gene_id = coalesce(gene_id, gene_ID)) |>
  left_join(select(gn, gene_id, gene_name), by = "gene_id") |>
  mutate(gene_name = coalesce(gene_name, drim_gene_name, associated_gene, gene_id)) |>
  left_join(dge, by = "gene_id")
stopifnot(all(sign(cand$lr_delta[cand$both_methods]) == sign(cand$drim_delta[cand$both_methods])))
message("candidates: ", nrow(cand), " isoforms (rMATS-long ", sum(cand$sig_rl), ", DRIMSeq ", sum(cand$sig_drim),
        ", both ", sum(cand$both_methods), ")")
# rule-8 positive control: U2af1-213 is up in Q157R in every source
u <- cand |> filter(transcript_ID == "ENSMUST00000468653.1")
if (nrow(u) && u$delta <= 0) stop("RULE 8 VIOLATION: U2af1-213 delta not > 0")

# gene grouping shared by IsoQuant-free measures: rMATS-long gene of every tested isoform, else the GTF gene
grp <- bind_rows(select(rl, transcript_ID, g = gene_id), transmute(ann, transcript_ID, g = gene_ID)) |>
  distinct(transcript_ID, .keep_all = TRUE)
gof <- set_names(grp$g, grp$transcript_ID)

# === 2) IsoQuant: same intron chain (gffcompare "=") and Δ of that model within its IsoQuant gene ===
iqd <- file.path(proj, "results/02_isoquant/LR")
tm <- read_tsv(file.path(proj, "results/03_sqanti/compare/vs_isoquant.R2_updated.gtf.tmap"), show_col_types = FALSE) |>
  filter(class_code == "=") |> select(transcript_ID = qry_id, iq_id = ref_id)
iqg <- read_tsv(file.path(iqd, "LR.transcript_models.gtf"), comment = "#", col_names = FALSE, show_col_types = FALSE) |>
  filter(X3 == "transcript") |>
  transmute(iq_gene = str_match(X9, 'gene_id "([^"]+)"')[, 2], iq_id = str_match(X9, 'transcript_id "([^"]+)"')[, 2])
iqc <- read_tsv(file.path(iqd, "LR.discovered_transcript_grouped_file_name_counts.tsv"), show_col_types = FALSE) |>
  rename(iq_id = 1) |> inner_join(iqg, by = "iq_id")
iqm <- as.matrix(iqc[samples$sample_id]); rownames(iqm) <- iqc$iq_id
iqp <- iqm / rowsum(iqm, iqc$iq_gene)[iqc$iq_gene, ]
iq <- tibble(iq_id = rownames(iqp), iq_delta = rowMeans(iqp[, mu, drop = FALSE]) - rowMeans(iqp[, wt, drop = FALSE]))
cand <- cand |> left_join(tm, by = "transcript_ID") |> left_join(iq, by = "iq_id") |>
  mutate(iq_agree = if_else(is.na(iq_delta), NA, sign(iq_delta) == sign(delta)))
message("IsoQuant: ", sum(!is.na(cand$iq_id)), " candidates with an identical IsoQuant model; Δ same sign ",
        sum(cand$iq_agree, na.rm = TRUE), " / ", sum(!is.na(cand$iq_agree)))

# === 3) short reads (separate LSK cohort): Salmon isoform share within the same gene grouping ===
sq <- map(set_names(ss$sample_id), \(s) read_tsv(file.path(proj, "results/07_shortread/salmon", s, "quant.sf"),
                                                 show_col_types = FALSE))
sc <- sapply(sq, \(x) x$NumReads); rownames(sc) <- sq[[1]]$Name
dup <- read_tsv(file.path(proj, "data/ref/expanded/salmon_index/duplicate_clusters.tsv"), show_col_types = FALSE)
seq_dup <- union(dup$RetainedRef, dup$DuplicateRef)
sc <- sc[rownames(sc) %in% names(gof), ]
sg <- rowsum(sc, gof[rownames(sc)])
cid <- intersect(cand$transcript_ID, rownames(sc))
sp <- sc[cid, , drop = FALSE] / sg[gof[cid], , drop = FALSE]
sr <- tibble(transcript_ID = cid, sr_gene_reads = rowMeans(sg[gof[cid], , drop = FALSE]),
             sr_prop_WT = rowMeans(sp[, swt, drop = FALSE]), sr_prop_Q157R = rowMeans(sp[, smu, drop = FALSE]),
             sr_up = sep(sp[, smu, drop = FALSE], sp[, swt, drop = FALSE]),
             sr_dn = sep(sp[, swt, drop = FALSE], sp[, smu, drop = FALSE])) |>
  mutate(sr_delta = sr_prop_Q157R - sr_prop_WT)
cand <- cand |> left_join(sr, by = "transcript_ID") |>
  mutate(sr_testable = !is.na(sr_delta) & sr_gene_reads >= 10 & !transcript_ID %in% seq_dup,
         sr_replicated = if_else(sr_testable, sign(sr_delta) == sign(delta) & abs(sr_delta) >= 0.05, NA),
         sr_separated = if_else(sr_testable, if_else(delta > 0, sr_up, sr_dn), NA)) |>
  select(-sr_up, -sr_dn)

# === 4) tiers + score ===
cand <- cand |>
  mutate(score = both_methods + lr_separated %in% TRUE + lr_reads %in% TRUE + iq_agree %in% TRUE +
                 sr_replicated %in% TRUE + sr_separated %in% TRUE,
         tier = case_when(both_methods & lr_separated %in% TRUE & sr_replicated %in% TRUE & !iq_agree %in% FALSE ~ "A",
                          both_methods & (sr_replicated %in% TRUE | iq_agree %in% TRUE) ~ "B",
                          TRUE ~ "C"),
         gene_level_DE = !is.na(gene_FDR) & gene_FDR < 0.05,
         direction = if_else(delta > 0, "up in Q157R", "down in Q157R")) |>
  arrange(tier, desc(score), desc(abs(delta)))
iso_cols <- c("tier", "score", "gene_name", "gene_id", "transcript_ID", "origin", "structural_category", "predicted_NMD",
              "direction", "delta", "lr_prop_WT", "lr_prop_Q157R", "lr_adjp", "drim_delta", "drim_FDR",
              "both_methods", "lr_separated", "lr_reads_WT", "lr_reads_Q157R", "lr_reads", "iq_id", "iq_delta",
              "iq_agree", "sr_prop_WT", "sr_prop_Q157R", "sr_delta", "sr_gene_reads", "sr_replicated", "sr_separated",
              "gene_logFC", "gene_FDR", "gene_level_DE")
iso_tab <- select(cand, all_of(iso_cols))
write_tsv(iso_tab, file.path(out, "phase8_candidates_isoforms.tsv"))
tiers <- iso_tab |> count(tier, name = "isoforms") |>
  left_join(iso_tab |> group_by(tier) |> summarise(genes = n_distinct(gene_id), novel = sum(origin == "novel", na.rm = TRUE),
                                                    sr_testable = sum(!is.na(sr_replicated)),
                                                    sr_replicated = sum(sr_replicated, na.rm = TRUE)), by = "tier")
write_tsv(tiers, file.path(out, "phase8_tier_counts.tsv")); print(tiers)

# gene-level list: best isoform per gene (tier, score, |Δ|) + its strongest opposite-direction partner
best <- iso_tab |> group_by(gene_id) |> slice(1) |> ungroup()
partner <- best |> select(gene_id, isoform = transcript_ID, delta) |>
  inner_join(select(rl, gene_id, partner_isoform = transcript_ID, partner_delta = lr_delta), by = "gene_id") |>
  filter(partner_isoform != isoform, sign(partner_delta) == -sign(delta)) |>
  group_by(gene_id) |> slice_max(abs(partner_delta), n = 1, with_ties = FALSE) |> ungroup() |>
  select(gene_id, partner_isoform, partner_delta)
genes <- best |> left_join(partner, by = "gene_id") |>
  transmute(tier, score, gene_name, gene_id, isoform = transcript_ID, origin, structural_category, predicted_NMD,
            delta, partner_isoform, partner_delta, both_methods, lr_separated, lr_reads,
            iq_agree, sr_delta, sr_replicated, sr_separated, gene_logFC, gene_FDR, gene_level_DE) |>
  arrange(tier, desc(score), desc(abs(delta)))
stopifnot(nrow(genes) == n_distinct(iso_tab$gene_id))
write_tsv(genes, file.path(out, "phase8_candidates_genes.tsv"))
message("genes: ", nrow(genes), " (tier A ", sum(genes$tier == "A"), ", B ", sum(genes$tier == "B"), ")")
print(genes |> filter(tier %in% c("A", "B")) |> select(tier, score, gene_name, isoform, origin, predicted_NMD, delta,
                                                      sr_delta, iq_agree, gene_level_DE) |> head(40), n = 40, width = Inf)

# === 5) figures ===
TIER_COLORS <- c(A = "#1B7837", B = "#E08214", C = "grey70")
d <- iso_tab |> filter(!is.na(sr_replicated))
lab <- d |> filter(tier == "A") |> group_by(gene_name) |> slice_max(abs(delta), n = 1, with_ties = FALSE) |> ungroup()
r_all <- cor(d$delta, d$sr_delta)
p1 <- ggplot(d, aes(delta, sr_delta)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = tier), size = 1.6, alpha = 0.8) +
  geom_text_repel(data = lab, aes(label = gene_name), size = 3, max.overlaps = 30, min.segment.length = 0) +
  scale_colour_manual(values = TIER_COLORS) +
  labs(x = "Long reads (LK): Δ isoform proportion, Q157R − WT", y = "Short reads (LSK): Δ isoform proportion, Q157R − WT",
       title = "Isoform-switch candidates: long vs short reads",
       subtitle = sprintf("%d candidates testable in short reads; r = %.2f; same sign %d%%", nrow(d), r_all,
                          round(100 * mean(sign(d$delta) == sign(d$sr_delta))))) + theme_lr()
save_plot(p1, file.path(out, "phase8_lr_vs_sr_delta"), 7, 6)

top <- genes |> filter(tier %in% c("A", "B")) |> head(30)
ev <- top |> transmute(gene = factor(paste0(gene_name, " (", tier, ")"), levels = rev(unique(paste0(gene_name, " (", tier, ")")))),
                       `rMATS-long + DRIMSeq` = both_methods, `LR replicates separate` = lr_separated,
                       `LR ≥ 10 reads` = lr_reads, `IsoQuant agrees` = iq_agree, `SR same sign` = sr_replicated,
                       `SR replicates separate` = sr_separated, `Predicted NMD` = predicted_NMD %in% TRUE,
                       `Gene-level DE` = gene_level_DE) |>
  distinct(gene, .keep_all = TRUE) |>
  pivot_longer(-gene, names_to = "evidence", values_to = "value") |>
  mutate(evidence = factor(evidence, levels = unique(evidence)),
         value = case_when(is.na(value) ~ "n/a", value ~ "yes", TRUE ~ "no"))
p2 <- ggplot(ev, aes(evidence, gene, fill = value)) + geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = c(yes = "#1B7837", no = "grey85", `n/a` = "#FBF3E4"), breaks = c("yes", "no", "n/a"), name = NULL) +
  labs(x = NULL, y = NULL, title = "Evidence for the top isoform-switch candidates (tier A/B)") +
  theme_lr() + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(),
                     axis.ticks = element_blank())
save_plot(p2, file.path(out, "phase8_evidence_top_genes"), 7.5, 0.22 * nrow(top) + 2.5)

ta <- genes |> filter(tier == "A") |> head(12)
if (nrow(ta)) {
  lr_long <- read_tsv(file.path(proj, "results/06_diff/rmats_long/differential_transcripts.tsv"), show_col_types = FALSE) |>
    filter(feature_id %in% ta$isoform) |> select(isoform = feature_id, all_of(paste0(samples$sample_id, "_proportion"))) |>
    pivot_longer(-isoform, names_to = "sample_id", values_to = "prop") |> mutate(sample_id = sub("_proportion$", "", sample_id)) |>
    left_join(select(samples, sample_id, genotype), by = "sample_id") |> mutate(data = "Long reads (LK)")
  sr_long <- as_tibble(sp[ta$isoform, , drop = FALSE], rownames = "isoform") |>
    pivot_longer(-isoform, names_to = "sample_id", values_to = "prop") |>
    left_join(select(ss, sample_id, genotype), by = "sample_id") |> mutate(data = "Short reads (LSK)")
  pd <- bind_rows(lr_long, sr_long) |> left_join(select(ta, isoform, gene_name), by = "isoform") |>
    mutate(panel = paste0(gene_name, "\n", isoform), genotype = factor(genotype, levels = GENO_LEVELS))
  p3 <- ggplot(pd, aes(genotype, prop, colour = genotype)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey40") +
    geom_point(size = 2, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    facet_grid(data ~ panel, scales = "free_y") + scale_colour_manual(values = GENO_COLORS, guide = "none") +
    labs(x = NULL, y = "Isoform proportion within gene", title = "Tier A switches: every replicate, both datasets") +
    theme_lr(9) + theme(strip.text.x = element_text(size = 6.5))
  save_plot(p3, file.path(out, "phase8_tierA_proportions"), 1.6 * nrow(ta) + 1.5, 4.2)
}
message("phase 8 candidates done")
