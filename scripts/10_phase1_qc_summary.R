# =============================================================================
# 10_phase1_qc_summary.R — Phase 1: combine per-sample QC, plot, genotype check
# Inputs : results/01_align/{qc,genotype}/*
# Outputs: results/01_align/phase1_qc_table.tsv, phase1_genotype_check.tsv, plots
# Stops with an error if any sample's Q157R call disagrees with samples.tsv.
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(ggrepel) })
proj <- Sys.getenv("PROJ")   # exported by config.sh (01b_qc_summary.slurm)
if (proj == "") stop("PROJ not set: run via scripts/submit.sh scripts/01b_qc_summary.slurm")
source(file.path(proj, "scripts/00_plot_theme.R"))
d_qc <- file.path(proj, "results/01_align/qc")
d_gt <- file.path(proj, "results/01_align/genotype")
d_out <- file.path(proj, "results/01_align")

samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))

# === seqkit read stats ===
message("reading seqkit stats")
sk <- map_dfr(samples$sample_id, \(s)
  read_tsv(file.path(d_qc, paste0(s, ".seqkit_stats.tsv")), show_col_types = FALSE) |>
    transmute(sample_id = s, flnc_reads = num_seqs, flnc_bases = sum_len,
              mean_len = avg_len, median_len = Q2, N50))

# === samtools stats (SN section) ===
message("reading samtools stats")
sn <- map_dfr(samples$sample_id, \(s) {
  x <- read_tsv(file.path(d_qc, paste0(s, ".samtools_stats_SN.tsv")),
                col_names = c("key", "value", "comment"), show_col_types = FALSE)
  v <- \(k) as.numeric(x$value[x$key == k])
  tibble(sample_id = s,
         seqs = v("sequences:"), reads_mapped = v("reads mapped:"),
         reads_unmapped = v("reads unmapped:"),
         error_rate = v("error rate:"))
})

qc <- samples |> select(sample_id, genotype, barcode) |>
  left_join(sk, by = "sample_id") |> left_join(sn, by = "sample_id") |>
  mutate(pct_mapped = 100 * reads_mapped / seqs)
write_tsv(qc, file.path(d_out, "phase1_qc_table.tsv"))
print(qc, width = Inf)

# === read-length distributions ===
lh <- map_dfr(samples$sample_id, \(s)
  read_tsv(file.path(d_qc, paste0(s, ".length_hist.tsv")), show_col_types = FALSE)) |>
  rename(sample_id = sample) |>
  left_join(select(samples, sample_id, genotype), by = "sample_id") |>
  group_by(sample_id) |> mutate(frac = n / sum(n)) |> ungroup()
# guard: every sample present and histogram totals equal seqkit read counts
hist_tot <- lh |> group_by(sample_id) |> summarise(n = sum(n)) |>
  left_join(select(qc, sample_id, flnc_reads), by = "sample_id")
if (anyNA(lh$sample_id) || nrow(hist_tot) != nrow(samples) || any(hist_tot$n != hist_tot$flnc_reads)) {
  print(hist_tot); stop("length histogram check failed (missing sample or read-count mismatch)")
}
p_len <- ggplot(lh |> filter(bin_start <= 8000),
                aes(bin_start, frac, group = sample_id, color = genotype)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = GENO_COLORS) +
  labs(x = "FLNC read length (bp, 100-bp bins)", y = "fraction of reads",
       title = "Read-length distribution per sample") + theme_lr()
save_plot(p_len, file.path(d_out, "phase1_read_length_distribution"), 7, 4)

p_depth <- ggplot(qc, aes(sample_id, flnc_reads / 1e6, fill = genotype)) +
  geom_col(width = 0.7) + scale_fill_manual(values = GENO_COLORS) +
  geom_text(aes(label = sprintf("%.1f M\n%.1f%% mapped", flnc_reads / 1e6, pct_mapped)),
            vjust = -0.2, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  facet_grid(~genotype, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "FLNC reads (millions)", title = "Depth and mapping rate") + theme_lr()
save_plot(p_depth, file.path(d_out, "phase1_depth_mapping"), 7, 4)

# === U2af1 Q157R genotype check ===
message("genotype check at chr17:31867169 (T>C on + strand = CAG>CGG)")
gt <- map_dfr(samples$sample_id, \(s)
  read_tsv(file.path(d_gt, paste0(s, ".u2af1_Q157R.tsv")), show_col_types = FALSE,
           col_types = cols(.default = "c"))) |>
  mutate(across(c(depth, ref_T, alt_C, other, del, alt_C_frac), as.numeric),
         call = case_when(depth < 10 ~ "low_coverage",
                          alt_C_frac >= 0.20 ~ "Q157R",
                          alt_C_frac <= 0.02 ~ "WT",
                          TRUE ~ "ambiguous"),
         matches = call == genotype_expected)
write_tsv(gt, file.path(d_out, "phase1_genotype_check.tsv"))
print(gt, width = Inf)

p_gt <- ggplot(gt |> mutate(genotype_expected = factor(genotype_expected, GENO_LEVELS)),
               aes(sample, alt_C_frac, fill = genotype_expected)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", alt_C, depth)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = GENO_COLORS, name = "expected") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "fraction C reads (Q157R allele)",
       title = "U2af1 Q157R (chr17:31,867,169 T>C)") + theme_lr()
save_plot(p_gt, file.path(d_out, "phase1_u2af1_Q157R_genotype"), 6, 4)

if (!all(gt$matches)) {
  stop("GENOTYPE MISMATCH: ", paste(gt$sample[!gt$matches], gt$call[!gt$matches], collapse = "; "))
}
message("all 6 genotypes confirmed")
