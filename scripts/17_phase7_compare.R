# =============================================================================
# 17_phase7_compare.R — Phase 7: short-read (in-house, PE150, WT1-3 vs KI1-3 = Q157R) vs long-read results
# All effects Q157R − WT (WT = reference). Different animals: this is replication, not matched validation.
# Inputs : results/07_shortread/{genotype,salmon,rmats/out}; long-read results from Phases 3, 5, 6
# Steps  : 1) genotype confirmation   2) detection of novel long-read isoforms in short reads (Salmon)
#          3) gene-level DE replication (edgeR on Salmon gene sums)
#          4) isoform-usage replication for the long-read high-confidence switches (rMATS-long ∩ DRIMSeq)
#          5) rMATS-turbo: event counts, and the U2AF1 Q157 +1 G signature at Q157R-favoured A3SS acceptors
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
o <- file.path(proj, "results/07_shortread"); out <- file.path(o, "compare"); dir.create(out, showWarnings = FALSE)
ss <- read_tsv(file.path(proj, "samples_shortread.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))
S <- ss$sample_id; wt <- S[ss$genotype == "WT"]; mu <- S[ss$genotype == "Q157R"]
design <- model.matrix(~ genotype, data = ss); stopifnot(colnames(design)[2] == "genotypeQ157R")
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)
ann <- read_tsv(file.path(proj, "results/03_sqanti/final/final_isoform_annotation.tsv.gz"), show_col_types = FALSE,
                guess_max = 2e5) |> mutate(gene = associated_gene)

# === 1) genotype ===
gt <- map_dfr(S, \(s) read_tsv(file.path(o, "genotype", paste0(s, ".u2af1_Q157R.tsv")), show_col_types = FALSE)) |>
  mutate(call = case_when(depth < 10 ~ "low_coverage", alt_C_frac >= 0.2 ~ "Q157R", alt_C_frac <= 0.02 ~ "WT", TRUE ~ "ambiguous"))
write_tsv(gt, file.path(out, "phase7_genotype_check.tsv")); print(gt)
if (any(gt$call != gt$genotype_expected)) stop("short-read genotype mismatch: ", paste(gt$sample[gt$call != gt$genotype_expected], collapse = ","))

# === 2) Salmon: detection of long-read isoforms ===
q <- map(set_names(S), \(s) read_tsv(file.path(o, "salmon", s, "quant.sf"), show_col_types = FALSE))
cnt <- sapply(q, \(x) x$NumReads); rownames(cnt) <- q[[1]]$Name
tpm <- sapply(q, \(x) x$TPM);      rownames(tpm) <- q[[1]]$Name
# Salmon collapses transcripts with identical spliced sequence into one (index duplicate_clusters.tsv): the retained
# isoform carries all reads. Collapsed isoforms get 0-count rows (gene sums unchanged) and are excluded from detection;
# both members of a collapsed pair are excluded from isoform-usage comparisons.
dup <- read_tsv(file.path(proj, "data/ref/expanded/salmon_index/duplicate_clusters.tsv"), show_col_types = FALSE)
seq_dup <- union(dup$RetainedRef, dup$DuplicateRef)
z <- matrix(0, length(dup$DuplicateRef), ncol(cnt), dimnames = list(dup$DuplicateRef, colnames(cnt)))
cnt <- rbind(cnt, z); tpm <- rbind(tpm, z)
stopifnot(all(ann$transcript_ID %in% rownames(cnt)))
message("Salmon: ", nrow(dup), " isoforms collapsed into an identical-sequence isoform (",
        sum(str_starts(dup$DuplicateRef, "ESPRESSO")), " novel); excluded from detection")
det <- ann |> filter(!transcript_ID %in% dup$DuplicateRef) |> transmute(transcript_ID, origin, structural_category,
                        detected = rowSums(cnt[transcript_ID, ] >= 1) >= 2 & rowSums(cnt[transcript_ID, ]) >= 5)
det_tab <- det |> group_by(origin) |> summarise(n = n(), pct = 100 * mean(detected), detected = sum(detected)) |> relocate(pct, .after = detected)
write_tsv(det_tab, file.path(out, "phase7_isoform_detection.tsv")); print(det_tab)
mi <- map_dfr(S, \(s) jsonlite::fromJSON(file.path(o, "salmon", s, "aux_info/meta_info.json"))[c("num_processed", "percent_mapped")] |>
                as_tibble() |> mutate(sample_id = s, lib = jsonlite::fromJSON(file.path(o, "salmon", s, "lib_format_counts.json"))$expected_format))
write_tsv(mi, file.path(out, "phase7_salmon_mapping.tsv")); print(mi)

# === 3) gene-level DE replication ===
ok <- str_starts(ann$gene, "ENSMUSG") & !str_detect(ann$gene, "_")
g_sr <- rowsum(round(cnt[ann$transcript_ID[ok], ]), ann$gene[ok])
y <- DGEList(g_sr); y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
y <- estimateDisp(y, design); f <- glmQLFit(y, design)
dge_sr <- topTags(glmQLFTest(f, coef = 2), n = Inf)$table |> rownames_to_column("gene") |>
  transmute(gene, sr_logFC = logFC, sr_FDR = FDR)
dge_lr <- read_tsv(file.path(proj, "results/06_diff/phase6_DGE_edgeR_Q157R_vs_WT.tsv"), show_col_types = FALSE) |>
  transmute(gene, gene_name, lr_logFC = logFC_Q157R_vs_WT, lr_FDR = FDR)
mg <- inner_join(dge_lr, dge_sr, by = "gene")
lrs <- mg |> filter(lr_FDR < 0.05)
dge_cmp <- tibble(n_shared = nrow(mg), r_all = cor(mg$lr_logFC, mg$sr_logFC),
                  lr_sig = nrow(lrs), lr_sig_same_sign_in_sr = sum(sign(lrs$lr_logFC) == sign(lrs$sr_logFC)),
                  lr_sig_also_sr_FDR05 = sum(lrs$sr_FDR < 0.05 & sign(lrs$lr_logFC) == sign(lrs$sr_logFC)),
                  sr_sig = sum(mg$sr_FDR < 0.05))
write_tsv(mg, file.path(out, "phase7_gene_DE_lr_vs_sr.tsv")); write_tsv(dge_cmp, file.path(out, "phase7_gene_DE_summary.tsv"))
print(dge_cmp, width = Inf)
ribo <- mg |> filter(str_detect(gene_name, "^Rp[sl]\\d"), lr_FDR < 0.05)
message("ribosomal genes DE in long reads: ", nrow(ribo), "; same direction in short reads: ", sum(sign(ribo$sr_logFC) == sign(ribo$lr_logFC)),
        "; median short-read logFC ", round(median(ribo$sr_logFC), 2))

# === 4) isoform-usage replication: long-read high-confidence switches (rMATS-long ∩ DRIMSeq) ===
rl  <- read_tsv(file.path(proj, "results/06_diff/rmats_long/differential_transcripts.tsv"), show_col_types = FALSE)
dtu <- read_tsv(file.path(proj, "results/06_diff/phase6_DTU_DRIMSeq_Q157R_vs_WT.tsv"), show_col_types = FALSE)
hc <- inner_join(rl |> filter(adj_pvalue < 0.05, abs(delta_isoform_proportion) >= 0.10) |>
                   transmute(transcript_ID = feature_id, gene_id, lr_delta = delta_isoform_proportion),
                 dtu |> filter(DTU_tx_FDR < 0.05, abs(dprop_Q157R_minus_WT) > 0.10) |> select(transcript_ID, gene_name),
                 by = "transcript_ID")
gid <- set_names(rl$gene_id, rl$feature_id)                                   # rMATS-long gene grouping
iso <- names(gid)[gid %in% hc$gene_id]
gsum <- rowsum(cnt[iso, ], gid[iso])
prop <- cnt[hc$transcript_ID, , drop = FALSE] / gsum[hc$gene_id, , drop = FALSE]
hc <- hc |> mutate(sr_prop_WT = rowMeans(prop[, wt, drop = FALSE]), sr_prop_Q157R = rowMeans(prop[, mu, drop = FALSE]),
                   sr_delta = if_else(transcript_ID %in% seq_dup, NA_real_, sr_prop_Q157R - sr_prop_WT), sr_gene_reads = rowMeans(gsum[gene_id, , drop = FALSE]),
                   same_sign = sign(sr_delta) == sign(lr_delta))
usage_cmp <- tibble(n_high_conf = nrow(hc), testable_sr = sum(!is.na(hc$sr_delta) & hc$sr_gene_reads >= 10),
                    same_sign = sum(hc$same_sign & hc$sr_gene_reads >= 10, na.rm = TRUE),
                    r_delta = cor(hc$lr_delta, hc$sr_delta, use = "complete.obs"),
                    binom_p_same_sign = binom.test(sum(hc$same_sign & hc$sr_gene_reads >= 10, na.rm = TRUE),
                                                   sum(!is.na(hc$sr_delta) & hc$sr_gene_reads >= 10), 0.5)$p.value)
write_tsv(hc, file.path(out, "phase7_usage_highconf_lr_vs_sr.tsv")); write_tsv(usage_cmp, file.path(out, "phase7_usage_summary.tsv"))
print(usage_cmp, width = Inf); print(hc |> filter(gene_name %in% c("Cd34", "Set", "Tmpo", "Usp1", "Dock10", "Atrx")), width = Inf)

# === 5) rMATS-turbo (short reads): events + Q157 +1 G at Q157R-favoured vs disfavoured A3SS acceptors ===
rm_dir <- file.path(o, "rmats/out")
ev <- map_dfr(c("SE", "A3SS", "A5SS", "RI", "MXE"), \(e) {
  x <- read_tsv(file.path(rm_dir, paste0(e, ".MATS.JCEC.txt")), show_col_types = FALSE)
  tibble(event = e, n = nrow(x), sig = sum(x$FDR < 0.05 & abs(x$IncLevelDifference) >= 0.1, na.rm = TRUE),
         sig_up_inclusion = sum(x$FDR < 0.05 & x$IncLevelDifference >= 0.1, na.rm = TRUE))
})
write_tsv(ev, file.path(out, "phase7_rmats_turbo_events.tsv")); print(ev)
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm39))
a3 <- read_tsv(file.path(rm_dir, "A3SS.MATS.JCEC.txt"), show_col_types = FALSE) |>
  filter(FDR < 0.05, abs(IncLevelDifference) >= 0.1)
# 1-based first exonic base (+1) of each acceptor; long exon = inclusion form. + strand: acceptor at exon start;
# - strand: acceptor at exon end. ΔPSI > 0 -> long-exon acceptor favoured in Q157R.
a3 <- a3 |> mutate(long_acc  = if_else(strand == "+", longExonStart_0base + 1, longExonEnd),
                   short_acc = if_else(strand == "+", shortES + 1, shortEE),
                   fav = if_else(IncLevelDifference > 0, long_acc, short_acc),
                   dis = if_else(IncLevelDifference > 0, short_acc, long_acc))
b1 <- function(chr, pos, strand) as.character(getSeq(BSgenome.Mmusculus.UCSC.mm39, GRanges(chr, IRanges(pos, width = 1), strand = strand)))
ag <- function(chr, pos, strand) as.character(getSeq(BSgenome.Mmusculus.UCSC.mm39,       # 2 intronic bases before +1
        GRanges(chr, IRanges(if_else(strand == "+", pos - 2L, pos + 1L), width = 2), strand = strand)))
if (nrow(a3)) {
  a3 <- a3 |> mutate(fav_p1 = b1(chr, fav, strand), dis_p1 = b1(chr, dis, strand),
                     fav_ag = ag(chr, fav, strand), dis_ag = ag(chr, dis, strand))
  message("A3SS sanity: AG before favoured/disfavoured acceptors: ", round(100 * mean(a3$fav_ag == "AG"), 1), "% / ",
          round(100 * mean(a3$dis_ag == "AG"), 1), "%")
  m <- matrix(c(sum(a3$fav_p1 == "G"), sum(a3$fav_p1 == "A"), sum(a3$dis_p1 == "G"), sum(a3$dis_p1 == "A")), 2,
              dimnames = list(c("G", "A"), c("Q157R_favoured", "Q157R_disfavoured")))
  ft <- fisher.test(m); print(m)
  sig_tab <- tibble(n_sig_A3SS = nrow(a3), favoured_G_frac = m["G", 1] / sum(m[, 1]), disfavoured_G_frac = m["G", 2] / sum(m[, 2]),
                    fisher_OR = unname(ft$estimate), fisher_p = ft$p.value)
  write_tsv(sig_tab, file.path(out, "phase7_rmats_A3SS_plus1G.tsv")); write_tsv(a3, file.path(out, "phase7_rmats_A3SS_significant.tsv"))
  print(sig_tab, width = Inf)
}
cd34 <- read_tsv(file.path(rm_dir, "A3SS.MATS.JCEC.txt"), show_col_types = FALSE) |> filter(geneSymbol %in% "Cd34" | str_detect(GeneID, "ENSMUSG00000016494"))
print(cd34 |> select(any_of(c("geneSymbol", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "IncLevel1", "IncLevel2",
                              "IncLevelDifference", "FDR"))), width = Inf)
message("phase 7 comparison done")
