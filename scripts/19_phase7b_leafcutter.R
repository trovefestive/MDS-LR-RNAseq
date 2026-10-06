# =============================================================================
# 19_phase7b_leafcutter.R — Phase 7b: LeafCutter results, short reads (LSK) and long reads (LK)
# deltapsi = Q157R − WT (groups file lists WT first = baseline; checked below on column names and U2af1-213).
# Inputs : results/07_leafcutter/{shortread,longread}/<set>_ds_{cluster_significance,effect_sizes}.txt (07e)
# Steps  : 1) import + coordinate check (GT-AG) + rule-8 checks   2) significant clusters / introns
#          3) U2AF1 Q157 test: +1 base after the AG at Q157R-favoured vs disfavoured acceptors sharing a donor;
#             −3 C/T as negative control   4) short vs long reads; long reads vs LRP2 (negated); rMATS-turbo;
#             Phase 8 tier A genes (Cd34)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(BSgenome.Mmusculus.UCSC.mm39) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
lc <- file.path(proj, "results/07_leafcutter"); out <- file.path(lc, "compare"); dir.create(out, showWarnings = FALSE)
G <- BSgenome.Mmusculus.UCSC.mm39
base_at <- function(chr, pos, strand) as.character(getSeq(G, GRanges(chr, IRanges(pos, width = 1), strand = strand)))
dinuc   <- function(chr, pos, strand) as.character(getSeq(G, GRanges(chr, IRanges(pos, width = 2), strand = strand)))

# === 1) import ===
read_set <- function(set) {
  d <- file.path(lc, set)
  cs <- read_tsv(file.path(d, paste0(set, "_ds_cluster_significance.txt")), show_col_types = FALSE)
  es <- read_tsv(file.path(d, paste0(set, "_ds_effect_sizes.txt")), show_col_types = FALSE)
  # rule 8: baseline column must be WT, perturbed Q157R, deltapsi = Q157R − WT
  if (!identical(colnames(es)[3:4], c("WT", "Q157R"))) stop("RULE 8: ", set, " effect-size columns are ", paste(colnames(es)[3:4], collapse = ","))
  stopifnot(max(abs(es$deltapsi - (es$Q157R - es$WT))) < 1e-6)
  es <- es |> separate(intron, c("chr", "A", "B", "clu"), sep = ":", remove = FALSE, convert = TRUE) |>
    mutate(strand = str_sub(clu, -1), cluster = paste0(chr, ":", clu))
  list(cs = cs, es = es |> left_join(select(cs, cluster, status, p.adjust, genes), by = "cluster"))
}
S <- list(shortread = read_set("shortread"), longread = read_set("longread"))

# coordinates: regtools clustering writes A = 0-based first intron base, B = 1-based first base of the next exon,
# i.e. intron (1-based) = A+1 .. B-1. Verified on GT…AG (strand-aware) before any sequence test.
for (set in names(S)) {
  e <- S[[set]]$es |> filter(strand %in% c("+", "-")) |> slice_sample(n = 5000)
  don <- if_else(e$strand == "+", e$A + 1L, e$B - 2L); acc <- if_else(e$strand == "+", e$B - 2L, e$A + 1L)
  gtag <- mean(dinuc(e$chr, don, e$strand) == "GT" & dinuc(e$chr, acc, e$strand) == "AG")
  message(set, ": GT-AG fraction under intron = A+1..B-1: ", round(gtag, 3))
  if (gtag < 0.9) stop("intron coordinate convention check failed for ", set)
}
S <- map(S, \(x) { x$es <- x$es |> mutate(istart = A + 1L, iend = B - 1L); x })

# positive control: U2af1-213-specific intron (created by the Q157R base) must have deltapsi > 0 where tested
gtf <- read_tsv(file.path(proj, "data/ref/expanded/expanded_final.gtf"), col_names = FALSE, comment = "#", show_col_types = FALSE) |>
  filter(X3 == "exon", str_detect(X9, 'gene_id "ENSMUSG00000061613')) |>
  transmute(chr = X1, start = X4, end = X5, strand = X7, tx = str_match(X9, 'transcript_id "([^"]+)"')[, 2]) |>
  arrange(tx, start) |> group_by(tx) |> transmute(chr, strand, istart = lag(end) + 1L, iend = start - 1L) |> filter(!is.na(istart)) |> ungroup()
# The Q157R base (chr17:31,867,169, − strand) creates the 5′SS of the U2af1-213 intron chr17:31,867,070–31,867,168
# (also annotated in ENSMUST00000563584 and used by novel ESPRESSO isoforms): that intron is the positive control.
u213 <- gtf |> filter(tx == "ENSMUST00000468653.1", chr == "chr17", abs(iend - 31867168L) <= 2L)
stopifnot(nrow(u213) == 1)
message("U2af1-213-specific introns: ", nrow(u213))
ctrl <- map_dfr(names(S), \(set) S[[set]]$es |> semi_join(u213, by = c("chr", "istart", "iend")) |>
                  mutate(set = set) |> select(set, intron, WT, Q157R, deltapsi, status, p.adjust))
print(ctrl, width = Inf); write_tsv(ctrl, file.path(out, "phase7b_U2af1_213_control.tsv"))
if (any(ctrl$deltapsi <= 0)) stop("RULE 8 VIOLATION: U2af1-213 junction deltapsi not > 0")
if (!nrow(ctrl)) message("U2af1-213 junction not in a tested cluster: rule 8 rests on the column-order check")

# === 2) significant clusters and introns ===
summ <- map_dfr(names(S), \(set) {
  cs <- S[[set]]$cs; es <- S[[set]]$es
  sig_cl <- cs$cluster[cs$status == "Success" & !is.na(cs$p.adjust) & cs$p.adjust < 0.05]
  si <- es |> filter(cluster %in% sig_cl, abs(deltapsi) >= 0.10)
  tibble(set, clusters_tested = sum(cs$status == "Success"), clusters_sig = length(sig_cl),
         genes_sig = n_distinct(na.omit(unlist(str_split(cs$genes[cs$cluster %in% sig_cl], ",")))),
         introns_sig_dpsi10 = nrow(si), up_in_Q157R = sum(si$deltapsi > 0))
})
write_tsv(summ, file.path(out, "phase7b_leafcutter_summary.tsv")); print(summ, width = Inf)
for (set in names(S)) write_tsv(S[[set]]$es |> filter(!is.na(p.adjust), p.adjust < 0.05) |> arrange(p.adjust, desc(abs(deltapsi))) |>
                                  select(cluster, genes, intron, chr, istart, iend, strand, WT, Q157R, deltapsi_Q157R_minus_WT = deltapsi, p.adjust),
                                file.path(out, paste0("phase7b_", set, "_significant_introns.tsv")))

# === 3) U2AF1 Q157 signature: alternative acceptors sharing a donor, inside significant clusters ===
# For each donor with ≥ 2 introns in a significant cluster: Q157R-favoured acceptor = max deltapsi (> 0),
# disfavoured = min deltapsi (< 0). +1 = first exonic base after the AG; −3 = intronic base before the AG.
sig3 <- map_dfr(names(S), \(set) {
  es <- S[[set]]$es |> filter(!is.na(p.adjust), p.adjust < 0.05, strand %in% c("+", "-")) |>
    mutate(donor = if_else(strand == "+", istart, iend), acc = if_else(strand == "+", iend, istart))
  pr <- es |> group_by(cluster, donor) |> filter(n_distinct(acc) >= 2, max(deltapsi) > 0, min(deltapsi) < 0) |>
    summarise(chr = chr[1], strand = strand[1], fav = acc[which.max(deltapsi)], dis = acc[which.min(deltapsi)],
              d_fav = max(deltapsi), d_dis = min(deltapsi), genes = genes[1], .groups = "drop")
  if (!nrow(pr)) return(tibble())
  p1 <- function(a, s) if_else(s == "+", a + 1L, a - 1L)          # first exonic base after the 3′SS
  m3 <- function(a, s) if_else(s == "+", a - 2L, a + 2L)          # −3: base before the AG
  pr |> mutate(set = set, fav_p1 = base_at(chr, p1(fav, strand), strand), dis_p1 = base_at(chr, p1(dis, strand), strand),
               fav_m3 = base_at(chr, m3(fav, strand), strand), dis_m3 = base_at(chr, m3(dis, strand), strand),
               fav_ag = dinuc(chr, if_else(strand == "+", fav - 1L, fav), strand),
               dis_ag = dinuc(chr, if_else(strand == "+", dis - 1L, dis), strand))
})
write_tsv(sig3, file.path(out, "phase7b_alt_acceptor_pairs.tsv"))
fisher_row <- function(d, set, test, b1, b2, colf, cold) {
  m <- matrix(c(sum(d[[colf]] == b1), sum(d[[colf]] == b2), sum(d[[cold]] == b1), sum(d[[cold]] == b2)), 2)
  ft <- if (all(colSums(m) > 0)) fisher.test(m) else list(estimate = NA, p.value = NA)
  tibble(set, test, n_pairs = nrow(d), fav_frac = m[1, 1] / sum(m[, 1]), dis_frac = m[1, 2] / sum(m[, 2]),
         OR = unname(ft$estimate), p = ft$p.value)
}
sig_tests <- map_dfr(unique(sig3$set), \(set) {
  d <- filter(sig3, set == !!set)
  message(set, ": AG at favoured / disfavoured acceptors ", round(100 * mean(d$fav_ag == "AG"), 1), "% / ",
          round(100 * mean(d$dis_ag == "AG"), 1), "%")
  bind_rows(fisher_row(d, set, "+1 G vs A (Q157 signature)", "G", "A", "fav_p1", "dis_p1"),
            fisher_row(d, set, "-3 C vs T (S34F control)", "C", "T", "fav_m3", "dis_m3"))
})
write_tsv(sig_tests, file.path(out, "phase7b_3ss_tests.tsv")); print(sig_tests, width = Inf)

# === 4) comparisons ===
key <- \(e) e |> mutate(k = paste(chr, istart, iend, strand, sep = ":"))
sr <- key(S$shortread$es); lr <- key(S$longread$es)
both <- inner_join(select(sr, k, sr_dpsi = deltapsi, sr_padj = p.adjust, genes),
                   select(lr, k, lr_dpsi = deltapsi, lr_padj = p.adjust), by = "k")
sig_either <- both |> filter((!is.na(sr_padj) & sr_padj < 0.05 & abs(sr_dpsi) >= 0.1) | (!is.na(lr_padj) & lr_padj < 0.05 & abs(lr_dpsi) >= 0.1))
se <- sig_either   # keep the data frame out of tibble()'s sequential scope
cmp <- tibble(shared_introns = nrow(both), r_all = cor(both$sr_dpsi, both$lr_dpsi),
              n_sig_either = nrow(se), r_sig_either = cor(se$sr_dpsi, se$lr_dpsi),
              same_sign_sig_either = mean(sign(se$sr_dpsi) == sign(se$lr_dpsi)),
              n_sig_both = sum(se$sr_padj < 0.05 & se$lr_padj < 0.05, na.rm = TRUE))
write_tsv(cmp, file.path(out, "phase7b_short_vs_long.tsv")); print(cmp, width = Inf)

# long reads vs the earlier LRP2 long-read LeafCutter (its deltapsi_WT = WT − Q157R -> negated)
old <- file.path(proj, "2026_04/260327_MDS_lrp2/S4_MULTISAMPLE_ANALYSIS/M1_LEAFCUTTER_LONGREAD/Q157R_vs_WT.lr_leafcutter.ds_effect_sizes.txt")
if (file.exists(old)) {
  o <- read_tsv(old, show_col_types = FALSE) |> transmute(intron, old_dpsi = -deltapsi_WT) |>
    separate(intron, c("chr", "A", "B", "clu"), sep = ":") |> mutate(A = as.integer(A), B = as.integer(B)) |> filter(!is.na(A))
  m <- inner_join(S$longread$es, o, by = c("chr", "A", "B"))
  lrp2 <- tibble(shared_introns = nrow(m), r = if (nrow(m) > 2) cor(m$deltapsi, m$old_dpsi) else NA_real_)
  write_tsv(lrp2, file.path(out, "phase7b_longread_vs_LRP2.tsv")); message("long reads vs LRP2 (negated):"); print(lrp2)
}

# rMATS-turbo (same short-read BAMs): genes with significant events vs genes with significant LeafCutter clusters
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)
rm_genes <- map(c("SE", "A3SS", "A5SS", "RI", "MXE"), \(e)
  read_tsv(file.path(proj, "results/07_shortread/rmats/out", paste0(e, ".MATS.JCEC.txt")), show_col_types = FALSE) |>
    filter(FDR < 0.05, abs(IncLevelDifference) >= 0.1) |> pull(GeneID)) |> unlist() |> str_remove_all('"') |> unique()
rm_names <- gn$gene_name[match(rm_genes, gn$gene_id)] |> na.omit() |> unique()
sig_names <- \(set) { cs <- S[[set]]$cs; unique(na.omit(unlist(str_split(cs$genes[cs$status == "Success" & cs$p.adjust < 0.05], ",")))) }
ov <- tibble(set = names(S), leafcutter_genes = map_int(set, \(s) length(sig_names(s))),
             also_rmats_turbo = map_int(set, \(s) length(intersect(sig_names(s), rm_names))), rmats_turbo_genes = length(rm_names))
write_tsv(ov, file.path(out, "phase7b_vs_rmats_turbo_genes.tsv")); print(ov)

# Phase 8 tier A genes: significant LeafCutter cluster in either dataset?
tA <- read_tsv(file.path(proj, "results/08_summary/phase8_candidates_genes.tsv"), show_col_types = FALSE) |> filter(tier == "A")
tA_lc <- tA |> transmute(gene_name, lr_delta = delta, sr_delta, leafcutter_SR = gene_name %in% sig_names("shortread"),
                         leafcutter_LR = gene_name %in% sig_names("longread"))
write_tsv(tA_lc, file.path(out, "phase7b_tierA_leafcutter.tsv"))
message("tier A genes with a significant LeafCutter cluster: SR ", sum(tA_lc$leafcutter_SR), ", LR ", sum(tA_lc$leafcutter_LR), " of ", nrow(tA_lc))
cd34 <- bind_rows(sr |> mutate(set = "shortread"), lr |> mutate(set = "longread")) |> filter(str_detect(coalesce(genes, ""), "\\bCd34\\b")) |>
  select(set, cluster, intron, istart, iend, WT, Q157R, deltapsi, p.adjust)
print(cd34, n = 40, width = Inf); write_tsv(cd34, file.path(out, "phase7b_Cd34_introns.tsv"))

# === 5) figure: short vs long ΔPSI for shared introns ===
p <- ggplot(both, aes(lr_dpsi, sr_dpsi)) + geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(colour = "grey75", size = 0.5, alpha = 0.5) +
  geom_point(data = sig_either, colour = "#C8553D", size = 1.2) +
  labs(x = "Long reads (LK): ΔPSI, Q157R − WT", y = "Short reads (LSK): ΔPSI, Q157R − WT",
       title = "LeafCutter intron usage: short vs long reads",
       subtitle = sprintf("%s shared introns (r = %.2f); red: significant in either (r = %.2f)", format(nrow(both), big.mark = ","),
                          cmp$r_all, cmp$r_sig_either)) + theme_lr()
save_plot(p, file.path(out, "phase7b_short_vs_long_dpsi"), 6, 5.5)
message("phase 7b LeafCutter summary done")
