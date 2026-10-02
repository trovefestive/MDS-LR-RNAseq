# =============================================================================
# 14_phase5_splice_signatures.R — Phase 5: genotype-associated novel isoforms + U2AF1 3'SS signatures
# Sign convention: Δ = Q157R − WT (WT = reference level, rule 8).
# Inputs : results/03_sqanti/final/{final_counts,final_isoform_annotation}.tsv.gz (R2 primary)
#          results/03_sqanti/qc/espresso_junctions.txt (SQANTI3; donor/acceptor known vs novel)
#          results/04_nmd/gene_names.tsv
# Steps  : A) specific / enriched novel isoforms per genotype
#          B) junction classes (paper: A3SS = known donor + novel acceptor, A5SS, NE, NJ) + intron retention
#          C) 3'SS bases [-3,-2,-1,+1] of A3SS junctions; Fisher tests (paper splicing_factor_tumors_signatures.R)
#          D) sequence logos (-20..+3) and ORA of genes with Q157R-associated novel isoforms
# Definitions:
#   specific_to_G : >=1 read in >=2 of 3 samples of G and 0 reads in every sample of the other genotype
#   enriched_in_G : isoform usage within its gene shifts toward G — edgeR diffSpliceDGE (isoform = "exon"),
#                   FDR < 0.05 and |Δprop| > 0.10 (paper used a Wilcoxon across many samples; impossible at n=3)
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
proj <- Sys.getenv("PROJ"); if (proj == "") stop("PROJ not set: run via scripts/submit.sh")
source(file.path(proj, "scripts/00_plot_theme.R"))
fin <- file.path(proj, "results/03_sqanti/final"); out <- file.path(proj, "results/05_splice")
dir.create(out, showWarnings = FALSE)
samples <- read_tsv(file.path(proj, "samples.tsv"), show_col_types = FALSE) |>
  mutate(genotype = factor(genotype, levels = GENO_LEVELS))
S <- samples$sample_id; wt <- S[samples$genotype == "WT"]; mu <- S[samples$genotype == "Q157R"]
design <- model.matrix(~ genotype, data = samples); stopifnot(colnames(design)[2] == "genotypeQ157R")
gn <- read_tsv(file.path(proj, "results/04_nmd/gene_names.tsv"), show_col_types = FALSE)

cnt <- read_tsv(file.path(fin, "final_counts.tsv.gz"), show_col_types = FALSE)
ann <- read_tsv(file.path(fin, "final_isoform_annotation.tsv.gz"), show_col_types = FALSE, guess_max = 2e5) |>
  mutate(gene = associated_gene, gene_ok = str_starts(gene, "ENSMUSG") & !str_detect(gene, "_"))
stopifnot(identical(cnt$transcript_ID, ann$transcript_ID))
M <- as.matrix(cnt[S]); rownames(M) <- cnt$transcript_ID

# === A1) specific isoforms (detection pattern) ===
message("== A) specific / enriched novel isoforms")
det_wt <- rowSums(M[, wt] >= 1); det_mu <- rowSums(M[, mu] >= 1)
any_wt <- rowSums(M[, wt]) > 0;  any_mu <- rowSums(M[, mu]) > 0
ann <- ann |> mutate(specific = case_when(det_mu >= 2 & !any_wt ~ "specific_to_Q157R",
                                          det_wt >= 2 & !any_mu ~ "specific_to_WT", TRUE ~ NA_character_))

# === A2) usage shift within gene: diffSpliceDGE on isoforms of genes with >=2 isoforms ===
keep_iso <- ann$gene_ok
g_iso <- ann$gene[keep_iso]
multi <- names(which(table(g_iso) >= 2))
sel <- keep_iso & ann$gene %in% multi
gtot <- rowsum(M[sel, ], ann$gene[sel])
g_keep <- rownames(gtot)[rowSums(gtot >= 10) >= 3]                     # prefilter: >=10 reads in >=3 samples
sel <- sel & ann$gene %in% g_keep & rowSums(M) > 0
y <- DGEList(M[sel, ], genes = data.frame(transcript_ID = ann$transcript_ID[sel], gene = ann$gene[sel]))
y <- calcNormFactors(y); y <- estimateDisp(y, design); fit <- glmQLFit(y, design)
ds <- diffSpliceDGE(fit, coef = 2, geneid = "gene", exonid = "transcript_ID", verbose = FALSE)
du <- topSpliceDGE(ds, test = "exon", number = Inf) |> as_tibble() |>
  transmute(transcript_ID, usage_logFC = logFC, usage_P = P.Value, usage_FDR = FDR)
prop <- M[sel, ] / gtot[ann$gene[sel], ]
dprop <- tibble(transcript_ID = rownames(prop), prop_WT = rowMeans(prop[, wt]), prop_Q157R = rowMeans(prop[, mu]),
                dprop_Q157R_minus_WT = prop_Q157R - prop_WT)
ann <- ann |> left_join(du, by = "transcript_ID") |> left_join(dprop, by = "transcript_ID") |>
  mutate(enriched = case_when(usage_FDR < 0.05 & dprop_Q157R_minus_WT >  0.10 ~ "enriched_in_Q157R",
                              usage_FDR < 0.05 & dprop_Q157R_minus_WT < -0.10 ~ "enriched_in_WT", TRUE ~ NA_character_),
         assoc = coalesce(specific, enriched),
         assoc = str_replace(assoc, "^(specific_to|enriched_in)_", ""))     # "Q157R" / "WT" / NA

# === rule-8 checks: mutant-specific isoforms with usage tested must have Δprop > 0 and usage logFC > 0 ===
chk <- ann |> filter(specific == "specific_to_Q157R", !is.na(usage_logFC))
message("Q157R-specific isoforms with usage tested: ", nrow(chk), "; NA dprop: ", sum(is.na(chk$dprop_Q157R_minus_WT)),
        "; NA logFC: ", sum(is.na(chk$usage_logFC)))
if (any(chk$dprop_Q157R_minus_WT <= 0, na.rm = TRUE) || any(chk$usage_logFC <= 0, na.rm = TRUE))
  stop("RULE 8 VIOLATION: Q157R-specific isoform with non-positive Δprop/logFC")
sig <- ann |> filter(!is.na(enriched))
if (any(sign(sig$usage_logFC) != sign(sig$dprop_Q157R_minus_WT), na.rm = TRUE)) stop("RULE 8 VIOLATION: usage logFC vs Δprop sign")
message("rule 8 OK (", nrow(chk), " Q157R-specific isoforms with usage tested; ", nrow(sig), " enriched calls)")
u <- ann |> filter(transcript_ID == "ENSMUST00000468653.1")
message("U2af1-213: specific=", u$specific, " Δprop=", round(u$dprop_Q157R_minus_WT, 3), " FDR=", signif(u$usage_FDR, 3))

assoc_tab <- ann |> filter(!is.na(specific) | !is.na(enriched)) |>
  left_join(gn, by = c("gene" = "gene_id")) |>
  select(transcript_ID, gene, gene_name, origin, structural_category, subcategory, coding, predicted_NMD,
         specific, enriched, assoc, prop_WT, prop_Q157R, dprop_Q157R_minus_WT, usage_logFC, usage_FDR) |>
  left_join(select(cnt, transcript_ID, all_of(S)), by = "transcript_ID")      # per-sample read counts
write_tsv(assoc_tab, file.path(out, "phase5_genotype_associated_isoforms.tsv"))
print(count(assoc_tab, origin, specific, enriched), n = 30)

# === B) junction classes (SQANTI3 donor/acceptor novelty, strand-aware) ===
message("== B) junction classes")
suppressPackageStartupMessages({ library(BSgenome.Mmusculus.UCSC.mm39); library(Biostrings) })
jx <- read_tsv(file.path(proj, "results/03_sqanti/qc/espresso_junctions.txt"), show_col_types = FALSE,
               col_select = c(isoform, chrom, strand, genomic_start_coord, genomic_end_coord, junction_category,
                              start_site_category, end_site_category, splice_site, canonical)) |>
  filter(isoform %in% ann$transcript_ID) |>
  mutate(donor_cat    = if_else(strand == "+", start_site_category, end_site_category),
         acceptor_cat = if_else(strand == "+", end_site_category,   start_site_category),
         jtype = case_when(donor_cat == "known" & acceptor_cat == "novel" ~ "A3SS",
                           donor_cat == "novel" & acceptor_cat == "known" ~ "A5SS",
                           donor_cat == "novel" & acceptor_cat == "novel" ~ "NE",
                           junction_category == "novel"                  ~ "NJ",     # known sites, new pairing
                           TRUE                                          ~ "known"),
         acc_pos = if_else(strand == "+", genomic_end_coord, genomic_start_coord),   # last intron base at 3'SS
         jid = paste(chrom, genomic_start_coord, genomic_end_coord, strand, sep = ":")) |>
  left_join(select(ann, transcript_ID, origin, assoc, subcategory), by = c("isoform" = "transcript_ID"))
ev <- jx |> filter(origin == "novel", !is.na(assoc), jtype != "known") |> distinct(jid, assoc, jtype) |>
  count(assoc, jtype, name = "n_unique_junctions")
ir <- ann |> filter(origin == "novel", !is.na(assoc)) |> group_by(assoc) |>
  summarise(n_isoforms = n(), n_intron_retention = sum(subcategory == "intron_retention", na.rm = TRUE))
write_tsv(ev, file.path(out, "phase5_event_types.tsv")); write_tsv(ir, file.path(out, "phase5_intron_retention.tsv"))
print(ev); print(ir)

# === C) 3'SS sequences (-20..+3) and U2AF1 position tests ===
message("== C) 3'SS motif tests")
ss_seq <- function(d) {                       # one row per unique acceptor; returns 23-nt strings, 5'->3' of RNA
  gr <- GRanges(d$chrom, IRanges(start = if_else(d$strand == "+", d$acc_pos - 19L, d$acc_pos - 3L), width = 23),
                strand = d$strand)
  as.character(getSeq(BSgenome.Mmusculus.UCSC.mm39, gr))        # getSeq reverse-complements "-" ranges
}
acc_sets <- list(
  Q157R      = jx |> filter(origin == "novel", assoc == "Q157R", jtype == "A3SS"),
  WT         = jx |> filter(origin == "novel", assoc == "WT",    jtype == "A3SS"),
  all_novel  = jx |> filter(origin == "novel", jtype == "A3SS"),
  annotated  = jx |> filter(jtype == "known")) |>
  map(\(d) d |> distinct(chrom, strand, acc_pos) |> mutate(seq = ss_seq(pick(everything()))))
# positions in the 23-mer: 1..20 = -20..-1 (18,19,20 = -3,-2,-1), 21..23 = +1..+3
base_at <- function(s, i) substr(s, i, i)
pos_tab <- imap_dfr(acc_sets, \(d, nm) tibble(set = nm, n = nrow(d),
  pct_AG = 100 * mean(substr(d$seq, 19, 20) == "AG"),
  m3_C = sum(base_at(d$seq, 18) == "C"), m3_T = sum(base_at(d$seq, 18) == "T"),
  p1_G = sum(base_at(d$seq, 21) == "G"), p1_A = sum(base_at(d$seq, 21) == "A")))
print(pos_tab)
if (any(pos_tab$pct_AG[pos_tab$n >= 20] < 95)) stop("3'SS AG sanity check failed (<95% AG): check acceptor coordinates")
fish <- function(a, b, x, y) {               # set a vs set b; x vs y counts (e.g. G vs A)
  ta <- pos_tab[pos_tab$set == a, ]; tb <- pos_tab[pos_tab$set == b, ]
  m <- matrix(c(ta[[x]], ta[[y]], tb[[x]], tb[[y]]), nrow = 2)
  if (any(colSums(m) < 2)) return(tibble(OR = NA_real_, p = NA_real_))
  f <- fisher.test(m); tibble(OR = unname(f$estimate), p = f$p.value)
}
tests <- crossing(comparison = c("Q157R vs WT", "Q157R vs annotated", "WT vs annotated", "all_novel vs annotated"),
                  signature = c("Q157 +1 G vs A", "S34 -3 C vs T (negative control)")) |>
  mutate(a = word(comparison, 1), b = word(comparison, 3),
         x = if_else(str_starts(signature, "Q157"), "p1_G", "m3_C"),
         y = if_else(str_starts(signature, "Q157"), "p1_A", "m3_T")) |>
  mutate(res = pmap(list(a, b, x, y), fish)) |> unnest(res) |>
  left_join(transmute(pos_tab, a = set, n_a = n, frac_x_a = NA_real_), by = "a") |>
  rowwise() |>
  mutate(frac_x_a = pos_tab[[x]][pos_tab$set == a] / (pos_tab[[x]][pos_tab$set == a] + pos_tab[[y]][pos_tab$set == a]),
         frac_x_b = pos_tab[[x]][pos_tab$set == b] / (pos_tab[[x]][pos_tab$set == b] + pos_tab[[y]][pos_tab$set == b])) |>
  ungroup() |> select(comparison, signature, n_a, frac_x_a, frac_x_b, OR, p)
write_tsv(pos_tab, file.path(out, "phase5_3ss_position_counts.tsv")); write_tsv(tests, file.path(out, "phase5_3ss_fisher_tests.tsv"))
print(tests, width = Inf)

# === D1) sequence logos around novel A3SS acceptors (-20..+3) ===
message("== D) logos + ORA")
suppressPackageStartupMessages(library(ggseqlogo))
set.seed(1)
logo_sets <- list(`Q157R-associated novel A3SS` = acc_sets$Q157R$seq, `WT-associated novel A3SS` = acc_sets$WT$seq,
                  `annotated 3'SS (5,000 sampled)` = sample(acc_sets$annotated$seq, min(5000, nrow(acc_sets$annotated))))
logo_sets <- logo_sets[lengths(logo_sets) >= 10]
names(logo_sets) <- paste0(names(logo_sets), " (n=", lengths(logo_sets), ")")
pl <- ggseqlogo(logo_sets, method = "prob", ncol = 1) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 18, 19, 20, 21, 23), labels = c(-20, -16, -11, -6, -3, -2, -1, "+1", "+3")) +
  labs(x = "position relative to 3' splice site (AG|)", title = "3' splice site composition") + theme_lr(9)
save_plot(pl, file.path(out, "phase5_3ss_logos"), 8, 2.2 * length(logo_sets) + 1)

# === D2) ORA: genes with Q157R-associated novel isoforms (GO BP; universe = testable genes) ===
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
strip_v <- function(x) sub("\\.\\d+$", "", x)
universe <- unique(strip_v(c(g_keep, ann$gene[!is.na(ann$specific) & ann$gene_ok])))
ora_genes <- function(g) unique(strip_v(ann$gene[ann$origin == "novel" & ann$assoc %in% g & ann$gene_ok]))
ora <- map_dfr(c("Q157R", "WT"), \(g) {
  genes <- ora_genes(g); message("ORA ", g, ": ", length(genes), " genes")
  if (length(genes) < 10) return(tibble())
  e <- enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP", universe = universe,
                pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  as_tibble(e@result) |> mutate(group = g, .before = 1)
})
if (nrow(ora)) {
  write_tsv(ora, file.path(out, "phase5_ora_go_bp.tsv"))
  print(ora |> filter(p.adjust < 0.05) |> count(group, name = "n_GO_BP_FDR<0.05"))
  print(ora |> group_by(group) |> slice_min(p.adjust, n = 8, with_ties = FALSE) |> ungroup() |>
          dplyr::select(group, Description, GeneRatio, BgRatio, p.adjust), n = 20, width = Inf)
}
message("phase 5 done")
