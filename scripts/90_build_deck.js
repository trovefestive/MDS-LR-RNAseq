// =============================================================================
// 90_build_deck.js — summary deck (results/reports/MDS_LR_RNAseq_analysis_summary.pptx) from the phase figures
//   usage (from the project root, node >= 18):  NODE_PATH=tmp/deck_node/node_modules node scripts/90_build_deck.js
//   deps: pptxgenjs (npm install --prefix tmp/deck_node pptxgenjs). Optional: APPLY_THEME_JS=<path to apply_theme.js>
//   writes the theme colours into the file (otherwise scheme colours fall back to Office defaults).
// =============================================================================
const pptxgen = require("pptxgenjs");
const path = require("path");
const PROJ = path.resolve(__dirname, "..");
const R = (p) => path.join(PROJ, "results", p);
const OUT = path.join(PROJ, "results/reports/MDS_LR_RNAseq_analysis_summary.pptx");

// Palette: project genotype colours (WT blue 3B6FB6, Q157R red C8553D) on a dark navy frame
const THEME = {
  name: "LR RNA-seq U2af1", headFontFace: "Cambria", bodyFontFace: "Calibri",
  colors: { dk1: "1B2233", lt1: "FFFFFF", dk2: "4A5568", lt2: "EEF1F6",
            accent1: "C8553D", accent2: "3B6FB6", accent3: "1B7837", accent4: "E08214", accent5: "8A94A6", accent6: "D9DEE7",
            hlink: "3B6FB6", folHlink: "8A94A6" },
};
const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";            // 13.33 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "PacBio long-read RNA-seq: U2af1 Q157R vs WT";
const C = pres.SchemeColor;
const W = 13.33, H = 7.5, M = 0.6;

// === layouts ===
pres.defineSlideMaster({
  title: "DARK", background: { color: THEME.colors.dk1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: M, y: 2.3, w: W - 2 * M, h: 1.4, fontSize: 40, bold: true, color: C.background1, align: "left", margin: 0 } } },
    { placeholder: { options: { name: "body", type: "body", x: M, y: 3.8, w: W - 2 * M, h: 2.2, fontSize: 18, color: C.background2, align: "left", margin: 0 } } },
  ],
});
pres.defineSlideMaster({
  title: "CONTENT", background: { color: THEME.colors.lt1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: M, y: 0.35, w: W - 2 * M, h: 0.95, fontSize: 26, bold: true, color: C.text1, align: "left", valign: "middle", margin: 0 } } },
    { text: { text: "U2af1 Q157R vs WT  |  all effects are Q157R − WT", options: { x: M, y: H - 0.45, w: 8, h: 0.3, fontSize: 10, color: C.accent5, margin: 0 } } },
  ],
  slideNumber: { x: W - 1.2, y: H - 0.45, w: 0.6, h: 0.3, fontSize: 10, color: THEME.colors.accent5 },
});

// === helpers ===
let n = 0;
const slide = (title, section) => { const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: section }); s.addText(title, { placeholder: "title" }); n++; return s; };
const stat = (s, x, y, w, big, label, color, fs) => {
  s.addText(big, { x, y, w, h: 0.9, fontSize: fs || 40, bold: true, color: color || C.accent1, isTextBox: true, margin: 0, valign: "bottom" });
  s.addText(label, { x, y: y + 0.92, w, h: 0.7, fontSize: 12, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
};
const card = (s, x, y, w, h, head, text, fill) => {
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h, rectRadius: 0.08, fill: { color: fill || C.background2 }, line: { color: fill || C.background2 } });
  s.addText(head, { x: x + 0.2, y: y + 0.12, w: w - 0.4, h: 0.4, fontSize: 14, bold: true, color: C.text1, isTextBox: true, margin: 0 });
  s.addText(text, { x: x + 0.2, y: y + 0.52, w: w - 0.4, h: h - 0.62, fontSize: 12, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
};
const img = (s, p, x, y, w, h) => s.addImage({ path: R(p), x, y, w, h, sizing: { type: "contain", w, h } });
const note = (s, t) => s.addNotes(t);
const tbl = (s, rows, x, y, w, colW, fs) => s.addTable(rows.map((r, i) => r.map((c) => ({ text: String(c), options: { bold: i === 0, color: i === 0 ? THEME.colors.lt1 : THEME.colors.dk1, fill: { color: i === 0 ? THEME.colors.dk2 : (i % 2 ? THEME.colors.lt2 : THEME.colors.lt1) } } }))),
  { x, y, w, colW, fontSize: fs || 11, fontFace: "Calibri", border: { type: "solid", color: "FFFFFF", pt: 1 }, margin: 0.05 });

// ===================== 1. Title =====================
pres.addSection({ title: "Overview" });
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Overview" });
  s.addText("PacBio long-read RNA-seq of U2af1 Q157R mouse LK cells", { placeholder: "title" });
  s.addText([{ text: "Isoform discovery, splice-site signature and differential isoform usage, with short-read and annotation-free replication", options: { breakLine: true } },
             { text: "Analysis summary, Phases 0–8  |  October 2026", options: { fontSize: 14, color: C.accent5 } }], { placeholder: "body" });
  s.addShape(pres.ShapeType.ellipse, { x: W - 3.5, y: 1.0, w: 1.0, h: 1.0, fill: { color: C.accent2 }, line: { color: C.accent2 } });
  s.addShape(pres.ShapeType.ellipse, { x: W - 2.2, y: 1.0, w: 1.0, h: 1.0, fill: { color: C.accent1 }, line: { color: C.accent1 } });
  s.addText("WT", { x: W - 3.5, y: 2.05, w: 1.0, h: 0.3, fontSize: 11, color: C.accent5, align: "center", isTextBox: true, margin: 0 });
  s.addText("Q157R", { x: W - 2.2, y: 2.05, w: 1.0, h: 0.3, fontSize: 11, color: C.accent5, align: "center", isTextBox: true, margin: 0 });
  note(s, "Deck summarises the reanalysis following Miller et al. 2026 (bioRxiv): ESPRESSO → SQANTI3 → rMATS-long, adapted for PacBio Kinnex FLNC reads.");
}

// ===================== 2. Design =====================
{
  const s = slide("Study design: 3 WT vs 3 U2af1 Q157R, Lin⁻ Kit⁺ cells", "Overview");
  stat(s, M, 1.5, 2.8, "6", "mice, n = 3 per genotype", C.accent2);
  stat(s, M + 2.9, 1.5, 2.8, "50 M", "full-length (FLNC) reads, PacBio Kinnex / Revio");
  stat(s, M + 5.8, 1.5, 2.8, "~2 kb", "median read length, polyA required", C.accent2);
  stat(s, M + 8.7, 1.5, 3.2, "GRCm39", "GENCODE vM39; refTSS CAGE, PolyASite");
  tbl(s, [["Sample", "Genotype", "FLNC reads", "Mean length (bp)"],
          ["V335", "WT", "10.27 M", "2,050"], ["V334", "WT", "5.68 M", "2,106"], ["A310", "WT", "9.15 M", "2,132"],
          ["X504", "Q157R", "7.16 M", "2,284"], ["A258", "Q157R", "9.83 M", "2,083"], ["A309", "Q157R", "8.30 M", "1,987"]],
      M, 3.5, 6.2, [1.4, 1.4, 1.7, 1.7], 12);
  card(s, 7.3, 3.5, 5.4, 1.45, "Reference paper", "Miller et al. 2026 (bioRxiv): long-read cDNA in AML/MDS; ESPRESSO assembly, SQANTI3 curation, rMATS-long. Adapted here for HiFi FLNC reads (no pychopper; minimap2 splice:hq).");
  card(s, 7.3, 5.1, 5.4, 1.45, "Hard rule on direction", "WT is always the reference: every logFC, ΔPSI and Δproportion is Q157R − WT. Each script asserts that mutant-exclusive isoforms get a positive effect.", C.background2);
  note(s, "V334 is the shallowest sample at about half the depth of the others. No matched short reads from these mice; a separate LSK short-read cohort is used in Phase 7.");
}

// ===================== 3. Pipeline =====================
{
  const s = slide("Pipeline: eight phases, one assembler, orthogonal checks", "Overview");
  const steps = [["1", "QC + alignment", "minimap2 splice:hq; genotype pileup"], ["2", "Assembly", "ESPRESSO per chromosome; IsoQuant cross-check"],
                 ["3", "Curation", "SQANTI3 rules filter; requantify; ≥2 samples, ≥5 reads"], ["4", "ORF / NMD", "coding + NMD share, novel vs known"],
                 ["5", "Splice signature", "+1 G after AG (Q157); −3 C/T control"], ["6", "Differential isoforms", "rMATS-long, DRIMSeq, edgeR"],
                 ["7", "Short-read replication", "Salmon, STAR, rMATS-turbo on expanded GTF; LeafCutter on both datasets"], ["8", "Candidates", "evidence tiers for wet-lab follow-up"]];
  const bw = 2.85, bh = 1.75, gx = 0.25, gy = 0.3, x0 = M + 0.1, y0 = 1.5;
  steps.forEach((st, i) => {
    const x = x0 + (i % 4) * (bw + gx), y = y0 + Math.floor(i / 4) * (bh + gy);
    s.addShape(pres.ShapeType.roundRect, { x, y, w: bw, h: bh, rectRadius: 0.08, fill: { color: C.background2 }, line: { color: C.background2 } });
    s.addShape(pres.ShapeType.ellipse, { x: x + 0.2, y: y + 0.2, w: 0.55, h: 0.55, fill: { color: i === 6 ? C.accent2 : C.accent1 }, line: { color: i === 6 ? C.accent2 : C.accent1 } });
    s.addText(st[0], { x: x + 0.2, y: y + 0.2, w: 0.55, h: 0.55, fontSize: 16, bold: true, color: C.background1, align: "center", valign: "middle", isTextBox: true, margin: 0 });
    s.addText(st[1], { x: x + 0.9, y: y + 0.22, w: bw - 1.05, h: 0.5, fontSize: 14, bold: true, color: C.text1, valign: "middle", isTextBox: true, margin: 0 });
    s.addText(st[2], { x: x + 0.2, y: y + 0.9, w: bw - 0.4, h: 0.8, fontSize: 12, color: C.text2, valign: "top", isTextBox: true, margin: 0 });
  });
  s.addText("Blue = separate short-read dataset (different mice, LSK cells). All other phases use the long reads.", { x: M, y: 5.85, w: 12, h: 0.4, fontSize: 12, italic: true, color: C.accent5, isTextBox: true, margin: 0 });
}

// ===================== 4. Phase 1 =====================
pres.addSection({ title: "Long reads" });
{
  const s = slide("Phase 1: clean alignment; every genotype confirmed from the reads", "Long reads");
  stat(s, M, 1.5, 2.6, "≥99.9%", "of FLNC reads mapped, every sample", C.accent2);
  stat(s, M, 3.2, 2.6, "0.21%", "per-base mismatch rate (HiFi)", C.accent2);
  stat(s, M, 4.9, 2.6, "40–42%", "Q157R allele in A258, A309, X504; 0–0.2% in WT", undefined, 36);
  img(s, "01_align/phase1_u2af1_Q157R_genotype.png", 3.6, 1.4, 6.2, 4.4);
  card(s, 10.1, 1.5, 2.7, 4.2, "Why this matters", "Genotype is read directly at chr17:31,867,169 (T>C, codon 157). Heterozygous allele fractions in the three mutants and none in WT, so no sample swap. The same check is repeated on the short-read cohort in Phase 7.");
  note(s, "minimap2 2.31 -ax splice:hq -uf with the GENCODE junction BED; samtools pileup at the Q157R base.");
}

// ===================== 5. Phase 2 =====================
{
  const s = slide("Phase 2: 253 k ESPRESSO isoforms; the U2af1-213 isoform is a cis effect", "Long reads");
  stat(s, M, 1.45, 2.4, "253 k", "isoforms before curation (127 k novel)");
  stat(s, M + 2.5, 1.45, 2.4, "81–85%", "of reads assigned per sample", C.accent2);
  stat(s, M + 5.0, 1.45, 2.4, "100%", "of U2af1-213 junction reads carry the C allele");
  img(s, "06_diff/phase6_dotplot_U2af1.png", M, 3.2, 7.4, 3.6);
  card(s, 8.3, 1.5, 4.4, 2.4, "U2af1-213 (ENSMUST00000468653)", "Present only in Q157R (CPM 4.6–11.4 vs 0–0.9). Every read across its new junction carries the mutant base: the Q157R mutation creates a 5′ splice site in U2af1 itself. It is UP in the mutant and is used as the positive control for direction in every later step.");
  card(s, 8.3, 4.1, 4.4, 2.7, "Cross-check and composition", "IsoQuant run on the same BAMs as an orthogonal assembler. A310 (WT) carries about 2× Mpo/Elane, i.e. more promyelocyte-like cells in that sort; granule-gene differences are read with caution.", C.background2);
  note(s, "ESPRESSO 1.6.0 S/C/Q per chromosome (-T 8 with retries under the cluster's strict overcommit). The old LRP2 run reported this isoform as 'down' because its A_vs_B naming uses A as reference; see the Phase 6 slide.");
}

// ===================== 6. Phase 3 =====================
{
  const s = slide("Phase 3: SQANTI3 curation keeps 171,224 isoforms, 44,907 new", "Long reads");
  stat(s, M, 1.45, 2.3, "126 k", "GENCODE isoforms kept (reference treated as truth)", C.accent2);
  stat(s, M + 2.4, 1.45, 2.3, "44.9 k", "new isoforms: ≥2 samples, ≥5 reads");
  stat(s, M + 4.8, 1.45, 2.3, "4.2–4.7%", "of expression from new isoforms, every sample", C.accent2, 30);
  img(s, "03_sqanti/phase3_novel_categories.png", M, 3.2, 6.6, 3.6);
  card(s, 7.5, 1.5, 5.2, 2.5, "Filter choice: rules, not ML", "The paper used SQANTI3's ML filter. On our data its default kept only 337 novel isoforms with a coding/NMD bias, and without ORF features it dropped CAGE- and polyA-supported short isoforms. The rules filter (CAGE, polyA, RT-switching, canonical junctions) is primary; ML is reported for comparison.");
  card(s, 7.5, 4.2, 5.2, 2.6, "Counts and the FSM-read filter", "Primary counts use all assigned reads (R2). The paper's FSM-only read filter removes 15–20% of reads here (vs ~47% for ONT), kept as a sensitivity check. IsoQuant rebuilt the identical intron chain for 26% of GENCODE and 6% of novel isoforms.", C.background2);
  note(s, "Final set: 32,736 novel-in-catalog, 11,076 novel-not-in-catalog, 490 fusion. 43,347 of the novel isoforms are supported in ≥2 replicates of the same genotype. 2,926 novel isoforms exceed 10% of their gene's expression.");
}

// ===================== 7. Phase 4 =====================
{
  const s = slide("Phase 4: new isoforms are NMD-prone; NMD load does not shift", "Long reads");
  img(s, "04_nmd/phase4_nmd_known_vs_novel.png", M, 1.4, 5.6, 4.6);
  stat(s, 6.7, 1.5, 2.9, "36.0%", "of new coding isoforms predicted NMD vs 14.0% of GENCODE (OR 3.4)");
  stat(s, 9.8, 1.5, 2.9, "1.88 vs 1.91%", "of expression from NMD isoforms, Q157R vs WT (n.s.)", C.accent2, 28);
  card(s, 6.7, 3.6, 6.0, 2.9, "What this says", "Novel-in-catalog isoforms are the most NMD-prone (42%). At the per-gene level, 0 of 5,907 genes show a significant change in NMD share (p-values uniform). Q157R does not raise global NMD burden in LK cells; the mutation's effect is on which isoforms are made, not on how many are degraded. Wilcoxon was avoided: with 3 vs 3 its minimum two-sided p is 0.10.");
}

// ===================== 8. Phase 5 =====================
{
  const s = slide("Phase 5: Q157R-gained 3′ splice sites carry the +1 G signature", "Long reads");
  img(s, "05_splice/phase5_3ss_logos.png", M, 1.35, 6.4, 5.5);
  stat(s, 7.4, 1.45, 2.6, "82.5%", "+1 G at Q157R-associated new 3′SS (n = 127)");
  stat(s, 10.2, 1.45, 2.6, "44.9%", "+1 G at WT-associated new 3′SS (n = 108)", C.accent2);
  card(s, 7.4, 3.3, 5.3, 1.6, "Fisher's exact test", "+1 G vs A: OR 5.70, p = 2.2 × 10⁻⁶. Annotated 3′SS background: 64% G. Q157R sites are enriched above background; WT-associated sites are depleted.");
  card(s, 7.4, 5.05, 5.3, 1.8, "Negative control holds", "The S34F-type −3 C vs T signature shows no difference (OR 1.04, p = 1). GO over-representation of genes with Q157R-enriched isoforms: nothing at FDR < 0.05 (RNA splicing FDR 0.14).", C.background2);
  note(s, "Q157 contacts the +1 position after the 3′SS AG; the paper tested Q157P. Associated isoforms are mostly detection-pattern specific (854 Q157R / 801 WT novel), 23 are usage-enriched.");
}

// ===================== 9. Phase 6 =====================
pres.addSection({ title: "Differential isoforms" });
{
  const s = slide("Phase 6: hundreds of usage switches, mostly without gene-level change", "Differential isoforms");
  stat(s, M, 1.45, 2.4, "866", "isoforms / 576 genes, rMATS-long (paper's tool)");
  stat(s, M + 2.5, 1.45, 2.4, "130", "isoforms, DRIMSeq (Dirichlet-multinomial)", C.accent2);
  stat(s, M + 5.0, 1.45, 2.4, "122", "found by both, 100% same sign, r = 0.90");
  stat(s, M + 7.5, 1.45, 2.4, "97", "genes DE at gene level (89 up; ribosomal/OXPHOS, flagged)", C.accent2);
  img(s, "06_diff/phase6_rmatslong_Cd34_structure.png", M, 3.35, 7.6, 2.4);
  s.addText("Top switch: Cd34. Cd34-201 (blue) +0.22 and Cd34-202 (red) −0.23 in Q157R; the two differ only by an alternative 3′SS 156 nt apart in the last exon. Total Cd34 unchanged.",
    { x: M, y: 5.85, w: 7.6, h: 0.9, fontSize: 12, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
  card(s, 8.5, 3.35, 4.2, 3.4, "Event classes (top 100 genes)", "Complex 40  |  exon skipping 35  |  alternative 3′SS 13  |  intron retention 10  |  alternative first exon 9  |  alternative ends 4  |  ALE 2  |  A5SS 1. Exon skipping and A3SS dominate the simple events, as expected for a 3′SS factor. 577 of 606 rMATS-long switch genes have no gene-level change. No Ybx1 exon-skipping shift (the paper's event was SRSF2-driven).");
  note(s, "rMATS-long run with group 1 = Q157R, group 2 = WT, --num-threads 1 (thread-join race). The ribosomal/OXPHOS DGE signal is treated as an observation to confirm; it was not replicated in Phase 7.");
}

// ===================== 10. LRP2 sign convention =====================
{
  const s = slide("Earlier LRP2 run agrees once its direction label is read correctly", "Differential isoforms");
  card(s, M, 1.45, 6.0, 2.6, "The puzzle that started this reanalysis", "The LRP2 output reported the mutant-exclusive U2af1 isoform as DOWN (logFC −8.80) in 'Q157R_vs_WT'. LRP2 names contrasts A_vs_B with A as the reference, so its values are WT relative to Q157R. Negating them gives Q157R − WT (confirmed with the LRP2 developer). Across 10,176 genes: slope −1.03, r = −0.999 against a recomputed logFC.");
  tbl(s, [["Level", "Shared features", "r (after negation)", "Significant here", "Also in LRP2", "Same sign"],
          ["Gene (DGE)", "9,900", "0.95", "94", "92", "100%"], ["Isoform (DTE)", "19,700", "0.95", "60", "60", "100%"], ["Isoform usage (DTU)", "13,061", "0.76", "102", "77", "100%"]],
      M, 4.3, 6.0, [1.6, 1.1, 1.1, 0.8, 0.7, 0.7], 11);
  stat(s, 7.2, 1.5, 2.6, "+8.86", "U2af1-213 logFC here (Q157R − WT)");
  stat(s, 10.0, 1.5, 2.6, "+8.80", "the same isoform in LRP2 after negation");
  card(s, 7.2, 3.4, 5.5, 2.9, "Conclusion", "The two pipelines agree closely. The sign disagreement was entirely the naming convention, not an error in either analysis. This is why every script in the new pipeline asserts direction on mutant-exclusive isoforms and every table states its sign convention.", C.background2);
}

// ===================== 11. Phase 7 design =====================
pres.addSection({ title: "Short-read replication" });
{
  const s = slide("Phase 7: independent short-read cohort (LSK) on the expanded GTF", "Short-read replication");
  stat(s, M, 1.45, 2.4, "3 vs 3", "WT vs Q157R knock-in mice, Illumina 2×150, stranded", C.accent2);
  stat(s, M + 2.5, 1.45, 2.4, "56–94 M", "read pairs per sample");
  stat(s, M + 5.0, 1.45, 2.4, "67–73%", "uniquely mapped to GRCm39 (STAR)", C.accent2);
  stat(s, M + 7.5, 1.45, 2.4, "39–43%", "Q157R allele in KI1–3; 0–0.2% in WT");
  tbl(s, [["Sample", "Expected", "Depth at codon 157", "Alt C fraction", "Call"],
          ["SR_WT1", "WT", "1,105", "0.001", "WT"], ["SR_WT2", "WT", "324", "0.000", "WT"], ["SR_WT3", "WT", "801", "0.002", "WT"],
          ["SR_KI1", "Q157R", "830", "0.394", "Q157R"], ["SR_KI2", "Q157R", "668", "0.431", "Q157R"], ["SR_KI3", "Q157R", "1,475", "0.411", "Q157R"]],
      M, 3.4, 6.4, [1.3, 1.1, 1.6, 1.3, 1.1], 11);
  card(s, 7.5, 3.4, 5.2, 3.3, "What was run", "Expanded GTF = the 171,224 final isoforms. Salmon (decoy-aware, stranded) for isoform shares; STAR for BAMs; rMATS-turbo (b1 = Q157R, b2 = WT, --novelSS) for events. Different animals and a more primitive population (LSK vs LK): this tests replication of mutation effects, not matched validation.");
  note(s, "rMATS-turbo segfaulted on the all-BAM run; it was split into per-BAM prep + post. Salmon's library auto-detect called 2 of 6 samples unstranded although 86–92% of fragments were ISR; library type was fixed to ISR.");
}

// ===================== 12. Phase 7 replication =====================
{
  const s = slide("Phase 7: splicing replicates across cohorts; gene expression does not", "Short-read replication");
  s.addText("Isoform usage (122 high-confidence long-read switches)", { x: M, y: 1.35, w: 6, h: 0.4, fontSize: 14, bold: true, color: C.text1, isTextBox: true, margin: 0 });
  stat(s, M, 1.8, 2.9, "83 / 121", "same direction in short reads (69%); binomial p = 5 × 10⁻⁵", C.accent3);
  stat(s, M + 3.1, 1.8, 2.9, "r = 0.44", "Δ proportion, long vs short reads", C.accent3);
  s.addText("Gene-level expression (10,854 shared genes)", { x: M, y: 3.75, w: 6, h: 0.4, fontSize: 14, bold: true, color: C.text1, isTextBox: true, margin: 0 });
  stat(s, M, 4.2, 2.9, "r = 0.18", "logFC, long vs short reads");
  stat(s, M + 3.1, 4.2, 2.9, "0", "DE genes in short reads; 9 of 34 ribosomal genes same direction");
  s.addText("Detection of long-read isoforms", { x: 7.3, y: 1.35, w: 5.5, h: 0.4, fontSize: 14, bold: true, color: C.text1, isTextBox: true, margin: 0 });
  stat(s, 7.3, 1.8, 2.6, "61%", "of 44,856 new isoforms detected (≥2 samples, ≥5 reads)", C.accent2);
  stat(s, 10.1, 1.8, 2.6, "64%", "of GENCODE isoforms detected", C.accent2);
  card(s, 7.3, 3.75, 5.4, 2.9, "Reading", "Usage switches behave like a mutation-intrinsic effect: they hold in different mice and in a different progenitor population. The ribosomal/OXPHOS increase seen in long reads is not supported. Whether that is the cell type (LSK vs LK) or a library/length effect in the long reads cannot be separated here.", C.background2);
}

// ===================== 13. Cd34 + +1G in short reads =====================
{
  const s = slide("Phase 7: Cd34 replicates by a third method; +1 G preference is weaker", "Short-read replication");
  img(s, "06_diff/phase6_dotplot_Cd34.png", M, 1.35, 7.0, 3.3);
  s.addText("Long reads: Cd34 isoform CPM per replicate (Phase 6).", { x: M, y: 4.65, w: 7, h: 0.3, fontSize: 11, italic: true, color: C.accent5, isTextBox: true, margin: 0 });
  card(s, 8.0, 1.45, 4.7, 2.3, "rMATS-turbo, same A3SS", "PSI of the upstream acceptor 0.07–0.12 in Q157R vs 0.23–0.28 in WT, ΔPSI −0.16, FDR ≈ 0, no overlap between replicates. Salmon shares agree (Cd34-201 +0.08, Cd34-202 −0.11). Same shift to the downstream acceptor, in a second cohort and cell type.");
  card(s, 8.0, 3.9, 4.7, 2.8, "+1 G at Q157R-favoured acceptors", "Among 1,624 significant A3SS events: favoured acceptor 61.3% G vs disfavoured 55.9% (OR 1.25, p = 0.013). Same direction as Phase 5 but much weaker than 82.5% vs 44.9% (OR 5.7). The tests differ: within-event acceptor pairs here, Q157R- vs WT-associated novel junctions there. 99% of acceptors have AG at −2/−1 (coordinate check).", C.background2);
  s.addText("rMATS-turbo: 89,547 events tested; 8,886 at FDR < 0.05 and |ΔPSI| ≥ 0.1 (SE 2,852, RI 2,671, A3SS 1,624, A5SS 1,144, MXE 595). Liberal at n = 3; read as an upper bound.",
    { x: M, y: 5.2, w: 7.0, h: 1.4, fontSize: 12, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
}

// ===================== 14. LeafCutter =====================
{
  const s = slide("Phase 7b: LeafCutter, annotation-free, reproduces the signature in both datasets", "Short-read replication");
  img(s, "07_leafcutter/compare/phase7b_short_vs_long_dpsi.png", M, 1.35, 5.4, 5.0);
  stat(s, 6.3, 1.45, 3.1, "OR 9.7", "+1 G, Q157R-favoured vs disfavoured acceptor, long reads (70.7% vs 19.8%; p = 5 × 10⁻³⁶)");
  stat(s, 9.6, 1.45, 3.1, "OR 3.5", "same test, short reads (69.1% vs 38.9%; p = 4 × 10⁻⁵)");
  tbl(s, [["Dataset", "Clusters tested", "Significant", "−3 C/T control"],
          ["Long reads (LK)", "13,543", "295", "OR 0.95, p = 0.75"], ["Short reads (LSK)", "7,270", "168", "OR 0.80, p = 0.42"]],
      6.3, 3.35, 6.4, [1.9, 1.5, 1.3, 1.7], 11);
  card(s, 6.3, 4.6, 6.4, 2.1, "No assembler, no annotation", "Introns are clustered straight from split reads (regtools → LeafCutter 0.2.9), independent of ESPRESSO, SQANTI3 and GENCODE. The U2af1-213 cis junction is up in both datasets (ΔPSI +0.12 / +0.10). Cd34 shows the same last-exon acceptor shift a fourth time (ΔPSI −0.24 long, −0.17 short). Among introns significant in either dataset, long and short reads agree in direction 75% of the time (r = 0.52).", C.background2);
  note(s, "Acceptor pairs = two acceptors sharing a donor inside a significant cluster; favoured = larger ΔPSI. n = 3 per group is below LeafCutter's calibrated range, so adjusted p-values are approximate. regtools RF/FR modes assume paired reads and mislabelled every long-read strand; fixed with an XS tag derived from the read flag (FLNC reads are transcript-oriented).");
}

// ===================== 15. Phase 8 tiers =====================
pres.addSection({ title: "Candidates" });
{
  const s = slide("Phase 8: 918 candidate isoforms scored on six kinds of evidence", "Candidates");
  img(s, "08_summary/phase8_evidence_top_genes.png", M, 1.3, 5.3, 5.7);
  tbl(s, [["Tier", "Rule", "Genes"],
          ["A", "both methods + LR replicates separate + SR same direction + IsoQuant not contradicting", "36"],
          ["B", "both methods + (SR same direction or IsoQuant agrees)", "35"],
          ["C", "everything else (mostly rMATS-long only)", "540"]],
      6.3, 1.4, 6.4, [0.7, 4.9, 0.8], 11);
  card(s, 6.3, 3.7, 6.4, 3.0, "What the tiers mean", "Evidence per isoform: rMATS-long + DRIMSeq agreement; long-read replicates fully separated; ≥10 reads; IsoQuant built the same intron chain with the same Δ sign; short-read share same direction; short-read replicates separated. None of the tier A/B isoforms is predicted NMD; 31 of 36 tier A genes have no gene-level change. Across all 907 short-read-testable candidates agreement is weak (r = 0.18), so rMATS-long-only calls are treated as weak.", C.background2);
}

// ===================== 16. Tier A =====================
{
  const s = slide("Tier A switches hold in every replicate of both datasets", "Candidates");
  img(s, "08_summary/phase8_tierA_proportions.png", M, 1.3, 12.1, 2.6);
  s.addText("Top: long reads (LK). Bottom: short reads (LSK). Isoform proportion within gene, every replicate. Absolute levels differ between methods when isoforms share most of their sequence; direction is the comparable quantity.",
    { x: M, y: 3.95, w: 12.1, h: 0.5, fontSize: 11, italic: true, color: C.accent5, isTextBox: true, margin: 0 });
  const cards = [["Cd34  |  A3SS, last exon", "Significant in rMATS-long, DRIMSeq, short-read rMATS-turbo and LeafCutter; two cohorts, two cell types. RT-PCR across the two acceptors gives products 156 bp apart."],
                 ["Score-6 switches", "Cdc14a, Mtm1, Erbin, Crlf3: every scored evidence type met. Crlf3 is exon skipping; Cdc14a also changes at gene level."],
                 ["A3SS trio", "Atn1, Tcf19, Thada: the event class expected for a U2AF1 mutation. Check whether the Q157R-favoured acceptor carries +1 G."],
                 ["OXPHOS / mito-ribosome", "Atp5if1, Atp5mg, Cox16, Ndufv3, Ndufa7, Ndufa13 (novel), Ndufb8, Mrpl23, Mrpl33. Usage shifts replicate; their gene-level increases did not."]];
  cards.forEach((c, i) => card(s, M + i * 3.07, 4.6, 2.9, 2.2, c[0], c[1], i === 0 ? undefined : C.background2));
}

// ===================== 17. Caveats / next =====================
{
  const s = pres.addSlide({ masterName: "DARK", sectionTitle: "Candidates" });
  s.addText("Where this leaves us", { x: M, y: 0.6, w: W - 2 * M, h: 1.0, fontSize: 40, bold: true, color: C.background1, isTextBox: true, margin: 0, valign: "middle" });
  const col = (x, head, items) => {
    s.addText(head, { x, y: 1.8, w: 3.9, h: 0.5, fontSize: 18, bold: true, color: C.accent1, isTextBox: true, margin: 0 });
    s.addText(items.map((t, i) => ({ text: t, options: { bullet: { indent: 14 }, breakLine: i < items.length - 1, paraSpaceAfter: 8 } })),
      { x, y: 2.4, w: 3.9, h: 4.3, fontSize: 13, color: C.background2, valign: "top", isTextBox: true, margin: 0 });
  };
  col(M, "Established", ["Q157R creates a cis splice site in U2af1 itself (U2af1-213)", "Gained 3′ splice sites prefer +1 G: OR 5.7 (ESPRESSO), 9.7 / 3.5 (LeafCutter, long / short); −3 control null", "Isoform-usage switches replicate in an independent LSK cohort", "Cd34 last-exon A3SS: four methods, two cohorts", "LRP2 and new pipeline agree once signs are read correctly"]);
  col(M + 4.2, "Caveats", ["n = 3 vs 3 in both datasets", "Short reads are different mice and a different population (LSK vs LK)", "Gene-level ribosomal/OXPHOS increase not replicated", "A310 (WT) has a higher promyelocyte signature", "rMATS-turbo, rMATS-long and LeafCutter are liberal at this n"]);
  col(M + 8.4, "Next", ["RT-PCR / RT-qPCR on Cd34 and the score-6 switches", "+1 G check on Atn1, Tcf19, Thada acceptors", "Optional: custom protein FASTA from novel ORFs if proteomics becomes available", "Code and small results are on GitHub; scratch outputs need syncing to project storage"]);
}

(async () => {
  await pres.writeFile({ fileName: OUT });
  if (process.env.APPLY_THEME_JS) { const { applyTheme } = require(process.env.APPLY_THEME_JS); await applyTheme(OUT, THEME); }
  console.log("wrote", OUT, "slides:", n + 2);
})();
