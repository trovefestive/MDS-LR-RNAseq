#!/usr/bin/env python3
# =============================================================================
# 11_phase2_summary.py — Phase 2: summarize the merged ESPRESSO matrix (pre-filter)
#   - per-sample assigned reads vs FLNC, known vs novel isoforms, novel read fraction
#   - preview of the paper's novel-isoform rule (>=2 samples, >=5 total reads)
#   - top isoforms over-represented in A310 (read-length spike seen in Phase 1)
#   - U2af1 isoform counts, and the allele at chr17:31867169 (Q157R) in reads using the
#     U2af1-213 junction (intron chr17:31867070-31867168) vs the canonical junction
#     (intron chr17:31867070-31867156)
# Run from the project root after `source config.sh` (needs PROJ; samtools on PATH). Small: login node OK.
# Outputs: results/02_espresso/phase2_*.tsv ; gzipped matrix copy for the repo.
# =============================================================================
import csv, gzip, os, re, shutil, subprocess, sys

proj = os.environ.get("PROJ") or sys.exit("PROJ not set: source config.sh first")
esp = os.path.join(proj, "results/02_espresso")
abund = os.path.join(esp, "merged/merged_N2_R0_abundance.esp")

# === load matrix ===
print("[phase2] reading", os.path.relpath(abund, proj))
with open(abund) as fh:
    r = csv.reader(fh, delimiter="\t"); h = next(r); S = h[3:]
    rows = [(x[0], x[1], x[2], [float(v) for v in x[3:]]) for x in r]
flnc = {}
with open(os.path.join(proj, "samples.tsv")) as fh:
    for d in csv.DictReader(fh, delimiter="\t"):
        st = os.path.join(proj, "results/01_align/qc", d["sample_id"] + ".seqkit_stats.tsv")
        with open(st) as s2: flnc[d["sample_id"]] = int(list(csv.DictReader(s2, delimiter="\t"))[0]["num_seqs"])
tot = [sum(x[3][i] for x in rows) for i in range(len(S))]
nov = [x for x in rows if x[0].startswith("ESPRESSO")]
known = [x for x in rows if not x[0].startswith("ESPRESSO")]
nov_pass = [x for x in nov if sum(c > 0 for c in x[3]) >= 2 and sum(x[3]) >= 5]

# === per-sample summary ===
with open(os.path.join(esp, "phase2_sample_summary.tsv"), "w") as out:
    out.write("sample\tflnc_reads\tassigned_reads\tpct_assigned\tpct_reads_on_novel\n")
    for i, s in enumerate(S):
        nv = sum(x[3][i] for x in nov)
        out.write(f"{s}\t{flnc[s]}\t{tot[i]:.0f}\t{100*tot[i]/flnc[s]:.1f}\t{100*nv/tot[i]:.2f}\n")
with open(os.path.join(esp, "phase2_isoform_summary.tsv"), "w") as out:
    out.write("metric\tvalue\n")
    for k, v in [("isoforms_total", len(rows)), ("isoforms_known_gencode", len(known)),
                 ("isoforms_novel_espresso", len(nov)), ("novel_ge2samples_ge5reads_prefilter", len(nov_pass))]:
        out.write(f"{k}\t{v}\n")

# === A310 over-representation (CPM) ===
cpm = lambda x, i: x[3][i] / tot[i] * 1e6
ia = S.index("A310"); oth = [i for i in range(len(S)) if i != ia]
top = sorted(rows, key=lambda x: cpm(x, ia) - sum(cpm(x, i) for i in oth) / len(oth), reverse=True)[:10]
with open(os.path.join(esp, "phase2_A310_overrepresented.tsv"), "w") as out:
    out.write("transcript_ID\ttranscript_name\tgene_ID\tcpm_A310\tmean_cpm_others\n")
    for x in top:
        out.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{cpm(x, ia):.1f}\t{sum(cpm(x, i) for i in oth)/len(oth):.1f}\n")

# === U2af1 isoforms ===
with open(os.path.join(esp, "phase2_u2af1_isoforms.tsv"), "w") as out:
    out.write("transcript_ID\ttranscript_name\t" + "\t".join(S) + "\n")
    for x in rows:
        if "ENSMUSG00000061613" in x[2] and sum(x[3]) >= 5:
            out.write(f"{x[0]}\t{x[1]}\t" + "\t".join(f"{c:.1f}" for c in x[3]) + "\n")

# === allele at Q157R position in reads using each junction ===
POS, NOVEL, CANON = 31867169, (31867070, 31867168), (31867070, 31867156)
with open(os.path.join(esp, "phase2_u2af1_junction_alleles.tsv"), "w") as out:
    out.write("sample\tjunction\tbase_T_ref\tbase_C_Q157R\tother\n")
    for s in S:
        bam = os.path.join(proj, "results/01_align/bam", s + ".sorted.bam")
        sam = subprocess.run(["samtools", "view", "-q", "20", bam, "chr17:31867100-31867200"],
                             stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout
        cnt = {"U2af1-213_novel_donor": {}, "canonical": {}}
        for l in sam.splitlines():
            f = l.split("\t"); r = int(f[3]); q = 0; base = None; juncs = []
            for n, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
                n = int(n)
                if op in "M=X":
                    if r <= POS < r + n: base = f[9][q + POS - r]
                    r += n; q += n
                elif op in "IS": q += n
                elif op == "D": r += n
                elif op == "N": juncs.append((r, r + n - 1)); r += n
            k = "U2af1-213_novel_donor" if NOVEL in juncs else ("canonical" if CANON in juncs else None)
            if k and base: cnt[k][base] = cnt[k].get(base, 0) + 1
        for k, c in cnt.items():
            out.write(f"{s}\t{k}\t{c.get('T',0)}\t{c.get('C',0)}\t{sum(v for b,v in c.items() if b not in 'TC')}\n")

# === gzipped matrix copy (small enough for the repo) ===
gz = os.path.join(esp, "merged_N2_R0_abundance.esp.gz")
with open(abund, "rb") as fi, gzip.open(gz, "wb") as fo: shutil.copyfileobj(fi, fo)

for f in ["phase2_sample_summary.tsv", "phase2_isoform_summary.tsv", "phase2_u2af1_junction_alleles.tsv"]:
    print(f"== {f}"); print(open(os.path.join(esp, f)).read())
print("[phase2] done")
