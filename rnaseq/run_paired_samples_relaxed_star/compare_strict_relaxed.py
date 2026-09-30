#!/usr/bin/env python3
"""Compare the strict (default STAR) and relaxed STAR Diatraea analyses.

Usage: python3 compare_strict_relaxed.py <strict_star_salmon> <strict_downstream_dir>
                                         <relaxed_star_salmon> <relaxed_downstream_dir> <out_prefix>
Writes <out_prefix>.tsv (all numbers) and <out_prefix>.md (readable tables).
"""
import csv, math, sys
from pathlib import Path

strict_ss, strict_dir, relax_ss, relax_dir, out = map(Path, sys.argv[1:6])
SAMPLES = ["control_rep1", "control_rep2", "control_rep3", "infected_rep1", "infected_rep2", "infected_rep3"]


def star_log(ss, s):
    vals = {}
    for line in open(ss / "log" / f"{s}.Log.final.out"):
        if "|" in line:
            k, v = line.split("|"); vals[k.strip()] = v.strip().rstrip("%")
    u, m = float(vals["Uniquely mapped reads %"]), float(vals["% of reads mapped to multiple loci"])
    return u, m, u + m, int(vals["Number of input reads"])


def salmon_reads(ss, s):
    rows = list(csv.reader(open(ss / s / "quant.sf"), delimiter="\t"))[1:]
    return sum(float(r[4]) for r in rows)


def read_table(fn, key=0):
    rows = list(csv.reader(open(fn)))
    header = rows[0]
    # DESeq2 tables have one column less in the header (row names)
    off = 1 if len(rows) > 1 and len(rows[1]) == len(header) + 1 else 0
    return {r[key]: dict(zip(header, r[off:])) for r in rows[1:]}


def overlap(a, b):
    a, b = set(a), set(b)
    return len(a), len(b), len(a & b), len(a - b), len(b - a)


tsv, md = [], []
tsv.append(["section", "item", "strict", "relaxed", "shared", "strict_only", "relaxed_only"])

# 1. mapping
md += ["## Mapping (per sample)", "",
       "| Sample | Input pairs | STAR total mapped % strict | relaxed | Unique % strict | relaxed | Salmon reads strict (M) | relaxed (M) | Salmon gain |",
       "|---|---|---|---|---|---|---|---|---|"]
for s in SAMPLES:
    us, ms, ts, n = star_log(strict_ss, s)
    ur, mr, tr, _ = star_log(relax_ss, s)
    qs, qr = salmon_reads(strict_ss, s), salmon_reads(relax_ss, s)
    tsv += [["STAR_total_mapped_pct", s, f"{ts:.2f}", f"{tr:.2f}", "", "", ""],
            ["STAR_unique_pct", s, f"{us:.2f}", f"{ur:.2f}", "", "", ""],
            ["Salmon_reads", s, f"{qs:.0f}", f"{qr:.0f}", "", "", ""]]
    md.append(f"| {s} | {n/1e6:.1f} M | {ts:.1f} | **{tr:.1f}** | {us:.1f} | {ur:.1f} | {qs/1e6:.2f} | {qr/1e6:.2f} | {qr/qs:.2f}x |")

# 2. DEGs
md += ["", "## Differentially expressed genes (RUVs k=1, |log2FC| > 1, padj < 0.05)", "",
       "| Set | Strict | Relaxed | Shared | Strict only | Relaxed only |", "|---|---|---|---|---|---|"]
for d in ["up", "down"]:
    a = read_table(strict_dir / f"{d}_regulated.csv"); b = read_table(relax_dir / f"{d}_regulated.csv")
    o = overlap(a, b)
    tsv.append(["DEG", d, *map(str, o)])
    md.append(f"| {d} ({'higher in control' if d == 'up' else 'higher in infected'}) | {o[0]} | {o[1]} | {o[2]} | {o[3]} | {o[4]} |")

# log2FC agreement over genes tested in both
A = read_table(strict_dir / "dea_all_genes.csv"); B = read_table(relax_dir / "dea_all_genes.csv")
shared = [g for g in A if g in B]
x = [float(A[g]["log2FoldChange"]) for g in shared]; y = [float(B[g]["log2FoldChange"]) for g in shared]
mx, my = sum(x) / len(x), sum(y) / len(y)
r = sum((i - mx) * (j - my) for i, j in zip(x, y)) / math.sqrt(sum((i - mx) ** 2 for i in x) * sum((j - my) ** 2 for j in y))
tsv.append(["genes_tested", "all", str(len(A)), str(len(B)), str(len(shared)), str(len(set(A) - set(B))), str(len(set(B) - set(A)))])
tsv.append(["log2FC_pearson_r", "shared_tested_genes", "", "", f"{r:.3f}", "", ""])
md += ["", f"Genes tested: strict {len(A)}, relaxed {len(B)}, shared {len(shared)}. "
           f"Pearson r of log2FC over shared genes: **{r:.3f}**."]

# 3. GO and KEGG
md += ["", "## Enrichment", "", "| Analysis | Strict | Relaxed | Shared | Strict only | Relaxed only |", "|---|---|---|---|---|---|"]
for kind, fn, key in [("GO", "GO_{}.csv", "GO.ID"), ("KEGG", "kegg_{}.csv", "ID")]:
    for d in ["up", "down"]:
        def terms(dirp):
            p = dirp / fn.format(d)
            if not p.exists() or p.stat().st_size == 0:
                return set()
            rows = list(csv.DictReader(open(p)))
            if kind == "KEGG":
                rows = [r for r in rows if float(r["p.adjust"]) < 0.05]
            return {r[key] for r in rows}
        o = overlap(terms(strict_dir), terms(relax_dir))
        tsv.append([kind, d, *map(str, o)])
        md.append(f"| {kind} {d} | {o[0]} | {o[1]} | {o[2]} | {o[3]} | {o[4]} |")

with open(f"{out}.tsv", "w") as f:
    csv.writer(f, delimiter="\t").writerows(tsv)
with open(f"{out}.md", "w") as f:
    f.write("# Strict vs relaxed STAR: Diatraea control vs infected\n\n"
            "Strict = default nf-core STAR filters; relaxed = `--outFilterScoreMinOverLread 0.3 "
            "--outFilterMatchNminOverLread 0.3`. Both analysed with the same scripts "
            "(RUVs k=1 DESeq2, topGO with numeric p-values, current KEGG).\n\n" + "\n".join(md) + "\n")
print("\n".join(md))
