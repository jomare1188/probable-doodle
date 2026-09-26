#!/usr/bin/env python3
"""Assign each Trinity contig of the unclassified reads to a category and
weight categories by salmon read counts. Run from maptest/uid/."""
import re, collections

MITO = re.compile(r"Cytochrome c oxidase subunit [123]|Cytochrome b(?!-)\b|NADH-ubiquinone oxidoreductase chain|"
                  r"ATP synthase subunit a\b|ATP synthase protein 8")

def fasta(fn):
    seqs, name = {}, None
    for l in open(fn):
        if l[0] == ">":
            name = l[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(l.strip())
    return {k: "".join(v) for k, v in seqs.items()}

contigs = fasta("contigs.fa")
reads = {}
for l in list(open("salmon_quant/quant.sf"))[1:]:
    f = l.split("\t"); reads[f[0]] = float(f[4])

genome = {}
for l in open("contigs_vs_genome.tsv"):
    q, s, pid, ln, ql, ev, bs = l.rstrip().split("\t")
    genome.setdefault(q, (float(pid), int(ln) / int(ql)))
sprot = {}
for l in open("contigs_vs_sprot.tsv"):
    f = l.rstrip().split("\t")
    sprot.setdefault(f[0], (float(f[1]), f[5]))
orf = {}
for l in open("contigs.fa.transdecoder.pep"):
    if l[0] == ">":
        m = re.search(r"type:(\S+)", l)
        orf.setdefault(l[1:].split(".p")[0], m.group(1) if m else "orf")
pfam = collections.defaultdict(set)
for l in open("contigs_pfam.domtbl"):
    if l[0] != "#":
        f = l.split(); pfam[f[0].split(".p")[0]].add(f[3])

rows = []
for c, seq in contigs.items():
    prot = sprot[c][1] if c in sprot else ""
    if MITO.search(prot):
        cat = "mitochondrial (absent from assembly)"
    elif c in genome:
        cat = "Diatraea genome, divergent"
    elif c in sprot or c in pfam:
        cat = "coding, no genome hit"
    elif c in orf:
        cat = "ORF only, no homology"
    else:
        cat = "no ORF, no homology"
    gp, gc = genome.get(c, ("", ""))
    sp_org = re.search(r"OS=(\S+ \S+)", prot).group(1) if prot else ""
    rows.append((c, len(seq), reads.get(c, 0.0), cat, gp, round(gc, 2) if gc != "" else "",
                 sprot[c][0] if c in sprot else "", re.sub(r" OS=.*", "", prot.split(" ", 1)[1]) if prot else "",
                 sp_org, ",".join(sorted(pfam.get(c, []))), orf.get(c, "")))
rows.sort(key=lambda r: -r[2])

hdr = ["contig", "length", "reads", "category", "genome_pident", "genome_qcov",
       "sprot_pident", "sprot_protein", "sprot_species", "pfam", "orf_type"]
with open("contig_annotation.tsv", "w") as o:
    o.write("\t".join(hdr) + "\n")
    for r in rows: o.write("\t".join(map(str, r)) + "\n")

tot = sum(r[2] for r in rows)
cat_reads, cat_n = collections.Counter(), collections.Counter()
for r in rows: cat_reads[r[3]] += r[2]; cat_n[r[3]] += 1
with open("category_summary.tsv", "w") as o:
    o.write("category\tcontigs\treads\tpercent_reads\n")
    for k, v in cat_reads.most_common():
        o.write(f"{k}\t{cat_n[k]}\t{v:.0f}\t{v / tot * 100:.1f}\n")

# most abundant contigs without any homology, for web BLAST
with open("top50_contigs_for_web_blast.fasta", "w") as o:
    n = 0
    for r in rows:
        if r[3] in ("ORF only, no homology", "no ORF, no homology") and r[1] >= 300:
            o.write(f">{r[0]}|len={r[1]}|reads={r[2]:.0f}|{r[3].replace(' ', '_').replace(',', '')}\n{contigs[r[0]]}\n")
            n += 1
            if n == 50: break

# 20 most abundant contigs of any category, for web BLAST
with open("top20_abundant_contigs.fasta", "w") as o:
    for r in rows[:20]:
        o.write(f">{r[0]}|len={r[1]}|reads={r[2]:.0f}|{r[3].replace(' ', '_').replace(',', '')}\n{contigs[r[0]]}\n")

# COX1 contigs (COI barcode) for species identity check in BOLD/NCBI
with open("cox1_contigs.fasta", "w") as o:
    for r in rows:
        if r[7].startswith("Cytochrome c oxidase subunit 1") and r[1] >= 400:
            o.write(f">{r[0]}|len={r[1]}|reads={r[2]:.0f}|sprot_{r[8].replace(' ', '_')}_{r[6]}\n{contigs[r[0]]}\n")

print(open("category_summary.tsv").read())
