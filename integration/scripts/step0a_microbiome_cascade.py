#!/usr/bin/env python3
"""
Step 0a - Hierarchical (cascade) taxonomic classification of the saliva nanopore amplicon reads.

Each read is classified by kraken2 against a sequence of databases and is counted ONCE,
at the first database that resolves it:

  stage 1  SILVA 138.2 SSU (custom kraken2 build, named with SILVA's OWN taxonomy)
  stage 2  NCBI targeted loci 16S/18S/28S/ITS (epi2me wf-metagenomics kraken2 db)
  stage 3  Kraken2 PlusPFP-16 (2026-06-26; RefSeq archaea, bacteria, viruses, plasmids,
           human, UniVec_Core, protozoa, fungi and plants)

A read is "resolved" at a stage when its kraken2 assignment has a genus (the assigned
taxon is a genus or lies below one, and the genus name is not a placeholder such as
"uncultured"), or when it is an organelle (SILVA Chloroplast / Mitochondria).
Unresolved reads (unclassified, or classified above genus) go to the next stage.
Reads never resolved keep the first classification they got ("<taxon>_uncl.").

A stricter variant, where only kraken2-UNCLASSIFIED reads go to the next stage, is derived
from the same runs (kraken2 classifies each read independently), as a sensitivity analysis.

Run from the repository root:  python3 scripts/step0a_microbiome_cascade.py
Requires kraken2, seqkit and taxonkit on PATH (conda env integration_diatrea) and the
environment variables set in config.sh:
  RAW_READS_DIR          one sub-folder per sample (control_1 ... infected_3) with the fastq.gz files
  KRAKEN2_DB_SILVA       SILVA 138.2 kraken2 database; its taxonomy/ holds the SILVA names/nodes.dmp
  KRAKEN2_DB_NCBI_LOCI   epi2me NCBI 16S/18S/28S/ITS kraken2 database
  NCBI_LOCI_TAXDUMP      NCBI taxonomy dump used with that database (new_taxdump 2025-01-01)
  KRAKEN2_DB_PLUSPFP16   Kraken2 PlusPFP-16 index 2026-06-26 (its own names/nodes.dmp)
  THREADS                optional, default 32
"""
import os, shutil, subprocess
from collections import Counter
from pathlib import Path
import pandas as pd

KRAKEN2, SEQKIT, TAXONKIT = (shutil.which(t) or t for t in ("kraken2", "seqkit", "taxonkit"))
for tool in (KRAKEN2, SEQKIT, TAXONKIT):
    assert Path(tool).exists(), f"{tool} not found on PATH"
ENVV = {k: os.environ.get(k, "") for k in ["RAW_READS_DIR", "KRAKEN2_DB_SILVA", "KRAKEN2_DB_NCBI_LOCI",
                                          "NCBI_LOCI_TAXDUMP", "KRAKEN2_DB_PLUSPFP16"]}
missing = [k for k, v in ENVV.items() if not v or not Path(v).is_dir()]
assert not missing, f"set these variables in config.sh (existing folders): {', '.join(missing)}"
RAW = Path(ENVV["RAW_READS_DIR"])
SAMPLES = ["control_1", "control_2", "control_3", "infected_1", "infected_2", "infected_3"]
MIN_LEN, MAX_LEN = 100, 2000          # same read-length filter as the original SILVA 138.2 run
CONFIDENCE = 0.0                      # kraken2 default, same as all original runs
THREADS = int(os.environ.get("THREADS", "32"))
STAGES = [
    dict(name="stage1_silva138.2", label="SILVA 138.2",
         db=Path(ENVV["KRAKEN2_DB_SILVA"]), taxdump=Path(ENVV["KRAKEN2_DB_SILVA"]) / "taxonomy"),
    dict(name="stage2_ncbi_16s_18s_28s_ITS", label="NCBI 16S/18S/28S/ITS",
         db=Path(ENVV["KRAKEN2_DB_NCBI_LOCI"]), taxdump=Path(ENVV["NCBI_LOCI_TAXDUMP"])),
    dict(name="stage3_pluspfp16_20260626", label="PlusPFP-16 2026-06-26",
         db=Path(ENVV["KRAKEN2_DB_PLUSPFP16"]), taxdump=Path(ENVV["KRAKEN2_DB_PLUSPFP16"])),
]
RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"]
PLACEHOLDER = ("uncultured", "unidentified", "unknown", "incertae sedis", "metagenome",
               "unclassified", "environmental")

OUT = Path("results/step0a_microbiome_cascade")
READS = OUT / "reads"
OUT.mkdir(parents=True, exist_ok=True); READS.mkdir(exist_ok=True)


def run(cmd, **kw):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(list(map(str, cmd)), check=True, **kw)


def taxonkit_dir(stage):
    """taxonkit needs names/nodes/merged/delnodes; the SILVA taxonomy has only the first two."""
    d = OUT / "taxdump" / stage["name"]
    d.mkdir(parents=True, exist_ok=True)
    for f in ["names.dmp", "nodes.dmp", "merged.dmp", "delnodes.dmp"]:
        src, dst = stage["taxdump"] / f, d / f
        if dst.exists() or dst.is_symlink():
            continue
        if src.exists():
            dst.symlink_to(src)
        else:
            dst.touch()
    return d


def lineage(taxids, stage):
    """taxid -> {rank: name} for the ranks in RANKS, using the stage's own taxonomy."""
    taxids = sorted({t for t in taxids if t != 0})
    if not taxids:
        return {}
    fmt = "\t".join("{" + r + "}" for r in ["domain|superkingdom"] + RANKS[1:])
    res = subprocess.run([str(TAXONKIT), "reformat2", "--data-dir", str(taxonkit_dir(stage)),
                          "-I", "1", "-f", fmt, "-r", ""],
                         input="\n".join(map(str, taxids)), capture_output=True, text=True, check=True)
    out = {}
    for line in res.stdout.rstrip("\n").split("\n"):
        f = line.split("\t")
        out[int(f[0])] = dict(zip(RANKS, f[1:] + [""] * (len(RANKS) - len(f[1:]))))
    return out


def is_placeholder(name):
    return name == "" or any(p in name.lower() for p in PLACEHOLDER)


def status(lin):
    names = " ".join(lin.values())
    if "Chloroplast" in names or "Mitochondria" in names:
        return "organelle"
    return "genus" if not is_placeholder(lin["genus"]) else "above_genus"


# ------------------------------------------------------------------ 1. read set
stats = []
for s in SAMPLES:
    fq = READS / f"{s}.fastq.gz"
    if not fq.exists():
        files = sorted((RAW / s).glob("*.fastq.gz"))
        p1 = subprocess.Popen(["zcat", *map(str, files)], stdout=subprocess.PIPE)
        subprocess.run([str(SEQKIT), "seq", "-m", str(MIN_LEN), "-M", str(MAX_LEN),
                        "-o", str(fq), "-"], stdin=p1.stdout, check=True)
        p1.wait()
    n_raw = int(subprocess.run(f"zcat {RAW / s}/*.fastq.gz | awk 'END{{print NR/4}}'", shell=True,
                               capture_output=True, text=True, check=True).stdout)
    n_kept = int(subprocess.run(f"zcat {fq} | awk 'END{{print NR/4}}'", shell=True,
                                capture_output=True, text=True, check=True).stdout)
    stats.append(dict(sample=s, raw_reads=n_raw, reads_100_2000bp=n_kept))
    print(s, n_raw, n_kept, flush=True)

# ------------------------------------------------------------------ 2. cascade
per_read = []                        # one row per read and stage
pending = {s: READS / f"{s}.fastq.gz" for s in SAMPLES}
for i, st in enumerate(STAGES, 1):
    sd = OUT / st["name"]; sd.mkdir(exist_ok=True)
    for s in SAMPLES:
        kout, krep = sd / f"{s}.kraken2.tsv.gz", sd / f"{s}.kraken2.report.txt"
        if not kout.exists():
            tmp = sd / f"{s}.kraken2.tsv"
            run([KRAKEN2, "--db", st["db"], "--threads", THREADS, "--confidence", CONFIDENCE,
                 "--gzip-compressed", "--output", tmp, "--report", krep, pending[s]])
            run(["gzip", "-f", tmp])
        k = pd.read_csv(kout, sep="\t", header=None, usecols=[0, 1, 2, 3],
                        names=["cu", "read_id", "taxid", "length"], compression="gzip")
        k["sample"], k["stage"], k["stage_n"] = s, st["name"], i
        per_read.append(k)
    # lineages for every taxid seen at this stage, in this stage's own taxonomy
    stage_df = pd.concat([d for d in per_read if d["stage_n"].iat[0] == i])
    lin = lineage(stage_df["taxid"].unique(), st)
    for d in per_read:
        if d["stage_n"].iat[0] != i:
            continue
        lin_rows = [lin.get(t, dict.fromkeys(RANKS, "")) for t in d["taxid"]]
        for r in RANKS:
            d[r] = [x[r] for x in lin_rows]
        d["status"] = ["unclassified" if c == "U" else status(x) for c, x in zip(d["cu"], lin_rows)]
        s = d["sample"].iat[0]
        unresolved = d.loc[d["status"].isin(["unclassified", "above_genus"]), "read_id"]
        if i < len(STAGES):
            ids = sd / f"{s}.unresolved_ids.txt"
            unresolved.to_csv(ids, index=False, header=False)
            nxt = sd / f"{s}.unresolved.fastq.gz"
            if not nxt.exists():
                run([SEQKIT, "grep", "-f", ids, "-o", nxt, pending[s]])
            pending[s] = nxt

per_read = pd.concat(per_read, ignore_index=True)


# ------------------------------------------------------------------ 3. final assignment per read
def label(row):
    if row["status"] == "unclassified":
        return "Unclassified"
    if row["status"] in ("genus", "organelle") and row["genus"]:
        return row["genus"]
    deepest = next((row[r] for r in reversed(RANKS) if row[r] and not is_placeholder(row[r])), "")
    return f"{deepest}_uncl." if deepest else "Unclassified"


def category(row):
    names = [row[r] for r in RANKS]
    if row["status"] == "unclassified" or not any(names):
        return "unclassified"
    if row["status"] == "organelle":
        return "organelle (chloroplast/mitochondria)"
    if row["superkingdom"] in ("Bacteria", "Archaea"):
        return row["superkingdom"]
    if "Fungi" in names:
        return "Fungi"
    if row["superkingdom"] == "Viruses":
        return "Viruses"
    # SILVA names eukaryote supergroups (e.g. Nucletmycea, Archaeplastida) at this rank
    return "Eukaryota (non-fungal)"


def assign(first_resolving):
    """first_resolving: DataFrame of per-read rows; returns one row per read."""
    first_resolving = first_resolving.copy()
    first_resolving["taxon"] = first_resolving.apply(label, axis=1)
    first_resolving["category"] = first_resolving.apply(category, axis=1)
    return first_resolving


resolved = per_read[per_read["status"].isin(["genus", "organelle"])]
first_res = resolved.sort_values("stage_n").drop_duplicates(["sample", "read_id"])
# reads never resolved: keep their first classification (any rank), else unclassified
never = per_read[~per_read.set_index(["sample", "read_id"]).index.isin(
    first_res.set_index(["sample", "read_id"]).index)]
never_cls = never[never["status"] != "unclassified"].sort_values("stage_n") \
    .drop_duplicates(["sample", "read_id"])
never_u = never[~never.set_index(["sample", "read_id"]).index.isin(
    never_cls.set_index(["sample", "read_id"]).index)].drop_duplicates(["sample", "read_id"])
final = assign(pd.concat([first_res, never_cls, never_u], ignore_index=True))
assert final.duplicated(["sample", "read_id"]).sum() == 0

# strict variant: a read moves on only if kraken2 left it UNCLASSIFIED
strict = per_read[per_read["cu"] == "C"].sort_values("stage_n").drop_duplicates(["sample", "read_id"])
strict_u = per_read[~per_read.set_index(["sample", "read_id"]).index.isin(
    strict.set_index(["sample", "read_id"]).index)].drop_duplicates(["sample", "read_id"])
strict = assign(pd.concat([strict, strict_u], ignore_index=True))

final[["sample", "read_id", "length", "stage", "taxid", "status", "category", "taxon"] + RANKS] \
    .to_csv(OUT / "per_read_final_assignment.tsv.gz", sep="\t", index=False, compression="gzip")


# ------------------------------------------------------------------ 4. genus tables
def genus_table(df):
    microbial = df[df["category"].isin(["Bacteria", "Archaea", "Fungi"])]
    counts = microbial.pivot_table(index="taxon", columns="sample", values="read_id",
                                   aggfunc="count", fill_value=0).reindex(columns=SAMPLES, fill_value=0)
    # lineage of each taxon: from the stage that contributed most of its reads
    top = (microbial.groupby(["taxon", "stage"]).size().rename("n").reset_index()
           .sort_values("n", ascending=False).drop_duplicates("taxon"))
    lin = microbial.drop_duplicates(["taxon", "stage"]).merge(top[["taxon", "stage"]])
    lin = lin.set_index("taxon")
    kingdom = microbial.groupby("taxon")["category"].agg(lambda x: x.value_counts().index[0])
    stages = microbial.groupby("taxon")["stage"].agg(lambda x: ";".join(
        f"{k}:{v}" for k, v in sorted(Counter(x).items())))
    tab = pd.DataFrame({
        "tax": [";".join([kingdom[t]] + [lin.loc[t, r] or "NA" for r in
                                          ["phylum", "class", "order", "family"]] + [t])
                for t in counts.index],
        "genus": counts.index, "kingdom": kingdom.reindex(counts.index).values,
        "phylum": [lin.loc[t, "phylum"] for t in counts.index],
        "reads_by_stage": stages.reindex(counts.index).values})
    tab = pd.concat([tab.reset_index(drop=True), counts.reset_index(drop=True)], axis=1)
    tab["total"] = tab[SAMPLES].sum(axis=1)
    return tab.sort_values("total", ascending=False)


tab = genus_table(final)
tab_strict = genus_table(strict)
tab.to_csv(OUT / "abundance_table_genus_cascade.tsv", sep="\t", index=False)
tab_strict.to_csv(OUT / "abundance_table_genus_cascade_strict.tsv", sep="\t", index=False)


# ------------------------------------------------------------------ 5. read accounting
def accounting(df, name):
    a = df.groupby(["sample", "category"]).size().unstack(fill_value=0).reindex(SAMPLES)
    a.insert(0, "variant", name)
    return a


acc = pd.concat([accounting(final, "genus-rule cascade (primary)"),
                 accounting(strict, "unclassified-only cascade (sensitivity)")])
by_stage = final.groupby(["sample", "stage", "status"]).size().unstack(fill_value=0)
st = pd.DataFrame(stats).set_index("sample")
acc.to_csv(OUT / "read_accounting_by_category.tsv", sep="\t")
by_stage.to_csv(OUT / "read_accounting_by_stage.tsv", sep="\t")
st.to_csv(OUT / "read_set_filtering.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("\n== read filtering\n", st)
print("\n== final read category per sample\n", acc)
print("\n== where reads were resolved (primary)\n", by_stage)
print("\n== genera: primary", len(tab), "| strict", len(tab_strict))
print(tab.groupby("kingdom").size())
print(tab.head(30).drop(columns="tax").to_string(index=False))
