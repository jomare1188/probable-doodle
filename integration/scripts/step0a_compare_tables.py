#!/usr/bin/env python3
"""
Step 0a (QC) - compare the cascade genus table with the two earlier tables
(abundance_table_ready.tsv used in the article, abundance_table_genus_CONSOLIDATED.tsv).
Descriptive only: log2 ratio of mean relative abundance (infected / control), pseudocount 0.5 read.
Run from the repository root:  python3 scripts/step0a_compare_tables.py
"""
import numpy as np
import pandas as pd

S = ["control_1", "control_2", "control_3", "infected_1", "infected_2", "infected_3"]
MB = "data/microbiome/"
OUT = "results/step0a_microbiome_cascade/"


def load(path, genus_col=None):
    t = pd.read_csv(path, sep="\t")
    t = t[t["tax"].notna()] if genus_col is None else t
    t["genus"] = t[genus_col] if genus_col else t["tax"].str.split(";").str[-1]
    return t.groupby("genus")[S].sum()


def log2_ratio(t):
    rel = (t + 0.5) / (t + 0.5).sum()
    return np.log2(rel[S[3:]].mean(axis=1) / rel[S[:3]].mean(axis=1))


tabs = {"cascade": load(OUT + "abundance_table_genus_cascade.tsv", "genus"),
        "consolidated": load(MB + "abundance_table_genus_CONSOLIDATED.tsv"),
        "article_ready": load(MB + "abundance_table_ready.tsv")}
cmp = pd.DataFrame({f"log2ratio_{k}": log2_ratio(v) for k, v in tabs.items()})
for k, v in tabs.items():
    cmp[f"reads_{k}"] = v.sum(axis=1)
casc = tabs["cascade"]
cmp["prevalent_cascade"] = ((casc > 0).sum(axis=1) >= 3).reindex(cmp.index)
cmp = cmp.sort_values("reads_cascade", ascending=False)
cmp.to_csv(OUT + "comparison_with_published_tables.tsv", sep="\t")

for k in ["consolidated", "article_ready"]:
    b = cmp.dropna(subset=["log2ratio_cascade", f"log2ratio_{k}"])
    b = b[b["prevalent_cascade"] == True]
    same = (np.sign(b["log2ratio_cascade"]) == np.sign(b[f"log2ratio_{k}"])).mean()
    print(f"cascade vs {k}: {len(b)} shared prevalent genera, same direction {same:.1%}, "
          f"Spearman r = {b['log2ratio_cascade'].rank().corr(b[f'log2ratio_{k}'].rank()):.2f}")
only_old = sorted(set(tabs["consolidated"].index) - set(casc.index))
print(len(only_old), "consolidated genera absent from the cascade, e.g.:", only_old[:25])
pd.set_option("display.width", 200)
print(cmp.head(30).round(2).to_string())
