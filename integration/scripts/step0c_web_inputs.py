#!/usr/bin/env python3
"""
Step 0c - Input files for the web tools (MicrobiomeAnalyst, MetOrigin) -> ./web_inputs/

Genus names for the web lists only (microbiome_DA.csv is unchanged):
  - SILVA strain suffixes are removed ("Coleofasciculus PCC-7420" -> "Coleofasciculus");
  - placeholder clades without a genus name (EUB33-2, UCG-004, "Erysipelotrichaceae UCG-003")
    are left out and listed in web_inputs/ma_taxa_excluded_placeholders.txt;
  - duplicates after cleaning are merged.
Metabolites: KEGG compound IDs (the ID type that mapped best in step 2); IDs never invented.
Run from the repository root:  python3 scripts/step0c_web_inputs.py
"""
import re
from pathlib import Path
import pandas as pd

OUT = Path("web_inputs"); OUT.mkdir(exist_ok=True)
Path("web_results").mkdir(exist_ok=True)
da = pd.read_csv("results/step0b_microbiome_DA/microbiome_DA.csv")
ids = pd.read_csv("results/step2_metabolite_ids/metabolite_ids_final.csv", dtype=str).fillna("")
ids["differential_published"] = ids.differential_published == "True"


def clean(g):
    """-> (web name or None, reason)"""
    if re.fullmatch(r"[A-Z]+\d+[-\w]*|UCG-\d+|.*\bUCG-\d+", g):
        return None, "placeholder clade (no genus name)"
    m = re.fullmatch(r"(Candidatus [A-Z][a-z]+|[A-Z][a-z]+)( .+)?", g)
    if m and m.group(2):
        return m.group(1), "strain suffix removed"
    return g, ""


bact = da[da.differential & (da.kingdom == "Bacteria")].copy()
bact["web_name"], bact["web_note"] = zip(*bact.genus.map(clean))
bact["tier12"] = bact.confidence_tier.str.match(r"^(1|2)")


def write_list(df, fname):
    names = list(dict.fromkeys(df.web_name.dropna()))       # keep order, drop duplicates
    (OUT / fname).write_text("\n".join(names) + "\n")
    return names


t12 = write_list(bact[bact.tier12].sort_values("aldex2_effect_size", key=abs, ascending=False),
                 "ma_taxa_list.txt")
alld = write_list(bact.sort_values("aldex2_effect_size", key=abs, ascending=False),
                  "ma_taxa_list_all_differential.txt")
excl = bact[bact.web_name.isna()][["genus", "confidence_tier", "direction", "web_note"]]
excl.to_csv(OUT / "ma_taxa_excluded_placeholders.txt", sep="\t", index=False)
bact[["genus", "web_name", "web_note", "confidence_tier", "direction",
      "log2_effect_infected_vs_control"]].to_csv(OUT / "ma_taxa_name_mapping.tsv", sep="\t", index=False)

# all bacterial genera detected in the study (prevalence filter), for MetOrigin's optional
# "bacteria detected in your study" field
det = da[(da.kingdom == "Bacteria") & (da.feature_type == "genus")].copy()
det["web_name"] = det.genus.map(lambda g: clean(g)[0])
detected = list(dict.fromkeys(det.sort_values("genus").web_name.dropna()))
(OUT / "metorigin_bacteria_detected.txt").write_text("\n".join(detected) + "\n")
print(f"metorigin_bacteria_detected.txt: {len(detected)} genera")

# metabolites
dif = ids[ids.differential_published]
(OUT / "ma_metabolite_list.txt").write_text("\n".join(dif.KEGG[dif.KEGG != ""]) + "\n")
dif[dif.KEGG == ""][["metabolite", "PubChem_CID", "ChEBI", "HMDB"]].to_csv(
    OUT / "ma_metabolite_list_unmapped.tsv", sep="\t", index=False)
dif[dif.KEGG != ""][["metabolite", "KEGG", "direction"]].to_csv(
    OUT / "ma_metabolite_list_key.tsv", sep="\t", index=False)
mo = pd.DataFrame({"Name": ids.metabolite, "KEGGID": ids.KEGG, "HMDBID": ids.HMDB,
                   "Diff": ids.differential_published.astype(int)})
mo.to_csv(OUT / "metorigin_input.csv", index=False)
# same content in the column order MetOrigin itself uses (its processed file 01_Data/InputCPD.csv
# is HMDBID, KEGGID, Name, Diff); rows with neither a KEGG nor an HMDB ID cannot be matched by
# MetOrigin and are left out of this version (they are handled from the literature in step 4e)
has_id = (mo.KEGGID != "") | (mo.HMDBID != "")
mo.loc[has_id, ["HMDBID", "KEGGID", "Name", "Diff"]].to_csv(OUT / "metorigin_input_MetOrigin_order.csv", index=False)
print("metorigin_input_MetOrigin_order.csv:", has_id.sum(), "metabolites (left out, no ID:",
      ", ".join(mo.Name[~has_id]) + ")")

print(f"ma_taxa_list.txt: {len(t12)} genera (from {bact.tier12.sum()} tier 1-2 bacterial genera)")
print(f"ma_taxa_list_all_differential.txt: {len(alld)} genera (from {len(bact)})")
print("excluded placeholders:", ", ".join(excl.genus))
print("renamed:", "; ".join(f"{a} -> {b}" for a, b in zip(bact.genus, bact.web_name)
                             if b and a != b))
print(f"ma_metabolite_list.txt: {(dif.KEGG != '').sum()} KEGG IDs; unmapped differential: "
      + ", ".join(dif.metabolite[dif.KEGG == ""]))
print("metorigin_input.csv:"); print(mo.to_string(index=False))
