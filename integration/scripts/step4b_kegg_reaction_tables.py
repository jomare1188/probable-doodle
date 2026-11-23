#!/usr/bin/env python3
"""
Step 4b - KEGG reaction knowledge for the study metabolites (cached, reproducible).
For every metabolite with a KEGG compound ID:
  compound -> reactions (KEGG link), reaction equation (substrate / product side as written),
  reaction -> KO (KEGG link), compound -> pathway (KEGG link, reference "map" pathways).
Output: results/step4b_kegg/{compound_reaction_ko.csv, compound_pathway.csv, pathway_names.csv}
KEGG REST responses are cached in resources/web_cache/step4b_kegg/ (retrieval date in
RETRIEVED_ON.txt); with the cache in place the script makes no web request.
Run from the repository root:  python3 scripts/step4b_kegg_reaction_tables.py
"""
import re, time, urllib.request
from pathlib import Path
import pandas as pd

OUT = Path("results/step4b_kegg"); OUT.mkdir(parents=True, exist_ok=True)
CACHE = Path("resources/web_cache/step4b_kegg"); CACHE.mkdir(parents=True, exist_ok=True)
if not (CACHE / "RETRIEVED_ON.txt").exists():
    (CACHE / "RETRIEVED_ON.txt").write_text(time.strftime("%Y-%m-%d") + "\n")


def get(path):
    f = CACHE / (re.sub(r"[^A-Za-z0-9]+", "_", path)[:150] + ".txt")
    if f.exists():
        return f.read_text()
    txt = urllib.request.urlopen(f"https://rest.kegg.jp/{path}", timeout=120).read().decode()
    f.write_text(txt); time.sleep(0.4)
    return txt


def pairs(txt):
    return [l.split("\t") for l in txt.strip().split("\n") if "\t" in l]


ids = pd.read_csv("results/step2_metabolite_ids/metabolite_ids_final.csv", dtype=str).fillna("")
cpds = ids[ids.KEGG != ""][["metabolite", "KEGG"]]

# compound -> reactions
cr = []
for _, r in cpds.iterrows():
    for a, b in pairs(get(f"link/reaction/{r.KEGG}")):
        cr.append((r.metabolite, r.KEGG, b.replace("rn:", "")))
cr = pd.DataFrame(cr, columns=["metabolite", "KEGG", "reaction"])

# reaction equations (10 per request) -> side of the compound
rxns = sorted(cr.reaction.unique()); eq = {}
for i in range(0, len(rxns), 10):
    txt = get("get/" + "+".join(rxns[i:i + 10]))
    for entry in txt.split("///"):
        m = re.search(r"ENTRY\s+(R\d{5})", entry); e = re.search(r"EQUATION\s+(.+)", entry)
        n = re.search(r"NAME\s+(.+)", entry)
        if m and e:
            eq[m.group(1)] = (e.group(1).strip(), n.group(1).strip() if n else "")


def side(c, r):
    e = eq.get(r, ("", ""))[0]
    if "<=>" not in e:
        return "unknown"
    left, right = e.split("<=>")
    return "substrate" if c in left else "product" if c in right else "unknown"


cr["equation"] = cr.reaction.map(lambda r: eq.get(r, ("", ""))[0])
cr["reaction_name"] = cr.reaction.map(lambda r: eq.get(r, ("", ""))[1])
cr["compound_side_as_written"] = [side(c, r) for c, r in zip(cr.KEGG, cr.reaction)]

# reaction -> KO
rk = []
for i in range(0, len(rxns), 50):
    for a, b in pairs(get("link/ko/" + "+".join(rxns[i:i + 50]))):
        rk.append((a.replace("rn:", ""), b.replace("ko:", "")))
rk = pd.DataFrame(rk, columns=["reaction", "KO"]).drop_duplicates()
crk = cr.merge(rk, on="reaction", how="left")
crk.to_csv(OUT / "compound_reaction_ko.csv", index=False)

# compound -> reference pathways (map), overview maps (map011xx, map012xx) excluded
cp = []
for _, r in cpds.iterrows():
    for a, b in pairs(get(f"link/pathway/{r.KEGG}")):
        p = b.replace("path:", "")
        if p.startswith("map") and not p.startswith(("map011", "map012")):
            cp.append((r.metabolite, r.KEGG, p))
pd.DataFrame(cp, columns=["metabolite", "KEGG", "map"]).to_csv(OUT / "compound_pathway.csv", index=False)
pn = pd.DataFrame(pairs(get("list/pathway")), columns=["map", "pathway_name"])
pn["map"] = pn["map"].str.replace("path:", "")
pn.to_csv(OUT / "pathway_names.csv", index=False)
(OUT / "RETRIEVED.txt").write_text(get("info/kegg").split("\n")[1] + "\n" + (CACHE / "RETRIEVED_ON.txt").read_text())

print(f"compounds: {len(cpds)} | reactions: {cr.reaction.nunique()} | reaction-KO pairs: {len(rk)} | "
      f"KOs: {rk.KO.nunique()}")
print(cr.groupby("metabolite").reaction.nunique().to_string())
print("compound-pathway links:", len(cp))
