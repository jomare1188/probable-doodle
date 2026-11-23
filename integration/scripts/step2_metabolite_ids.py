#!/usr/bin/env python3
"""
Step 2 - Map the 18 published metabolite names to KEGG compound, HMDB and ChEBI IDs.

Route A (primary): name -> PubChem CID (PUG REST, exact name/synonym search)
                   -> PubChem cross-references (RegistryID) -> KEGG C#####, HMDB#######, CHEBI:#
Route B (cross-ref):   KEGG REST find/compound with the same query, keeping only entries whose
                   name list contains the query EXACTLY (case-insensitive)
KEGG entries found by either route are compared with KEGG's own ChEBI / PubChem links.
IDs are never invented: a metabolite with no verified match gets empty cells and goes to the
manual table. All web responses are cached in resources/web_cache/step2_pubchem_kegg/ (retrieval date in
RETRIEVED_ON.txt); with the cache in place the script makes no web request.

Query names: the published names, with stereochemistry made explicit where the biology and
the annotation imply it (documented in QUERY_NOTES). Run from the repository root:
    python3 scripts/step2_metabolite_ids.py
"""
import json, re, time, urllib.parse, urllib.request
from datetime import date
from pathlib import Path
import pandas as pd

OUT = Path("results/step2_metabolite_ids"); OUT.mkdir(parents=True, exist_ok=True)
CACHE = Path("resources/web_cache/step2_pubchem_kegg"); CACHE.mkdir(parents=True, exist_ok=True)
if not (CACHE / "RETRIEVED_ON.txt").exists():
    (CACHE / "RETRIEVED_ON.txt").write_text(date.today().isoformat() + "\n")

# published column name -> (query name, note)
QUERIES = {
    "Fusaric acid": ("fusaric acid", ""),
    "Gln-C18:3": ("N-linolenoyl-L-glutamine", "N-acyl glutamine, C18:3 acyl = linolenoyl (assumed alpha-linolenic, 9Z,12Z,15Z)"),
    "Gln-C16:3": ("N-hexadecatrienoyl-L-glutamine", "N-acyl glutamine, C16:3 acyl; double-bond positions not reported"),
    "Gln-C17:3": ("N-heptadecatrienoyl-L-glutamine", "N-acyl glutamine, C17:3 acyl; double-bond positions not reported"),
    "kynurenic acid": ("kynurenic acid", ""),
    "pyroglutamic acid": ("L-pyroglutamic acid", "assumed L (5-oxo-L-proline); LC-MS does not resolve stereochemistry"),
    "soyasaponin Ba": ("soyasaponin Ba", ""),
    "7,3',4'-Trihydroxyflavone": ("3',4',7-trihydroxyflavone", "same compound, flavone numbering"),
    "L-Tryptophan": ("L-tryptophan", ""),
    "p-coumaric acid": ("p-coumaric acid", "trans isomer is the common natural form"),
    "ferulic acid": ("ferulic acid", "trans isomer is the common natural form"),
    "desferrioxamine D2": ("desferrioxamine D2", ""),
    "desferrioxamine E": ("desferrioxamine E", "synonym: nocardamine"),
    "spermidine": ("spermidine", ""),
    "fumonisin B1": ("fumonisin B1", ""),
    "proline": ("L-proline", "assumed L; LC-MS does not resolve stereochemistry"),
    "sinapic acid": ("sinapic acid", ""),
    "pipecolic acid": ("L-pipecolic acid", "assumed L (lysine catabolism); LC-MS does not resolve stereochemistry"),
}


def get(url, key):
    f = CACHE / (re.sub(r"[^A-Za-z0-9]+", "_", key)[:120] + ".txt")
    if f.exists():
        return f.read_text()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                txt = r.read().decode()
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                txt = ""
                break
            time.sleep(2 + attempt * 3)
    else:
        raise RuntimeError(url)
    f.write_text(txt); time.sleep(0.35)       # PubChem asks for <= 5 requests/s
    return txt


def pubchem_cids(name):
    q = urllib.parse.quote(name)
    txt = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/cids/TXT", f"pc_cid_{name}")
    return [int(x) for x in txt.split()]


def pubchem_props(cid):
    txt = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
              "Title,MolecularFormula,MonoisotopicMass,InChIKey/JSON", f"pc_prop_{cid}")
    return json.loads(txt)["PropertyTable"]["Properties"][0] if txt else {}


def pubchem_xrefs(cid):
    txt = get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON",
              f"pc_xref_{cid}")
    if not txt:
        return []
    return json.loads(txt)["InformationList"]["Information"][0].get("RegistryID", [])


def kegg_find(name):
    # KEGG find rejects quotes/commas: search the longest word, then keep exact full-name matches
    word = max(re.split(r"[^A-Za-z0-9]+", name), key=len)
    txt = get(f"https://rest.kegg.jp/find/compound/{urllib.parse.quote(word)}", f"kegg_find_{word}")
    hits = []
    for line in txt.strip().split("\n"):
        if not line:
            continue
        cid, names = line.split("\t")
        if name.lower() in [n.strip().lower() for n in names.split(";")]:
            hits.append(cid.replace("cpd:", ""))
    return hits


def kegg_links(cpd):
    txt = get(f"https://rest.kegg.jp/get/{cpd}", f"kegg_get_{cpd}")
    name = re.search(r"^NAME\s+(.+)$", txt, re.M)
    chebi = re.search(r"ChEBI:\s+(.+)$", txt, re.M)
    pubchem_sid = re.search(r"PubChem:\s+(\d+)", txt)
    return dict(kegg_name=name.group(1).rstrip(";") if name else "",
                kegg_chebi=chebi.group(1).split() if chebi else [],
                kegg_pubchem_sid=pubchem_sid.group(1) if pubchem_sid else "")


rows = []
for published, (query, note) in QUERIES.items():
    cids = pubchem_cids(query)
    cid = cids[0] if cids else None
    props = pubchem_props(cid) if cid else {}
    xr = pubchem_xrefs(cid) if cid else []
    kegg_pc = sorted({x for x in xr if re.fullmatch(r"C\d{5}", x)})
    hmdb = sorted({("HMDB00" + x[4:]) if re.fullmatch(r"HMDB\d{5}", x) else x
                   for x in xr if re.fullmatch(r"HMDB\d{5}(\d{2})?", x)})
    chebi = sorted({x for x in xr if re.fullmatch(r"CHEBI:\d+", x)}, key=lambda s: int(s.split(":")[1]))
    kegg_nm = kegg_find(query)
    kegg_all = sorted(set(kegg_pc) | set(kegg_nm))
    kl = {k: kegg_links(k) for k in kegg_all}
    rows.append(dict(
        metabolite=published, query_name=query, query_note=note,
        pubchem_cid=cid or "", pubchem_n_cids_for_name=len(cids),
        pubchem_title=props.get("Title", ""), formula=props.get("MolecularFormula", ""),
        monoisotopic_mass=props.get("MonoisotopicMass", ""), inchikey=props.get("InChIKey", ""),
        KEGG_from_pubchem=",".join(kegg_pc), KEGG_from_kegg_name_search=",".join(kegg_nm),
        KEGG_names=" | ".join(f"{k}: {v['kegg_name']}" for k, v in kl.items()),
        KEGG_chebi_links=" | ".join(f"{k}: {','.join(v['kegg_chebi'])}" for k, v in kl.items()),
        HMDB_from_pubchem=",".join(hmdb), ChEBI_from_pubchem=",".join(chebi)))

df = pd.DataFrame(rows)
# agreement between routes
def status(r):
    a, b = set(filter(None, r.KEGG_from_pubchem.split(","))), set(filter(None, r.KEGG_from_kegg_name_search.split(",")))
    if not a and not b:
        return "no KEGG"
    if a and b:
        return "agree" if a & b else "CONFLICT"
    return "PubChem only" if a else "KEGG name only"
df["KEGG_route_agreement"] = df.apply(status, axis=1)
df["retrieved"] = (CACHE / "RETRIEVED_ON.txt").read_text().strip()
df.to_csv(OUT / "metabolite_id_auto_mapping.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(df[["metabolite", "query_name", "pubchem_cid", "pubchem_n_cids_for_name", "formula",
          "KEGG_from_pubchem", "KEGG_from_kegg_name_search", "KEGG_route_agreement",
          "HMDB_from_pubchem", "ChEBI_from_pubchem"]].to_string(index=False))
print()
print(df[["metabolite", "pubchem_title", "KEGG_names", "KEGG_chebi_links"]].to_string(index=False))
