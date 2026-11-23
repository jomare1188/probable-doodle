#!/usr/bin/env python3
"""
Step 2b - Final metabolite ID table + manual table for the metabolites without a KEGG compound.

Acceptance rules (applied to the output of step2_metabolite_ids.py):
  KEGG  : accepted if found by both routes (PubChem xref and exact KEGG name search), or by one
          route when KEGG's own ChEBI link matches a PubChem ChEBI ID. PubChem xrefs that fail
          this rule are rejected and reported (e.g. kynurenic acid -> C14645, an estradiol ester).
  HMDB  : PubChem cross-reference of the accepted PubChem record.
  ChEBI : PubChem cross-reference; otherwise ChEBI search (exact name), noted as such.
Record-level overrides (documented in OVERRIDES) choose a more specific PubChem record.
KEGG pathway links are read from the step 2 web cache (fetched only if missing).
Run from the repository root:  python3 scripts/step2b_finalize_metabolite_ids.py
"""
import urllib.request
from pathlib import Path
import pandas as pd

OUT = Path("results/step2_metabolite_ids")
CACHE = Path("resources/web_cache/step2_pubchem_kegg")
auto = pd.read_csv(OUT / "metabolite_id_auto_mapping.csv", dtype=str).fillna("")
summary = pd.read_csv("results/step0/metabolome_summary_infected_vs_control.csv")

# ------------------------------------------------------------------ record-level overrides
# (each verified by hand against PubChem / KEGG / ChEBI on the retrieval date; formula and
# InChIKey are those of the chosen PubChem record)
OVERRIDES = {
    "sinapic acid": dict(
        pubchem_cid="637775", HMDB="HMDB0032616", ChEBI="CHEBI:15714", KEGG="C00482",
        formula="C11H12O5", InChIKey="PCMORTLOPMLEFB-ONEGZZNKSA-N",
        why="PubChem CID 10743 (stereo unspecified) has no KEGG/HMDB xrefs; the trans record "
            "CID 637775 cross-references KEGG C00482, HMDB0032616 and CHEBI:15714, consistent with "
            "the trans records used for p-coumaric (CID 637542) and ferulic acid (CID 445858)."),
    "Gln-C18:3": dict(
        pubchem_cid="5470079", ChEBI="CHEBI:186604",
        formula="C23H38N2O4", InChIKey="FEFAAEYMDHYUFL-PDBXOOCHSA-N",
        why="CID 5470078 returned by the name search is the ammonium salt (C23H41N3O4); its parent "
            "compound CID 5470079 (C23H38N2O4) is used. ChEBI 'N-linolenoyl glutamine' "
            "(CHEBI:186604, C23H38N2O4) from a ChEBI name search. No KEGG or HMDB entry."),
    "desferrioxamine D2": dict(
        why="PubChem CID 14671620 (title 'desmethylenylnocardamine') lists 'Desferrioxamine D2' and "
            "'Proferrioxamine-D2' as synonyms; formula C26H46N6O9. No KEGG or HMDB entry."),
    "fumonisin B1": dict(
        KEGG="C19241",
        why="KEGG C19241 found by exact name; PubChem has no KEGG xref, but KEGG's ChEBI link "
            "(CHEBI:38221) equals PubChem's ChEBI -> accepted."),
    "pipecolic acid": dict(
        KEGG="C00408",
        why="KEGG C00408 (L-pipecolate) from PubChem xref; KEGG's ChEBI links include PubChem's "
            "CHEBI:30913 -> accepted."),
    "kynurenic acid": dict(
        KEGG="C01717",
        why="PubChem xrefs list C01717 and C14645; C14645 is estradiol-17-phenylpropionate (wrong "
            "xref) and was rejected. C01717 is confirmed by exact KEGG name and ChEBI:18344."),
}

rows = []
for _, r in auto.iterrows():
    m = r.metabolite
    kegg_pc = set(filter(None, r.KEGG_from_pubchem.split(",")))
    kegg_nm = set(filter(None, r.KEGG_from_kegg_name_search.split(",")))
    chebi_pc = set(filter(None, r.ChEBI_from_pubchem.split(",")))
    links = {k.split(": ")[0]: set("CHEBI:" + x for x in k.split(": ")[1].split(",") if x)
             for k in r.KEGG_chebi_links.split(" | ") if ": " in k}
    accepted = {k for k in kegg_pc | kegg_nm
                if (k in kegg_pc and k in kegg_nm) or (links.get(k, set()) & chebi_pc)}
    rejected = (kegg_pc | kegg_nm) - accepted
    o = OVERRIDES.get(m, {})
    kegg = o.get("KEGG", ",".join(sorted(accepted)))
    chebi = o.get("ChEBI", sorted(chebi_pc, key=lambda s: int(s.split(":")[1]))[0] if chebi_pc else "")
    rows.append(dict(
        metabolite=m, query_name=r.query_name, query_note=r.query_note,
        PubChem_CID=o.get("pubchem_cid", r.pubchem_cid), formula=o.get("formula", r.formula),
        InChIKey=o.get("InChIKey", r.inchikey),
        KEGG=kegg, HMDB=o.get("HMDB", r.HMDB_from_pubchem.split(",")[0] if r.HMDB_from_pubchem else ""),
        ChEBI=chebi, KEGG_rejected=",".join(sorted(rejected - set(kegg.split(",")))),
        curation_note=o.get("why", ""), retrieved=r.retrieved))
ids = pd.DataFrame(rows)

# KEGG pathways of each accepted compound (map IDs)
def kegg_link_pathway(c):
    f = CACHE / f"kegg_link_pathway_{c}.txt"
    if not f.exists():
        f.write_text(urllib.request.urlopen(f"https://rest.kegg.jp/link/pathway/{c}", timeout=60).read().decode())
    return f.read_text()


def pathways(c):
    if not c:
        return ""
    txt = kegg_link_pathway(c)
    return ",".join(sorted({l.split("\t")[1].replace("path:", "") for l in txt.strip().split("\n")
                            if "\tpath:map" in l and not l.split("\tpath:")[1].startswith(("map011", "map012"))}))
ids["KEGG_pathways"] = [pathways(k.split(",")[0]) for k in ids.KEGG]

ids = ids.merge(summary[["metabolite", "differential_published", "direction",
                         "detected_only_in_infected", "log2FC_mean_infected_vs_control"]], on="metabolite")
ids.to_csv(OUT / "metabolite_ids_final.csv", index=False)

# ------------------------------------------------------------------ manual table
MANUAL = {
    "Gln-C18:3": dict(
        compound_class="N-acyl glutamine (fatty acid-amino acid conjugate, FAC)",
        probable_origin="co-metabolism: plant-derived linolenic acid conjugated to glutamine in the "
                        "caterpillar gut; gut bacteria can also synthesise N-acyl glutamines",
        related_KEGG="C06427 alpha-linolenic acid; C00064 L-glutamine; C16345 volicitin "
                     "(N-(17-hydroxylinolenoyl)-L-glutamine), the best-known insect FAC",
        related_pathway="map00592 alpha-linolenic acid metabolism; map00250 Ala/Asp/Glu metabolism",
        handling="literature-based (step 4e)"),
    "Gln-C16:3": dict(
        compound_class="N-acyl glutamine (FAC), C16:3 acyl chain",
        probable_origin="co-metabolism (as Gln-C18:3); C16:3 fatty acids are typical of plant "
                        "chloroplast lipids (e.g. 7Z,10Z,13Z-hexadecatrienoic acid)",
        related_KEGG="C00064 L-glutamine (no KEGG entry for the C16:3 acid or conjugate)",
        related_pathway="map00250 Ala/Asp/Glu metabolism",
        handling="literature-based (step 4e)"),
    "Gln-C17:3": dict(
        compound_class="N-acyl glutamine (FAC), odd-chain C17:3 acyl",
        probable_origin="uncertain; odd-chain acyl groups point to microbial fatty acids or to "
                        "oxidation of C18:3 - annotation should be treated as putative",
        related_KEGG="C00064 L-glutamine",
        related_pathway="map00250 Ala/Asp/Glu metabolism",
        handling="literature-based (step 4e)"),
    "soyasaponin Ba": dict(
        compound_class="triterpenoid saponin (group B soyasaponin)",
        probable_origin="external: diet (sugarcane; larvae fed only sugarcane with or without the "
                        "fungus). Soyasaponins are typical of legumes, not sugarcane, so the "
                        "annotation is doubtful (possibly an isobaric triterpenoid saponin)",
        related_KEGG="C08980 soyasapogenol B (aglycone); C08983 soyasaponin I (group B)",
        related_pathway="map00909 sesquiterpenoid and triterpenoid biosynthesis (plant)",
        handling="literature-based (step 4e); diet-derived, identity uncertain"),
    "7,3',4'-Trihydroxyflavone": dict(
        compound_class="flavone (plant flavonoid)",
        probable_origin="external: diet (sugarcane flavonoid); possible microbial/host modification",
        related_KEGG="no KEGG entry; isomers in KEGG (e.g. C10037 3,7,4'-trihydroxyflavone) are "
                     "different compounds and are NOT used",
        related_pathway="map00944 flavone and flavonol biosynthesis (plant)",
        handling="literature-based (step 4e); diet (sugarcane) marker"),
    "desferrioxamine D2": dict(
        compound_class="hydroxamate siderophore (desferrioxamine family)",
        probable_origin="microbial: desferrioxamines are made by actinobacteria (e.g. Streptomyces); "
                        "some fungi and bacteria take them up as xenosiderophores",
        related_KEGG="C23094 desferrioxamine E; C06940 desferrioxamine B",
        related_pathway="map00975 biosynthesis of siderophore group nonribosomal peptides",
        handling="linked to desferrioxamine E / siderophore pathway; literature-based (step 4e)"),
}
man = ids[ids.metabolite.isin(MANUAL)][["metabolite", "PubChem_CID", "formula", "ChEBI", "HMDB",
                                        "differential_published", "direction"]].copy()
for k in ["compound_class", "probable_origin", "related_KEGG", "related_pathway", "handling"]:
    man[k] = man.metabolite.map(lambda m: MANUAL[m][k])
# the published metabolomics gives no MS/MS or library-level confirmation for these annotations
man["annotation_confidence"] = "putative (no MS/MS or library confirmation reported)"
man.to_csv(OUT / "metabolite_manual_table_unmapped.csv", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(ids[["metabolite", "PubChem_CID", "KEGG", "HMDB", "ChEBI", "KEGG_rejected", "KEGG_pathways",
           "differential_published"]].to_string(index=False))
n = len(ids)
for col in ["KEGG", "HMDB", "ChEBI", "PubChem_CID"]:
    ok = (ids[col] != "").sum()
    dif = ((ids[col] != "") & ids.differential_published).sum()
    print(f"{col:12s}: {ok}/{n} mapped ({ok/n:.0%}); differential metabolites mapped: "
          f"{dif}/{ids.differential_published.sum()}")
print("KEGG compounds with >=1 KEGG pathway:", (ids.KEGG_pathways != "").sum())
