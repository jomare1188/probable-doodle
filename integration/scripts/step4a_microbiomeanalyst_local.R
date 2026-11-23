#!/usr/bin/env Rscript
# Step 4a - Local reproduction of the MicrobiomeAnalyst microbiome-metabolomics (list input)
# results, so the tables behind the website's images can be saved.
# Uses the SAME reference libraries as the web server (downloaded from
# https://www.microbiomeanalyst.ca/resources/libs/mmp/, see resources/microbiomeanalyst_libs/MD5SUMS.txt) and re-implements
# the functions of MicrobiomeAnalystR R/mmp_utils.R (xia-lab/MicrobiomeAnalystR, master, 2026-09-29):
#   PerformMicNameMap / PerformMetNameMap   -> name matching (validated against the web export)
#   doGemPrediction (list input)            -> AGORA / EMBL genome-scale-model "potential" scores
#   PerformMicNameMap(mtcls) + CreatPathwayLib -> KEGG knowledge network (bacterial pathways)
#   PerformMetListEnrichment                -> network enrichment (adds the pathway names that
#                                              the web export drops: it writes row.names = FALSE)
# Run from the repository root with the integration_diatrea env:
#   Rscript scripts/step4a_microbiomeanalyst_local.R
suppressPackageStartupMessages({ library(qs); library(dplyr); library(tidyr) })
L   <- "resources/microbiomeanalyst_libs/"
out <- "results/step4a_microbiomeanalyst_local"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
web <- "web_results/microbiomeanalyst/microbiomeanalyst_network"
rd  <- function(f) qread(paste0(L, f))
taxalvl <- "genus"; psc <- 0.5          # server defaults (doGemPrediction psc = 0.5)

mic <- readLines("web_inputs/ma_taxa_list.txt"); mic <- mic[mic != ""]
met <- readLines("web_inputs/ma_metabolite_list.txt"); met <- met[met != ""]
key <- read.delim("web_inputs/ma_metabolite_list_key.tsv")
da  <- read.csv("results/step0b_microbiome_DA/microbiome_DA.csv")
map_names <- read.delim("web_inputs/ma_taxa_name_mapping.tsv")

## ------------------------------------------------ 1. name matching (PerformMicNameMap)
mic.map <- data.frame(Query = mic, agora = NA, embl = NA, kegg = NA, ncbi = NA)
ta <- rd("agora_tax.qs")[[taxalvl]]; names(ta)[1] <- "taxa"
mic.map$agora <- ta$taxa[match(mic.map$Query, ta$taxa)]
mic.map$ncbi  <- ta$ncbi_id[match(mic.map$agora, ta$taxa)]
te <- rd("embl_tax.qs")[[taxalvl]]; names(te)[1] <- "taxa"
mic.map$embl  <- te$taxa[match(mic.map$Query, te$taxa)]
mic.map$ncbi[is.na(mic.map$ncbi)] <- te$ncbi_id[match(mic.map$embl[is.na(mic.map$ncbi)], te$taxa)]
tk <- rd("taxMapKEGG.qs")[[taxalvl]]; tki <- tk[["info"]]; names(tki)[1] <- "taxa"
mic.map$kegg  <- tki$taxa[match(mic.map$Query, tki$taxa)]
mic.map$ncbi[is.na(mic.map$ncbi)] <- tki$id[match(mic.map$kegg[is.na(mic.map$ncbi)], tki$taxa)]
webm <- read.csv(file.path(web, "taxa_match_result.csv"), row.names = 1, na.strings = c("", "NA"))
same_mic <- all(mapply(function(a, b) identical(is.na(a), is.na(b)), mic.map[, 2:4], webm[, 2:4]))
cat("Taxa mapping reproduces the web export (agora/embl/kegg columns):", same_mic, "\n")

## ------------------------------------------------ 2. metabolite matching (PerformMetNameMap, KEGG IDs)
g2k <- rd("gem2kegg.qs"); syn <- rd("synonymGem.qs"); k2n <- rd("general_kegg2name.qs")
met.map <- data.frame(Query = met, kegg = met)
for (db in c("agora", "embl")) {
  metdb <- rd(paste0(db, ".met.qs"))
  id <- g2k$metID[match(met, g2k$KEGG)]; id[!id %in% metdb] <- NA
  met.map[[paste0(db, "_id")]] <- id
  met.map[[db]] <- syn$Name[match(id, syn$metID)]
}
met.map$name <- k2n$Name[match(met.map$kegg, k2n$ID)]
met.map$node <- k2n$node[match(met.map$kegg, k2n$ID)]
webmm <- read.csv(file.path(web, "metabolite_match_result.csv"), row.names = 1)
cat("Metabolite mapping reproduces the web export (agora_id/embl_id/node):",
    identical(is.na(met.map$agora_id), is.na(webmm$agora_id) | webmm$agora_id == "") &&
    identical(is.na(met.map$embl_id),  is.na(webmm$embl_id)  | webmm$embl_id == ""), "\n")

## ------------------------------------------------ 3. GEM "potential" scores (AGORA, EMBL)
gem <- lapply(c("agora", "embl"), function(db) {
  m2m <- rd(paste0(db, ".qs"))[[taxalvl]]; names(m2m)[1] <- "taxa"
  ids <- met.map[[paste0(db, "_id")]]
  m2m %>% filter(taxa %in% mic.map[[db]], metID %in% ids) %>%
    mutate(db = toupper(db), kegg = met.map$kegg[match(metID, ids)],
           metabolite_study = key$metabolite[match(kegg, key$KEGG)],
           passes_server_threshold = potential >= psc) %>%
    dplyr::select(db, genus = taxa, kegg, metabolite_study, gem_metID = metID, gem_name = metabolite,
                  potential, passes_server_threshold)
})
gem <- bind_rows(gem)
gem <- gem %>% left_join(da %>% dplyr::select(genus, genus_direction = direction, confidence_tier,
                                              genus_log2_effect = log2_effect_infected_vs_control), by = "genus")
write.csv(gem, file.path(out, "gem_potential_long.csv"), row.names = FALSE)
for (db in c("AGORA", "EMBL")) {
  w <- gem %>% filter(db == !!db) %>% dplyr::select(genus, metabolite_study, potential) %>%
    pivot_wider(names_from = metabolite_study, values_from = potential, values_fill = 0)
  write.csv(w, file.path(out, paste0("gem_potential_matrix_", tolower(db), ".csv")), row.names = FALSE)
  cat(db, ": genera x metabolites =", nrow(w), "x", ncol(w) - 1, "| pairs with potential >=", psc, ":",
      sum(gem$db == db & gem$passes_server_threshold), "\n")
}

## ------------------------------------------------ 4. KEGG knowledge network (bacterial pathways)
# KEGG organisms of each queried genus (taxMapKEGG), their pathways and compounds (bacpathway.met)
bacpath <- rd("bacpathway.met.qs")
orgs <- lapply(setNames(mic, mic), function(g) unique(unlist(tk[[g]])))
orgs <- orgs[lengths(orgs) > 0]
links <- bind_rows(lapply(names(orgs), function(g) {
  bp <- bacpath[intersect(orgs[[g]], names(bacpath))]
  bind_rows(lapply(names(bp), function(o) {
    hits <- lapply(bp[[o]], function(cpds) intersect(cpds, met))
    hits <- hits[lengths(hits) > 0]
    if (!length(hits)) return(NULL)
    data.frame(genus = g, org = o, pathway = rep(names(hits), lengths(hits)), kegg = unlist(hits))
  }))
}))
kn <- links %>% group_by(genus, kegg) %>%
  summarise(n_kegg_organisms = n_distinct(org), n_organisms_of_genus = length(orgs[[first(genus)]]),
            pathways = paste(sort(unique(pathway)), collapse = "; "), .groups = "drop") %>%
  mutate(fraction_organisms = n_kegg_organisms / n_organisms_of_genus,
         metabolite_study = key$metabolite[match(kegg, key$KEGG)]) %>%
  left_join(da %>% dplyr::select(genus, genus_direction = direction, confidence_tier), by = "genus")
write.csv(kn, file.path(out, "kegg_knowledge_links_genus_metabolite.csv"), row.names = FALSE)
cat("KEGG knowledge links: genera with KEGG organisms:", length(orgs), "| genus-metabolite links:", nrow(kn),
    "| metabolites linked:", n_distinct(kn$kegg), "\n")

# network edges as the web viewer builds them (CreatPathwayLib, contain = "usrbac")
mtcls <- unique(unlist(tk[which(names(tk) %in% mic)]))
bp <- bacpath[intersect(mtcls, names(bacpath))]
paths <- unique(unlist(lapply(bp, names)))
current.lib <- setNames(lapply(paths, function(p) unique(unlist(lapply(bp, function(x) x[[p]])))), paths)
edges.bc <- rd("edge.bac.qs")
nodes <- unique(unlist(current.lib))
e <- edges.bc[edges.bc$id_rxn %in% edges.bc$id_rxn[edges.bc$met %in% nodes], ]
e_study <- e[e$id_rxn %in% e$id_rxn[e$met %in% met], ]
write.csv(as.data.frame(e_study), file.path(out, "kegg_network_reaction_edges_study_metabolites.csv"), row.names = FALSE)
cat("KEGG network (usrbac): pathways", length(current.lib), "| compounds", length(nodes),
    "| reaction edges touching study metabolites:", length(unique(e_study$rnIDs)), "\n")

## ------------------------------------------------ 5. network enrichment (PerformMetListEnrichment)
ora <- function(current.set) {
  universe <- unique(unlist(current.set)); q <- intersect(met.map$kegg, universe)
  size <- lengths(current.set); hits <- sapply(current.set, function(x) sum(q %in% x))
  p <- phyper(hits - 1, size, length(universe) - size, length(q), lower.tail = FALSE)
  r <- data.frame(Pathway = names(current.set), Total = size, Expected = length(q) * size / length(universe),
                  Hits = hits, Pval = p, FDR = p.adjust(p, "fdr"),
                  Hit_compounds = sapply(current.set, function(x) paste(intersect(q, x), collapse = ";")))
  r <- r[r$Hits > 0, ]; r[order(r$Pval), ]
}
sets <- list(bac = "kegg_bac_mummichog.qs", hsabac = "kegg_hsa_bac_mummichog.qs",
             hsa = "kegg_hsa_mummichog.qs", all = "kegg_all_mummichog.qs")
webe <- read.csv(file.path(web, "network_enrichment_metabolite_0.csv"))
found <- NA
for (nm in c(names(sets), "usrbac")) {
  cs <- if (nm == "usrbac") current.lib else {
    x <- rd(sets[[nm]])$pathways; setNames(x$cpds, x$name) }
  r <- ora(cs)
  top <- head(r, nrow(webe))
  ok <- nrow(top) == nrow(webe) && all(top$Total == webe$Total) && all(top$Hits == webe$Hits) &&
        isTRUE(all.equal(signif(top$Pval, 3), webe$Pval, tolerance = 1e-2))
  cat(sprintf("enrichment with library '%s': %d pathways with hits; matches web export: %s\n", nm, nrow(r), ok))
  write.csv(r, file.path(out, paste0("network_enrichment_", nm, ".csv")), row.names = FALSE)
  if (ok && is.na(found)) found <- nm
}
cat("Library used by the web run:", found, "\n")
if (!is.na(found)) {
  r <- read.csv(file.path(out, paste0("network_enrichment_", found, ".csv")))
  print(head(r, 10), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
