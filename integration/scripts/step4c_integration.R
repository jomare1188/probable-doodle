#!/usr/bin/env Rscript
# Step 4 - Knowledge-based, summary-level integration (independent cohorts).
#   a) metabolite-centric evidence table
#   b) pathway x omics matrix + heatmap + UpSet (descriptive; NO combined p with the microbiome)
#   c) direction agreement along gene - metabolite - taxon chains
#   d) Cytoscape node / edge tables (edges only from KEGG reactions, MicrobiomeAnalyst and,
#      in a separate edge type, the literature) + static network figure
#   e) literature-based handling of unmapped metabolites and fungi (positive control:
#      Fusarium + fumonisin B1 + fusaric acid)
# Inputs: steps 0-3, step 4a (MicrobiomeAnalyst local), 4b (KEGG reactions), MetOrigin exports.
# Run from the repository root with the integration_diatrea env:  Rscript scripts/step4c_integration.R
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(openxlsx); library(ggplot2); library(igraph); library(ggraph)
  library(ComplexHeatmap); library(UpSetR); library(grid)
})
set.seed(20260930)
out <- "results/step4_integration"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
mo_dir <- "web_results/metorigin/fruit_fly_host"           # primary MetOrigin run (D. melanogaster host, as in step 3)
mo_dir_h <- "web_results/metorigin/human_host_default"
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
pasteu <- function(x) { x <- unique(na.omit(x)); if (length(x)) paste(x, collapse = "; ") else "" }

## =============================================================== inputs
deseq <- read.csv("results/step0/deseq2_all_genes_infected_vs_control.csv")
genes <- read.delim("results/step1_gene_annotation/gene_annotation_all_genes.tsv") %>% filter(in_deseq2)
msum  <- read.csv("results/step0/metabolome_summary_infected_vs_control.csv", check.names = FALSE)
mids  <- read.csv("results/step2_metabolite_ids/metabolite_ids_final.csv", colClasses = "character")
da    <- read.csv("results/step0b_microbiome_DA/microbiome_DA.csv")
mg    <- read.csv("results/step3_multigsea/multigsea_kegg_dmel.csv")
crk   <- read.csv("results/step4b_kegg/compound_reaction_ko.csv")
cp    <- read.csv("results/step4b_kegg/compound_pathway.csv")
pn    <- read.csv("results/step4b_kegg/pathway_names.csv")
gem   <- read.csv("results/step4a_microbiomeanalyst_local/gem_potential_long.csv")
kkn   <- read.csv("results/step4a_microbiomeanalyst_local/kegg_knowledge_links_genus_metabolite.csv")
maen  <- read.csv("results/step4a_microbiomeanalyst_local/network_enrichment_usrbac.csv")

met <- msum %>% dplyr::select(metabolite, log2FC = log2FC_mean_infected_vs_control, direction,
                              only_infected = detected_only_in_infected, differential = differential_published,
                              published_p_adj) %>%
  left_join(mids %>% dplyr::select(metabolite, KEGG, HMDB, ChEBI, PubChem_CID), by = "metabolite")
met$KEGG[is.na(met$KEGG)] <- ""
sgn <- function(d) ifelse(d == "up_in_infected", 1, -1)
met$sign <- sgn(met$direction)
key_kegg <- setNames(met$metabolite, met$KEGG)[met$KEGG != ""]

da$tier12 <- da$differential & grepl("^(1|2)", da$confidence_tier)
dif_genera <- da %>% filter(differential) %>%
  dplyr::select(genus, kingdom, genus_direction = direction, genus_log2_effect = log2_effect_infected_vs_control,
                confidence_tier, tier12, contaminant_flag)
mapnum <- function(x) sub(".*?(\\d{5})$", "\\1", x)

## =============================================================== MetOrigin
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
read_origin <- function(d) read.csv(file.path(d, "02_Origin_Analysis/Metabolites_Origin.csv")) %>%
  mutate(key = norm(Name))
mo  <- read_origin(mo_dir); moh <- read_origin(mo_dir_h)
met$key <- norm(met$metabolite)
met <- met %>% left_join(mo %>% dplyr::select(key, MetOrigin_origin_fly = Origin), by = "key") %>%
  left_join(moh %>% dplyr::select(key, MetOrigin_origin_human = Origin), by = "key")
mpea <- bind_rows(lapply(c("Host", "Microbiota", "Co-Metabolism"), function(k) {
  f <- file.path(mo_dir, "03_Pathway_Analysis", paste0("MPEA_Results_", k, ".csv"))
  d <- read.csv(f, check.names = FALSE)
  if (!nrow(d)) return(NULL)
  data.frame(category = k, map = mapnum(d$ID), name = d$Name, p = as.numeric(d$Pvalue),
             metabolites = d$`Metabolite ID`)
}))
sank_files <- list.files(file.path(mo_dir, "05_Sankey_Network"), pattern = "xlsx$", recursive = TRUE, full.names = TRUE)
sank <- bind_rows(lapply(sank_files, function(f) read.xlsx(f) %>%
                           mutate(category = basename(dirname(f)), across(everything(), as.character))))
sank_long <- bind_rows(lapply(names(key_kegg), function(k) {
  s <- sank %>% filter(grepl(k, SubstrateName, fixed = TRUE) | grepl(k, ProductName, fixed = TRUE))
  if (!nrow(s)) return(NULL)
  s %>% mutate(KEGG = k, role = ifelse(grepl(k, SubstrateName, fixed = TRUE), "substrate", "product"))
}))
sank_genus <- sank_long %>% filter(!is.na(genus)) %>%
  group_by(genus, kingdom, KEGG, map = mapnum(PathwayID), reaction = ReactionName, ReactionType, role, KO) %>%
  summarise(n_species = n_distinct(species), .groups = "drop") %>%
  mutate(metabolite = key_kegg[KEGG])
sank_dif <- sank_genus %>% inner_join(dif_genera, by = "genus")

## =============================================================== reaction-linked Diatraea genes
gk <- genes %>% separate_rows(KEGG_ko, sep = ",") %>% filter(!is.na(KEGG_ko), KEGG_ko != "")
rx_genes <- crk %>% filter(!is.na(KO)) %>%
  inner_join(gk %>% dplyr::select(gene, KO = KEGG_ko, preferred_name, dmel_symbol, log2FC_infected_vs_control, padj, DE),
             by = "KO", relationship = "many-to-many") %>%
  mutate(gene_direction = ifelse(log2FC_infected_vs_control > 0, "up_in_infected", "down_in_infected"))
write.csv(rx_genes, file.path(out, "gene_metabolite_kegg_reaction_links.csv"), row.names = FALSE)

## =============================================================== (e) literature-based links
# Curated by the analyst from the publications cited in the "reference" column
# (first author, year, journal volume:pages); compared with bibliographic records on 2026-11-23.
# Lait et al. 2003 found no evidence of microbial N-linolenoyl-glutamine synthesis in Manduca sexta, so the
# bacterial contribution to the N-acyl glutamines is "possible", not established.
# No sugarcane report exists for 7,3',4'-trihydroxyflavone (class-level only) or for soyasaponin Ba.
lit <- tribble(
  ~source, ~source_type, ~target, ~relation, ~reference, ~note,
  "Fusarium", "fungus", "fumonisin B1", "produces",
  "Gelderblom et al. 1988, Appl Environ Microbiol 54:1806-1811",
  "Fumonisins were first isolated from F. moniliforme (= F. verticillioides). POSITIVE CONTROL",
  "Fusarium", "fungus", "Fusaric acid", "produces",
  "Bacon et al. 1996, Appl Environ Microbiol 62:4039-4043",
  "Fusaric acid is produced by many Fusarium species, including F. moniliforme (= F. verticillioides). POSITIVE CONTROL",
  "Streptomyces", "bacterium", "desferrioxamine E", "produces",
  "Barona-Gomez et al. 2004, J Am Chem Soc 126:16282-16283",
  "Desferrioxamine E is the main desferrioxamine of Streptomyces coelicolor; des gene cluster (desA-D)",
  "Streptomyces", "bacterium", "desferrioxamine D2", "produces",
  "Bickel et al. 1960, Helv Chim Acta 43:2118-2128; Keller-Schierlein et al. 1965, Helv Chim Acta 48:710-723",
  "Ferrioxamine D2 was isolated as a minor siderophore of Streptomyces cultures and its structure determined",
  "Microbacterium", "bacterium", "Gln-C18:3", "may produce (N-acyl glutamines)",
  "Spiteller et al. 2000, Biol Chem 381:755-762; Ping et al. 2007, Environ Microbiol 9:1572-1583",
  "Lepidopteran gut bacteria synthesise N-acyl glutamines from added precursors; the Spodoptera exigua isolate Microbacterium arborescens SE14 carries an enzyme that hydrolyses and synthesises N-acyl amino acids",
  "Microbacterium", "bacterium", "Gln-C16:3", "may produce (N-acyl glutamines, class-level)",
  "Spiteller et al. 2000, Biol Chem 381:755-762; Ping et al. 2007, Environ Microbiol 9:1572-1583",
  "class-level inference",
  "Microbacterium", "bacterium", "Gln-C17:3", "may produce (N-acyl glutamines, class-level)",
  "Spiteller et al. 2000, Biol Chem 381:755-762; Ping et al. 2007, Environ Microbiol 9:1572-1583",
  "class-level inference; C17:3 annotation putative",
  "Diatraea (host)", "host", "Gln-C18:3", "conjugates plant fatty acids with glutamine",
  "Lait et al. 2003, PNAS 100:7027-7032; Yoshinaga et al. 2008, PNAS 105:18058-18063; Alborn et al. 1997, Science 276:945-949",
  "N-linolenoyl-L-glutamine (= Gln-C18:3) is synthesised by caterpillar gut tissue from plant linolenic acid and midgut glutamine (no evidence of microbial synthesis in Manduca sexta); volicitin is its 17-hydroxy analogue",
  "sugarcane (diet)", "diet", "p-coumaric acid", "cell-wall hydroxycinnamate (released)",
  "Masarin et al. 2011, Biotechnol Biofuels 4:55; del Rio et al. 2015, Biomass Bioenergy 81:322-338; release: Topakas et al. 2003, J Biotechnol 102:33-44; Topakas et al. 2003, Enzyme Microb Technol 33:729-737",
  "Main sugarcane cell-wall hydroxycinnamate, ester-linked to lignin and hemicellulose; Fusarium (F. oxysporum) feruloyl esterases hydrolyse p-coumarate esters; consistent with the increase in infected larvae",
  "sugarcane (diet)", "diet", "ferulic acid", "cell-wall hydroxycinnamate (released)",
  "Masarin et al. 2011, Biotechnol Biofuels 4:55; release: Topakas et al. 2003, J Biotechnol 102:33-44; Topakas et al. 2003, Enzyme Microb Technol 33:729-737",
  "Minor sugarcane cell-wall hydroxycinnamate, ester-linked to arabinoxylan; released by Fusarium (F. oxysporum) feruloyl esterases acting with xylanases",
  "sugarcane (diet)", "diet", "sinapic acid", "sugarcane phenolic acid",
  "Duarte-Almeida et al. 2006, Plant Foods Hum Nutr 61:187-192",
  "Reported in sugarcane juice (soluble phenolic); no report as a sugarcane cell-wall component",
  "sugarcane (diet)", "diet", "7,3',4'-Trihydroxyflavone", "plant flavone (class-level)",
  "class-level only: sugarcane is rich in flavones (Duarte-Almeida et al. 2006, Plant Foods Hum Nutr 61:187-192); this compound has not been reported in sugarcane",
  "putatively diet-derived; decreased in infected larvae",
  "sugarcane (diet)", "diet", "soyasaponin Ba", "triterpenoid saponin (identity uncertain)",
  "none: soyasaponins are legume compounds, not reported in sugarcane; annotation doubtful",
  "diet-derived, identity uncertain")
lit$evidence <- "literature-based"
write.csv(lit, file.path(out, "literature_links.csv"), row.names = FALSE)

## =============================================================== (a) metabolite evidence table
fly_paths <- mg %>% mutate(map = mapnum(kegg_id))
evid <- lapply(seq_len(nrow(met)), function(i) {
  m <- met[i, ]; k <- m$KEGG
  maps <- cp$map[cp$KEGG == k]
  hp <- fly_paths %>% filter(map %in% mapnum(maps))
  host_paths <- if (nrow(hp)) paste(sprintf("%s (NES %+.2f, padj %.2g)", hp$pathway, hp$transcriptome_NES,
                                            hp$transcriptome_padj), collapse = "; ") else ""
  rg <- rx_genes %>% filter(KEGG == k)
  rg_u <- rg %>% distinct(gene, .keep_all = TRUE)
  ma_k <- kkn %>% filter(kegg == k) %>% inner_join(dif_genera %>% dplyr::select(genus), by = "genus")
  ma_g <- gem %>% filter(kegg == k, passes_server_threshold) %>% inner_join(dif_genera %>% dplyr::select(genus), by = "genus")
  sk <- sank_dif %>% filter(KEGG == k)
  bact <- unique(c(ma_k$genus, ma_g$genus, sk$genus))
  bd <- dif_genera %>% filter(genus %in% bact)
  fmt <- function(d) if (nrow(d)) paste(sprintf("%s(%s)", d$genus, ifelse(d$genus_direction == "up_in_infected", "+", "-")),
                                        collapse = ", ") else ""
  data.frame(
    metabolite = m$metabolite, KEGG = k, log2FC_infected_vs_control = m$log2FC,
    detected_only_in_infected = m$only_infected, direction = m$direction, differential = m$differential,
    host_pathways_dmel_KEGG = host_paths,
    n_tested_genes_linked_by_KEGG_reaction = n_distinct(rg$gene),
    n_DE_genes_linked_by_KEGG_reaction = n_distinct(rg$gene[rg$DE]),
    genes_linked_by_KEGG_reaction = if (nrow(rg_u)) paste(sprintf("%s[%s|%s, log2FC %+.2f, padj %.2g]", rg_u$gene,
      coalesce(rg_u$preferred_name, rg_u$KO), coalesce(rg_u$dmel_symbol, "-"), rg_u$log2FC_infected_vs_control, rg_u$padj),
      collapse = "; ") else "",
    n_DA_genera_linked = nrow(bd),
    n_DA_genera_linked_up = sum(bd$genus_direction == "up_in_infected"),
    n_DA_genera_linked_down = sum(bd$genus_direction == "down_in_infected"),
    DA_genera_tier12_linked = fmt(bd %>% filter(tier12) %>% arrange(genus)),
    link_sources = paste(c(if (nrow(ma_k)) "MicrobiomeAnalyst KEGG", if (nrow(ma_g)) "MicrobiomeAnalyst GEM",
                           if (nrow(sk)) "MetOrigin KEGG reaction"), collapse = "; "),
    MetOrigin_origin_fly = m$MetOrigin_origin_fly %||% NA, MetOrigin_origin_human = m$MetOrigin_origin_human %||% NA)
})
evid <- bind_rows(evid)
n_tier12 <- sum(da$tier12)
evid <- evid %>% left_join(bind_rows(
  sank_dif %>% filter(tier12) %>% distinct(KEGG, genus),
  kkn %>% inner_join(dif_genera %>% filter(tier12), by = "genus") %>% distinct(KEGG = kegg, genus),
  gem %>% filter(passes_server_threshold) %>% inner_join(dif_genera %>% filter(tier12), by = "genus") %>%
    distinct(KEGG = kegg, genus)) %>% distinct() %>% count(KEGG, name = "n_tier12_genera_linked"), by = "KEGG") %>%
  mutate(n_tier12_genera_linked = coalesce(n_tier12_genera_linked, 0L),
         fraction_tier12_genera_linked = round(n_tier12_genera_linked / n_tier12, 3),
         microbial_link_specificity = ifelse(n_tier12_genera_linked == 0, "no link",
           ifelse(fraction_tier12_genera_linked > 0.25, "non-specific (ubiquitous metabolite)", "specific")))

# final class: database evidence first, literature for metabolites outside the databases
host_ev <- grepl("Host", evid$MetOrigin_origin_fly) | evid$n_tested_genes_linked_by_KEGG_reaction > 0
mic_ev  <- grepl("Microbiota", evid$MetOrigin_origin_fly) | evid$n_DA_genera_linked > 0
evid$host_evidence <- host_ev; evid$microbial_evidence <- mic_ev
evid$final_class <- ifelse(host_ev & mic_ev, "co-metabolism", ifelse(host_ev, "host",
                    ifelse(mic_ev, "microbial", "unassigned")))
evid$class_basis <- "database (MetOrigin / KEGG reaction / MicrobiomeAnalyst)"
lit_class <- c("fumonisin B1" = "external: fungus (Fusarium)", "Fusaric acid" = "external: fungus (Fusarium)",
               "desferrioxamine E" = "microbial (bacterial siderophore)", "desferrioxamine D2" = "microbial (bacterial siderophore)",
               "Gln-C18:3" = "co-metabolism", "Gln-C16:3" = "co-metabolism", "Gln-C17:3" = "co-metabolism",
               "p-coumaric acid" = "external: diet (sugarcane)", "ferulic acid" = "external: diet (sugarcane)",
               "sinapic acid" = "external: diet (sugarcane)", "7,3',4'-Trihydroxyflavone" = "external: diet (sugarcane; putative)",
               "soyasaponin Ba" = "external: diet (sugarcane; identity uncertain)")
evid$database_class <- evid$final_class
idx <- evid$metabolite %in% names(lit_class)
evid$final_class[idx] <- lit_class[evid$metabolite[idx]]
evid$class_basis[idx] <- ifelse(evid$database_class[idx] == "unassigned", "literature-based",
                                "literature-based (database class kept in database_class)")
evid$secondary_contribution <- ""
hc <- evid$metabolite %in% c("p-coumaric acid", "ferulic acid")
evid$secondary_contribution[hc] <- "microbial / fungal release (MetOrigin: Microbiota; literature: Fusarium feruloyl esterases)"
evid$secondary_contribution[evid$metabolite == "sinapic acid"] <- "microbial (MetOrigin: Microbiota)"
evid <- evid %>% left_join(lit %>% group_by(metabolite = target) %>%
                             summarise(literature = paste(sprintf("%s %s [%s]", source, relation, reference), collapse = " | ")),
                           by = "metabolite")
evid$literature[is.na(evid$literature)] <- ""
write.csv(evid, file.path(out, "4a_metabolite_evidence_table.csv"), row.names = FALSE)

## =============================================================== (b) pathway x omics matrix
dif_kegg <- met$KEGG[met$differential & met$KEGG != ""]
cp_d <- cp %>% filter(KEGG %in% dif_kegg) %>% mutate(map = mapnum(map), sign = met$sign[match(KEGG, met$KEGG)])
maen$map <- mapnum(pn$map[match(maen$Pathway, pn$pathway_name)])
pathways <- unique(c(fly_paths$map[fly_paths$transcriptome_padj < 0.05], cp_d$map, mpea$map, maen$map[maen$Pval < 0.05]))
pm <- data.frame(map = pathways) %>%
  mutate(pathway = coalesce(pn$pathway_name[match(paste0("map", map), pn$map)],
                            fly_paths$pathway[match(map, fly_paths$map)]))
tx <- fly_paths %>% dplyr::select(map, genes_ranked, transcriptome_NES, transcriptome_padj, layers_concordant,
                                  metabolome_NES, metabolome_pval, combined_padj_genes_metabolites = combined_padj)
pm <- pm %>% left_join(tx, by = "map")
mt <- cp_d %>% group_by(map) %>% summarise(n_diff_metabolites = n_distinct(KEGG),
  diff_metabolites = paste(sort(unique(key_kegg[KEGG])), collapse = "; "),
  metabolite_direction_score = mean(sign), .groups = "drop")
pm <- pm %>% left_join(mt, by = "map")
mp <- mpea %>% group_by(map) %>% summarise(MetOrigin_min_p = min(p),
  MetOrigin_categories = paste(sprintf("%s(p=%.2g)", category, p), collapse = "; "),
  MetOrigin_micro_sig = paste(category[p < 0.05], collapse = ";"), .groups = "drop")
pm <- pm %>% left_join(mp, by = "map") %>%
  left_join(maen %>% filter(!is.na(map)) %>% group_by(map) %>% summarise(MA_enrich_p = min(Pval)), by = "map")
# microbial links inside each pathway: reaction-level (MetOrigin) links of DA genera to differential metabolites
mlink <- sank_dif %>% filter(KEGG %in% dif_kegg) %>% group_by(map) %>%
  summarise(n_DA_genera_reaction_linked = n_distinct(genus),
            n_DA_tier12_reaction_linked = n_distinct(genus[tier12]),
            microbe_direction_score = mean(sgn(genus_direction[!duplicated(genus)])), .groups = "drop")
pm <- pm %>% left_join(mlink, by = "map") %>%
  mutate(across(c(n_diff_metabolites, n_DA_genera_reaction_linked, n_DA_tier12_reaction_linked), ~coalesce(.x, 0L)))
# support rules
pm <- pm %>% mutate(
  support_transcriptome = !is.na(transcriptome_padj) & transcriptome_padj < 0.05,
  support_metabolome = n_diff_metabolites > 0 & ((!is.na(metabolome_pval) & metabolome_pval < 0.05) |
                                                 (!is.na(MetOrigin_min_p) & MetOrigin_min_p < 0.05)),
  # a reaction link alone is not enough (ubiquitous metabolites link to most genera);
  # the pathway must also be enriched for microbe-associated metabolites
  support_microbiome = n_DA_tier12_reaction_linked > 0 &
    ((!is.na(MA_enrich_p) & MA_enrich_p < 0.05) | grepl("Microbiota|Co-Metabolism", MetOrigin_micro_sig)),
  n_omics_support = support_transcriptome + support_metabolome + support_microbiome,
  # sensitivity only: GSEA convention FDR < 0.25 for the transcriptome
  support_transcriptome_FDR0.25 = !is.na(transcriptome_padj) & transcriptome_padj < 0.25,
  n_omics_support_FDR0.25 = support_transcriptome_FDR0.25 + support_metabolome + support_microbiome) %>%
  filter(n_omics_support > 0 | n_omics_support_FDR0.25 > 0) %>% arrange(desc(n_omics_support), transcriptome_padj)
write.csv(pm, file.path(out, "4b_pathway_by_omics_matrix.csv"), row.names = FALSE)

# heatmap: direction (colour) and support (asterisk); descriptive only
hm <- as.matrix(data.frame(Genes = pm$transcriptome_NES / 2.5, Metabolites = pm$metabolite_direction_score,
                           Microbes = pm$microbe_direction_score))
hm <- pmax(pmin(hm, 1), -1); rownames(hm) <- pm$pathway
sup <- cbind(pm$support_transcriptome, pm$support_metabolome, pm$support_microbiome)
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2166AC", "#F7F7F7", "#B2182B"))
ht <- Heatmap(hm, name = "direction\n(infected vs\ncontrol)", col = col_fun, na_col = "grey92",
  cluster_rows = FALSE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 8),
  column_names_rot = 0, column_names_centered = TRUE,
  cell_fun = function(j, i, x, y, w, h, fill) if (isTRUE(sup[i, j])) grid.text("*", x, y, gp = gpar(fontsize = 12)),
  right_annotation = rowAnnotation(`omics\nsupport` = anno_barplot(pm$n_omics_support, width = unit(1.5, "cm"))),
  column_title = "Pathway support by omic\n(descriptive; * = support rule met)", column_title_gp = gpar(fontsize = 10),
  row_names_max_width = max_text_width(rownames(hm), gp = gpar(fontsize = 8)),
  heatmap_legend_param = list(at = c(-1, 0, 1), labels = c("down", "0", "up")))
for (ext in c("png", "pdf")) {
  h <- max(5, 0.2 * nrow(hm) + 2)
  if (ext == "png") png(file.path(out, "4b_pathway_omics_heatmap.png"), width = 9.5, height = h + 1, units = "in", res = 300)
  else pdf(file.path(out, "4b_pathway_omics_heatmap.pdf"), width = 9.5, height = h + 1)
  draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom"); dev.off()
}
ul <- list(Genes = pm$pathway[pm$support_transcriptome], Metabolites = pm$pathway[pm$support_metabolome],
           Microbes = pm$pathway[pm$support_microbiome])
for (ext in c("png", "pdf")) {
  if (ext == "png") png(file.path(out, "4b_pathway_support_upset.png"), width = 6.5, height = 4.5, units = "in", res = 300)
  else pdf(file.path(out, "4b_pathway_support_upset.pdf"), width = 6.5, height = 4.5)
  print(upset(fromList(ul), order.by = "freq", sets.bar.color = "grey40", mainbar.y.label = "Pathways",
              sets.x.label = "Pathways supported", text.scale = 1.2))
  dev.off()
}

## =============================================================== (c) direction agreement along chains
gm <- rx_genes %>% filter(KEGG %in% dif_kegg) %>%
  distinct(KEGG, gene, DE, gene_direction, gene_log2FC = log2FC_infected_vs_control, gene_role = compound_side_as_written, reaction)
tm <- bind_rows(
  sank_dif %>% filter(KEGG %in% dif_kegg) %>% distinct(KEGG, genus, genus_direction, tier12, taxon_role = role,
                                                       reaction_type = ReactionType, link = "MetOrigin KEGG reaction"),
  kkn %>% filter(kegg %in% dif_kegg) %>% inner_join(dif_genera, by = c("genus", "genus_direction", "confidence_tier")) %>%
    distinct(KEGG = kegg, genus, genus_direction, tier12, link = "MicrobiomeAnalyst KEGG pathway"),
  gem %>% filter(kegg %in% dif_kegg, passes_server_threshold) %>% inner_join(dif_genera %>% dplyr::select(genus, tier12), by = "genus") %>%
    distinct(KEGG = kegg, genus, genus_direction, tier12, link = "MicrobiomeAnalyst GEM"))
chains <- met %>% filter(differential, KEGG != "") %>% dplyr::select(metabolite, KEGG, met_direction = direction) %>%
  left_join(gm, by = "KEGG", relationship = "many-to-many") %>%
  left_join(tm %>% filter(tier12), by = "KEGG", relationship = "many-to-many") %>%
  mutate(gene_vs_metabolite = ifelse(is.na(gene), NA, ifelse(gene_direction == met_direction, "concordant", "discordant")),
         taxon_vs_metabolite = ifelse(is.na(genus), NA, ifelse(genus_direction == met_direction, "concordant", "discordant")),
         # role-adjusted: only for irreversible MetOrigin reactions; a producer should change
         # in the same direction as the metabolite, a consumer in the opposite direction
         taxon_vs_metabolite_role_adjusted = ifelse(is.na(genus), NA,
           ifelse(is.na(reaction_type) | reaction_type != "irreversible", "role ambiguous",
             ifelse((taxon_role == "product") == (genus_direction == met_direction), "concordant", "discordant"))),
         chain = ifelse(is.na(gene) | is.na(genus), "incomplete chain",
                 ifelse(gene_vs_metabolite == "concordant" & taxon_vs_metabolite == "concordant", "concordant",
                        "discordant")))
write.csv(chains, file.path(out, "4c_direction_chains.csv"), row.names = FALSE)
csum <- bind_rows(
  chains %>% filter(!is.na(gene)) %>% distinct(metabolite, gene, gene_vs_metabolite, DE) %>%
    count(metabolite, link_type = "gene-metabolite (KEGG reaction)", result = gene_vs_metabolite),
  chains %>% filter(!is.na(genus)) %>% distinct(metabolite, genus, taxon_vs_metabolite) %>%
    count(metabolite, link_type = "taxon-metabolite (tier 1-2)", result = taxon_vs_metabolite),
  chains %>% filter(chain != "incomplete chain") %>% distinct(metabolite, gene, genus, chain) %>%
    count(metabolite, link_type = "gene-metabolite-taxon chain", result = chain)) %>%
  pivot_wider(names_from = result, values_from = n, values_fill = 0)
write.csv(csum, file.path(out, "4c_direction_summary.csv"), row.names = FALSE)
lit_dir <- lit %>% left_join(da %>% dplyr::select(source = genus, source_direction = direction, confidence_tier), by = "source") %>%
  left_join(met %>% dplyr::select(target = metabolite, target_direction = direction, differential), by = "target") %>%
  mutate(direction_agreement = ifelse(is.na(source_direction), "source not a DA taxon",
                                  ifelse(source_direction == target_direction, "concordant", "discordant")))
write.csv(lit_dir, file.path(out, "4e_literature_links_direction.csv"), row.names = FALSE)

## =============================================================== (d) Cytoscape tables + figure
cap <- max(abs(met$log2FC), na.rm = TRUE) + 1
nodes_met <- met %>% transmute(id = metabolite, type = "metabolite", label = metabolite,
  log2FC = ifelse(only_infected, cap, log2FC), log2FC_note = ifelse(only_infected, "detected only in infected (capped)", ""),
  differential, direction, final_class = evid$final_class[match(metabolite, evid$metabolite)])
e_gene <- rx_genes %>% filter(KEGG != "") %>%
  distinct(source = gene, target = key_kegg[KEGG], reaction, KO, role = compound_side_as_written) %>%
  group_by(source, target) %>% summarise(edge_type = "KEGG reaction (host gene)", detail = paste(unique(paste(reaction, KO, role)), collapse = "; "),
                                         .groups = "drop")
e_tax <- bind_rows(
  sank_dif %>% filter(tier12) %>% group_by(source = genus, target = metabolite) %>%
    summarise(edge_type = "KEGG reaction (MetOrigin)", detail = paste(unique(paste(reaction, role, ReactionType)), collapse = "; "), .groups = "drop"),
  kkn %>% inner_join(dif_genera %>% filter(tier12) %>% dplyr::select(genus), by = "genus") %>%
    transmute(source = genus, target = metabolite_study, edge_type = "MicrobiomeAnalyst KEGG pathway", detail = pathways),
  gem %>% filter(passes_server_threshold) %>% inner_join(dif_genera %>% filter(tier12) %>% dplyr::select(genus), by = "genus") %>%
    transmute(source = genus, target = metabolite_study, edge_type = paste("MicrobiomeAnalyst GEM", db),
              detail = sprintf("potential %.2f", potential)))
e_lit <- lit %>% transmute(source, target, edge_type = "literature-based", detail = reference)
edges <- bind_rows(e_gene, e_tax, e_lit) %>% filter(target %in% met$metabolite)
nodes_gene <- genes %>% filter(gene %in% e_gene$source) %>% transmute(id = gene, type = "host gene",
  label = ifelse(!is.na(preferred_name) & preferred_name != "", preferred_name,
                 ifelse(!is.na(dmel_symbol) & dmel_symbol != "", paste0("Dm:", dmel_symbol), gene)),
  log2FC = log2FC_infected_vs_control, log2FC_note = "",
  differential = DE, direction, final_class = NA)
tax_ids <- unique(c(e_tax$source, e_lit$source[e_lit$source %in% da$genus]))
nodes_tax <- da %>% filter(genus %in% tax_ids) %>% transmute(id = genus, type = ifelse(kingdom == "Fungi", "fungus", "bacterium"),
  label = genus, log2FC = log2_effect_infected_vs_control, log2FC_note = "ALDEx2 log2 CLR difference",
  differential, direction, final_class = confidence_tier)
nodes_other <- data.frame(id = setdiff(unique(e_lit$source), c(nodes_tax$id)), type = "external source") %>%
  mutate(label = id, log2FC = NA, log2FC_note = "", differential = NA, direction = NA, final_class = NA)
nodes <- bind_rows(nodes_met, nodes_gene, nodes_tax, nodes_other)
edges <- edges %>% filter(source %in% nodes$id)
write.csv(nodes, file.path(out, "4d_cytoscape_nodes.csv"), row.names = FALSE)
write.csv(edges, file.path(out, "4d_cytoscape_edges.csv"), row.names = FALSE)

# static figure: reaction-level and literature edges (GEM / KEGG-pathway edges are in the table only)
ubiq <- evid$metabolite[evid$microbial_link_specificity == "non-specific (ubiquitous metabolite)"]
fe <- edges %>% filter(edge_type %in% c("KEGG reaction (host gene)", "KEGG reaction (MetOrigin)", "literature-based")) %>%
  group_by(source, target) %>% summarise(edge_type = first(edge_type), .groups = "drop")
tax_ids_all <- nodes$id[nodes$type %in% c("bacterium", "fungus")]
keep_tax <- unique(c(fe$source[fe$edge_type != "KEGG reaction (MetOrigin)"],          # literature taxa
                     fe$source[fe$edge_type == "KEGG reaction (MetOrigin)" & !fe$target %in% ubiq]))
agg <- fe %>% filter(edge_type == "KEGG reaction (MetOrigin)", target %in% ubiq) %>%
  left_join(da %>% dplyr::select(source = genus, direction), by = "source") %>%
  group_by(target) %>% summarise(source = sprintf("%s: %d tier 1-2 genera\n(%d up / %d down)", first(target), n(),
    sum(direction == "up_in_infected"), sum(direction == "down_in_infected")), .groups = "drop") %>%
  mutate(edge_type = "KEGG reaction (MetOrigin)")
fe <- bind_rows(fe %>% filter(!(edge_type == "KEGG reaction (MetOrigin)" & target %in% ubiq)), agg) %>%
  filter(!(source %in% tax_ids_all & !source %in% keep_tax))
fn <- bind_rows(nodes %>% filter(id %in% c(fe$source, fe$target)),
                data.frame(id = unique(agg$source), type = "genus group (ubiquitous metabolite)", label = unique(agg$source)))
g <- graph_from_data_frame(fe, vertices = fn, directed = FALSE)
V(g)$lfc <- pmax(pmin(V(g)$log2FC, 8), -8)
p <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(linetype = edge_type, colour = edge_type), alpha = 0.6, width = 0.4) +
  geom_node_point(aes(shape = type, fill = lfc, size = type == "metabolite"), colour = "grey20") +
  geom_node_text(aes(label = label), size = 2.2, repel = TRUE, max.overlaps = 40) +
  scale_shape_manual(values = c(metabolite = 22, `host gene` = 21, bacterium = 24, fungus = 23, `external source` = 25,
                                `genus group (ubiquitous metabolite)` = 8)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", na.value = "grey80", name = "log2FC / effect\n(capped at 8)") +
  scale_size_manual(values = c(`TRUE` = 4, `FALSE` = 2.5), guide = "none") +
  scale_edge_colour_manual(values = c(`KEGG reaction (host gene)` = "grey30", `KEGG reaction (MetOrigin)` = "#7B3294",
                                      `literature-based` = "#E08214"), name = "edge") +
  scale_edge_linetype_manual(values = c(`KEGG reaction (host gene)` = "solid", `KEGG reaction (MetOrigin)` = "solid",
                                        `literature-based` = "dashed"), name = "edge") +
  theme_graph(base_family = "sans") + labs(title = "Knowledge-based network (no correlation edges)",
       caption = paste0("Taxa: differential genera of confidence tiers 1-2; edges: KEGG reactions (host genes, MetOrigin) and literature.\n",
                        "Genera linked to ubiquitous metabolites (>25% of tier 1-2 genera) are collapsed into one node per metabolite; ",
                        "all edges are in 4d_cytoscape_edges.csv")) + theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(out, "4d_network.png"), p, width = 12, height = 10, dpi = 300)
ggsave(file.path(out, "4d_network.pdf"), p, width = 12, height = 10)

## =============================================================== report
options(width = 220)
cat("== (a) metabolite evidence\n")
print(evid %>% dplyr::select(metabolite, direction, n_tested_genes_linked_by_KEGG_reaction, n_DE_genes_linked_by_KEGG_reaction,
                             n_DA_genera_linked, n_tier12_genera_linked, microbial_link_specificity,
                             database_class, final_class), row.names = FALSE)
cat("\n== (b) pathways:", nrow(pm), "| supported by 1/2/3 omics:", paste(table(factor(pm$n_omics_support, 1:3)), collapse = "/"),
    "| sensitivity (transcriptome FDR<0.25):", paste(table(factor(pm$n_omics_support_FDR0.25, 1:3)), collapse = "/"), "\n")
print(pm %>% filter(n_omics_support >= 2 | n_omics_support_FDR0.25 >= 2) %>% dplyr::select(pathway, transcriptome_NES, transcriptome_padj, diff_metabolites,
                                                            MetOrigin_min_p, MA_enrich_p, n_DA_tier12_reaction_linked,
                                                            n_omics_support, n_omics_support_FDR0.25) %>% mutate(across(where(is.numeric), ~signif(.x, 3))), row.names = FALSE)
cat("\n== (c) direction summary\n"); print(csum, row.names = FALSE)
cat("\n== (e) literature links / positive control\n")
print(lit_dir %>% dplyr::select(source, target, source_direction, confidence_tier, target_direction, direction_agreement), row.names = FALSE)
cat("\n== (d) network: nodes", nrow(nodes), "(", paste(names(table(nodes$type)), table(nodes$type), collapse = ", "),
    ") | edges", nrow(edges), "(", paste(names(table(edges$edge_type)), table(edges$edge_type), collapse = ", "), ")\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
