#!/usr/bin/env Rscript
# Step 5 - Deliverables: one Excel workbook with all tables, the final figures, and the methods.
# Run from the repository root with the integration_diatrea env:  Rscript scripts/step5b_deliverables.R
suppressPackageStartupMessages({ library(openxlsx); library(dplyr) })
out <- "deliverables"; fig <- file.path(out, "figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
rc <- function(f, ...) read.csv(f, check.names = FALSE, ...)
rt <- function(f, ...) read.delim(f, check.names = FALSE, ...)
mo <- "web_results/metorigin"

sheets <- list(
  # ---- step 0
  "T01_DE_genes" = list(rc("results/step1_gene_annotation/de_genes_annotated.csv"),
    "474 DE genes (DESeq2, padj<0.05), log2FC sign flipped to infected vs control, with eggNOG KO / name / D. melanogaster ortholog (steps 0, 1)"),
  "T02_metabolome" = list(rc("results/step0/metabolome_summary_infected_vs_control.csv"),
    "18 metabolites: log2FC infected vs control (mean, median), detected-only-in-infected, published p / p_adj (step 0)"),
  "T03_metab_outliers" = list(rc("results/step0/metabolome_long_with_flags.csv") %>% filter(flag_outlier),
    "Values > 5x from their group median (reported, NOT changed) (step 0)"),
  "T04_micro_reads" = list(rt("results/step0a_microbiome_cascade/read_accounting_by_category.tsv"),
    "Microbiome read accounting per sample and category, cascade (primary) and strict variant (step 0a)"),
  "T05_micro_genus_table" = list(rt("results/step0a_microbiome_cascade/abundance_table_genus_cascade.tsv"),
    "Genus counts from the cascade classification SILVA 138.2 -> NCBI 16S/18S/28S/ITS -> PlusPFP-16 (step 0a)"),
  "T06_micro_vs_published" = list(rt("results/step0a_microbiome_cascade/comparison_with_published_tables.tsv"),
    "Cascade vs earlier tables (log2 ratio of mean relative abundance, descriptive) (step 0a)"),
  "T07_microbiome_DA" = list(rc("results/step0b_microbiome_DA/microbiome_DA.csv"),
    "Microbiome differential abundance: ALDEx2 (primary) + ANCOM-BC2 + strict cascade; rule and confidence tiers (step 0b)"),
  "T08_PERMANOVA" = list(rc("results/step0b_microbiome_DA/permanova_aitchison.csv"),
    "PERMANOVA on Aitchison distance, all 720 permutations; descriptive only (min attainable p = 0.1) (step 0b)"),
  "T09_annot_coverage" = list(rt("results/step1_gene_annotation/annotation_coverage_percent.tsv"),
    "Gene annotation coverage (% of genes) (step 1)"),
  "T10_metabolite_IDs" = list(rc("results/step2_metabolite_ids/metabolite_ids_final.csv"),
    "Metabolite identifiers (KEGG, HMDB, ChEBI, PubChem) with curation notes (step 2)"),
  "T11_metab_unmapped" = list(rc("results/step2_metabolite_ids/metabolite_manual_table_unmapped.csv"),
    "Manual table: metabolites without a KEGG compound (step 2)"),
  "T12_multiGSEA" = list(rc("results/step3_multigsea/multigsea_kegg_dmel.csv"),
    "multiGSEA: KEGG D. melanogaster pathways, transcriptome + metabolome, Stouffer (genes + metabolites only) (step 3)"),
  "T13_MetOrigin_origin" = list(full_join(
      rc(file.path(mo, "fruit_fly_host/02_Origin_Analysis/Metabolites_Origin.csv")) %>% rename(Origin_fruit_fly = Origin, Links_fruit_fly = Links),
      rc(file.path(mo, "human_host_default/02_Origin_Analysis/Metabolites_Origin.csv")) %>% dplyr::select(Name, Origin_human = Origin),
      by = "Name"), "MetOrigin (Simple MetOrigin Analysis) origin of each metabolite; fruit-fly host (primary) and human host"),
  "T14_MetOrigin_MPEA" = list(bind_rows(lapply(c("Host", "Microbiota", "Co-Metabolism"), function(k) {
      d <- rc(file.path(mo, "fruit_fly_host/03_Pathway_Analysis", paste0("MPEA_Results_", k, ".csv")), colClasses = "character")
      if (nrow(d)) d %>% mutate(category = k, .before = 1) })), "MetOrigin pathway enrichment (MPEA), fruit-fly host"),
  "T15_MA_GEM" = list(rc("results/step4a_microbiomeanalyst_local/gem_potential_long.csv"),
    "MicrobiomeAnalyst AGORA/EMBL genome-scale-model potential, genus x metabolite (reproduced locally from the server libraries) (step 4a)"),
  "T16_MA_KEGG_links" = list(rc("results/step4a_microbiomeanalyst_local/kegg_knowledge_links_genus_metabolite.csv"),
    "MicrobiomeAnalyst KEGG knowledge links: genus - metabolite via KEGG pathways of the genus' organisms (step 4a)"),
  "T17_MA_enrichment" = list(rc("results/step4a_microbiomeanalyst_local/network_enrichment_usrbac.csv"),
    "MicrobiomeAnalyst network enrichment (usrbac library) with pathway names restored (step 4a)"),
  "T18_4a_metabolites" = list(rc("results/step4_integration/4a_metabolite_evidence_table.csv"),
    "4a Metabolite-centric evidence table and final class (step 4)"),
  "T19_4b_pathways" = list(rc("results/step4_integration/4b_pathway_by_omics_matrix.csv"),
    "4b Pathways x omics matrix: direction and support per omic (descriptive; no combined p with the microbiome) (step 4)"),
  "T20_4c_direction" = list(rc("results/step4_integration/4c_direction_summary.csv"),
    "4c Direction agreement summary per metabolite (step 4)"),
  "T21_4c_chains" = list(rc("results/step4_integration/4c_direction_chains.csv"),
    "4c Gene - metabolite - taxon chains (tier 1-2 taxa), same-sign and role-adjusted direction agreement (step 4)"),
  "T22_4d_nodes" = list(rc("results/step4_integration/4d_cytoscape_nodes.csv"), "4d Cytoscape node table (step 4)"),
  "T23_4d_edges" = list(rc("results/step4_integration/4d_cytoscape_edges.csv"),
    "4d Cytoscape edge table: KEGG reactions, MicrobiomeAnalyst, literature; no correlation edges (step 4)"),
  "T24_4e_literature" = list(rc("results/step4_integration/4e_literature_links_direction.csv"),
    "4e Literature-based links and their direction agreement; positive control Fusarium -> fumonisin B1 / fusaric acid"),
  "T25_validation" = list(rc("results/validation_steps0-3.csv"), "Independent validation of steps 0-3 (scripts/step5a_validate_outputs.R)"))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom")
addWorksheet(wb, "README")
readme <- data.frame(sheet = names(sheets), description = sapply(sheets, `[[`, 2),
                     rows = sapply(sheets, function(x) nrow(x[[1]])))
writeData(wb, "README", data.frame(x = c(
  "Diatraea saccharalis - knowledge-based, summary-level integration of omics from independent cohorts",
  "Direction everywhere: infected (F. verticillioides-contaminated sugarcane) vs control. Methods and run instructions: README.md",
  "")), colNames = FALSE)
writeData(wb, "README", readme, startRow = 4, headerStyle = hs)
setColWidths(wb, "README", 1:3, c(24, 140, 8))
for (nm in names(sheets)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, sheets[[nm]][[2]], startRow = 1)
  writeData(wb, nm, sheets[[nm]][[1]], startRow = 3, headerStyle = hs, withFilter = TRUE)
  freezePane(wb, nm, firstActiveRow = 4)
  setColWidths(wb, nm, seq_len(ncol(sheets[[nm]][[1]])), "auto")
}
saveWorkbook(wb, file.path(out, "Diatraea_integration_tables.xlsx"), overwrite = TRUE)

figs <- c(
  "Fig_microbiome_DA_effect" = "results/step0b_microbiome_DA/microbiome_DA_effect_plot",
  "Fig_multiGSEA_pathways" = "results/step3_multigsea/multigsea_top_pathways",
  "Fig_pathway_omics_heatmap" = "results/step4_integration/4b_pathway_omics_heatmap",
  "Fig_pathway_support_upset" = "results/step4_integration/4b_pathway_support_upset",
  "Fig_knowledge_network" = "results/step4_integration/4d_network")
for (n in names(figs)) for (ext in c("png", "pdf"))
  file.copy(paste0(figs[[n]], ".", ext), file.path(fig, paste0(n, ".", ext)), overwrite = TRUE)
file.copy("results/step4_integration/4d_cytoscape_nodes.csv", file.path(out, "cytoscape_nodes.csv"), overwrite = TRUE)
file.copy("results/step4_integration/4d_cytoscape_edges.csv", file.path(out, "cytoscape_edges.csv"), overwrite = TRUE)
cat("Workbook:", length(sheets) + 1, "sheets;", sum(readme$rows), "rows | figures:", length(list.files(fig)), "\n")
print(readme[, c("sheet", "rows")], row.names = FALSE)
