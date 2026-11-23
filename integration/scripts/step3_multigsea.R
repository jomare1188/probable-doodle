#!/usr/bin/env Rscript
# Step 3 - Transcriptome + metabolome pathway analysis with multiGSEA (KEGG, D. melanogaster).
#   Transcriptome: published DESeq2 Wald statistic (sign flipped to infected vs control), mapped to
#     D. melanogaster orthologs (step 1); per fly gene, the Diatraea gene with the largest |stat|.
#   Metabolome: published per-metabolite p-values, rank = sign(infected - control) * -log10(p),
#     KEGG compound IDs from step 2.
#   Pathways: graphite KEGG pathways for dmelanogaster; gene nodes = FlyBase CG IDs, metabolite
#     nodes = KEGG compound IDs (no ID conversion needed). Global/overview maps are excluded.
#   Per-layer fgsea (multiGSEA), combined with Stouffer's method (multiGSEA::combinePvalues);
#     BH adjustment of each layer and of the combined p.
# Run from the repository root with the integration_diatrea env:  Rscript scripts/step3_multigsea.R
suppressPackageStartupMessages({
  library(multiGSEA); library(graphite); library(org.Dm.eg.db); library(dplyr); library(ggplot2)
})
set.seed(20260930)
out_dir <- "results/step3_multigsea"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------ pathways (graphite KEGG, D. melanogaster)
kegg <- pathways("dmelanogaster", "kegg")
retrieved <- format(kegg[[1]]@timestamp)
overview <- c("Metabolic pathways", "Biosynthesis of amino acids", "Biosynthesis of cofactors",
              "Carbon metabolism", "2-Oxocarboxylic acid metabolism", "Fatty acid metabolism",
              "Degradation of aromatic compounds", "Biosynthesis of nucleotide sugars")
ids <- sapply(kegg, function(p) p@id)
is_overview <- names(kegg) %in% overview | grepl("^dme:01[12]", ids)
cat("graphite", as.character(packageVersion("graphite")), "- KEGG dmelanogaster:", length(kegg),
    "pathways, retrieved on", retrieved, "| overview maps excluded:", sum(is_overview), "\n")
kegg <- kegg[!is_overview]
strip <- function(x) sub("^[A-Z]+:", "", x)
sets <- list(
  transcriptome = lapply(as.list(kegg), function(p) unique(strip(nodes(p, "proteins")))),
  metabolome    = lapply(as.list(kegg), function(p) unique(strip(nodes(p, "metabolites")))))
names(sets$transcriptome) <- names(sets$metabolome) <- names(kegg)

## ------------------------------------------------ transcriptome ranks
deseq <- read.csv("results/step0/deseq2_all_genes_infected_vs_control.csv")
orth  <- read.delim("results/step1_gene_annotation/gene2dmel_orthologs_long.tsv", colClasses = "character")
cg <- suppressMessages(AnnotationDbi::select(org.Dm.eg.db, keys = unique(orth$dmel_entrez),
                                              keytype = "ENTREZID", columns = "FLYBASECG"))
tx <- orth %>% inner_join(cg, by = c("dmel_entrez" = "ENTREZID"), relationship = "many-to-many") %>%
  filter(!is.na(FLYBASECG)) %>%
  inner_join(deseq %>% dplyr::select(gene = gene_id, stat = stat_infected_vs_control), by = "gene") %>%
  group_by(FLYBASECG) %>% slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
write.csv(tx, file.path(out_dir, "transcriptome_rank_dmel_CG.csv"), row.names = FALSE)
rank_tx <- setNames(tx$stat, tx$FLYBASECG)
cat("Transcriptome: Diatraea genes with a fly ortholog in DESeq2:", length(unique(tx$gene)),
    "| fly genes ranked:", length(rank_tx), "\n")

## ------------------------------------------------ metabolome ranks
mids <- read.csv("results/step2_metabolite_ids/metabolite_ids_final.csv")
pub  <- read.csv("results/step0/metabolome_summary_infected_vs_control.csv")
met  <- mids %>% filter(KEGG != "") %>% dplyr::select(metabolite, KEGG) %>%
  inner_join(pub %>% dplyr::select(metabolite, published_p, direction), by = "metabolite") %>%
  mutate(rank = ifelse(direction == "up_in_infected", 1, -1) * -log10(published_p))
write.csv(met, file.path(out_dir, "metabolome_rank_kegg.csv"), row.names = FALSE)
rank_met <- setNames(met$rank, met$KEGG)
cat("Metabolome: KEGG-mapped metabolites ranked:", length(rank_met), "\n")

## ------------------------------------------------ coverage of each pathway
cov <- data.frame(pathway = names(kegg), kegg_id = sapply(kegg, function(p) p@id),
  genes_in_pathway = lengths(sets$transcriptome),
  genes_ranked = sapply(sets$transcriptome, function(s) sum(s %in% names(rank_tx))),
  metabolites_in_pathway = lengths(sets$metabolome),
  metabolites_measured = sapply(sets$metabolome, function(s) sum(s %in% names(rank_met))),
  measured_metabolites = sapply(sets$metabolome, function(s)
    paste(met$metabolite[match(intersect(s, names(rank_met)), met$KEGG)], collapse = "; ")))

## ------------------------------------------------ multiGSEA (fgsea per layer) + Stouffer
ranks <- list(transcriptome = rank_tx, metabolome = rank_met)
# keep only sets with >= 1 measured feature in the layer (fgsea drops the rest anyway)
es <- multiGSEA(sets, ranks, eps = 0)
res <- extractPvalues(enrichmentScores = es, pathwayNames = names(kegg))
res$combined_pval <- combinePvalues(res, method = "stouffer")
res$combined_padj <- p.adjust(res$combined_pval, method = "BH")
nes <- function(layer) { d <- es[[layer]]; setNames(d$NES, d$pathway)[names(kegg)] }
size <- function(layer) { d <- es[[layer]]; setNames(d$size, d$pathway)[names(kegg)] }
lead <- function(layer) { d <- es[[layer]]
  setNames(sapply(d$leadingEdge, paste, collapse = ","), d$pathway)[names(kegg)] }
out <- cov %>% mutate(
  transcriptome_NES = nes("transcriptome"), transcriptome_size = size("transcriptome"),
  transcriptome_pval = res$transcriptome_pval, transcriptome_padj = res$transcriptome_padj,
  metabolome_NES = nes("metabolome"), metabolome_size = size("metabolome"),
  metabolome_pval = res$metabolome_pval, metabolome_padj = res$metabolome_padj,
  combined_pval_stouffer = res$combined_pval, combined_padj = res$combined_padj,
  transcriptome_leading_edge_CG = lead("transcriptome"),
  n_layers_tested = (!is.na(transcriptome_pval)) + (!is.na(metabolome_pval)),
  transcriptome_direction = ifelse(is.na(transcriptome_NES), NA,
                                   ifelse(transcriptome_NES > 0, "up_in_infected", "down_in_infected")),
  metabolome_direction = ifelse(is.na(metabolome_NES), NA,
                                ifelse(metabolome_NES > 0, "up_in_infected", "down_in_infected")),
  # Stouffer on fgsea p-values ignores direction: flag pathways where the two layers disagree
  layers_concordant = ifelse(n_layers_tested < 2, NA, sign(transcriptome_NES) == sign(metabolome_NES))) %>%
  arrange(combined_pval_stouffer)
# leading-edge fly genes -> Diatraea genes
cg2dia <- tapply(tx$gene, tx$FLYBASECG, paste, collapse = ",")
out$transcriptome_leading_edge_Diatraea <- sapply(strsplit(out$transcriptome_leading_edge_CG, ","),
  function(v) if (all(is.na(v))) NA else paste(unique(cg2dia[v]), collapse = ","))
write.csv(out, file.path(out_dir, "multigsea_kegg_dmel.csv"), row.names = FALSE)

## ------------------------------------------------ report
options(width = 220)
cat("\nPathways tested - transcriptome:", sum(!is.na(out$transcriptome_pval)),
    "| metabolome:", sum(!is.na(out$metabolome_pval)), "| both:", sum(out$n_layers_tested == 2), "\n")
cat("transcriptome padj<0.05:", sum(out$transcriptome_padj < 0.05, na.rm = TRUE),
    "| metabolome padj<0.05:", sum(out$metabolome_padj < 0.05, na.rm = TRUE),
    "| combined (Stouffer) padj<0.05:", sum(out$combined_padj < 0.05, na.rm = TRUE), "\n")
cat("\nPathways with measured metabolites:\n")
print(out %>% filter(n_layers_tested == 2) %>%
        dplyr::select(pathway, genes_ranked, transcriptome_NES, transcriptome_padj, measured_metabolites,
                      metabolome_NES, metabolome_pval, combined_pval_stouffer, combined_padj) %>%
        mutate(across(where(is.numeric), ~signif(.x, 3))), row.names = FALSE)
cat("\nTop 25 by combined p:\n")
print(head(out %>% dplyr::select(pathway, genes_ranked, transcriptome_NES, transcriptome_padj,
                                 metabolites_measured, metabolome_NES, combined_padj) %>%
             mutate(across(where(is.numeric), ~signif(.x, 3))), 25), row.names = FALSE)

## ------------------------------------------------ figure: significant pathways
sig <- out %>% filter(combined_padj < 0.05 | transcriptome_padj < 0.05) %>% head(30)
if (nrow(sig)) {
  sig$pathway <- factor(sig$pathway, levels = rev(sig$pathway))
  p <- ggplot(sig, aes(transcriptome_NES, pathway, size = genes_ranked, colour = -log10(combined_padj))) +
    geom_vline(xintercept = 0, colour = "grey60") + geom_point() +
    geom_point(data = subset(sig, n_layers_tested == 2), shape = 21, colour = "black", stroke = 1) +
    scale_colour_viridis_c(name = "-log10 combined padj") +
    labs(x = "Transcriptome NES (infected vs control)", y = NULL, size = "Fly genes ranked",
         caption = "Black ring: pathway also contains measured metabolites (combined with Stouffer)") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_dir, "multigsea_top_pathways.png"), p, width = 8.5, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "multigsea_top_pathways.pdf"), p, width = 8.5, height = 7)
}
writeLines(c(paste("graphite KEGG dmelanogaster retrieved on", retrieved),
             capture.output(sessionInfo())), file.path(out_dir, "sessionInfo.txt"))
