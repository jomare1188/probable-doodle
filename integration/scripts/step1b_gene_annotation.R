#!/usr/bin/env Rscript
# Step 1b - Gene-level functional annotation of Diatraea saccharalis genes (DIATSA_LOCUS*).
#   KEGG KO, EC, KEGG reactions/pathways, gene name: published eggNOG-mapper annotation
#   D. melanogaster orthologs: eggNOG orthologs of the same seed hits (step 1a), FlyBase protein
#   IDs mapped to FlyBase gene / Entrez / symbol with org.Dm.eg.db
# Proteins are collapsed to genes with gene2protein.txt (union over isoforms).
# Run from the repository root with the integration_diatrea env:
#   Rscript scripts/step1b_gene_annotation.R
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(org.Dm.eg.db) })

out_dir <- "results/step1_gene_annotation"
# the header line of this file is "query ..." (no leading #); only the trailer lines start with ##
annot   <- read.delim("data/annotation/proteins.emapper.emapper.annotations.gz",
                      comment.char = "#", header = TRUE, quote = "")
g2p  <- read.table("data/annotation/gene2protein.txt", col.names = c("gene", "protein"))
# the .orthologs file starts with "##" run-information lines, then the "#query" header
# (fread stops, as intended, at the "##" trailer lines at the end of the file)
orth_gz <- file.path(out_dir, "dmel_orthologs.emapper.orthologs.gz")
orth <- as.data.frame(suppressWarnings(data.table::fread(cmd = paste("gzip -dc", shQuote(orth_gz)),
                                        skip = "#query", sep = "\t", quote = "")))
names(orth) <- c("query", "orth_type", "species", "orthologs")
deseq <- read.csv("results/step0/deseq2_all_genes_infected_vs_control.csv")
de    <- read.csv("results/step0/de_genes_infected_vs_control.csv")

na_dash <- function(x) ifelse(x %in% c("-", ""), NA, x)
collapse_set <- function(x) {
  v <- unique(unlist(strsplit(na.omit(x), ",")))
  if (length(v)) paste(sort(v), collapse = ",") else NA_character_
}

## ------------------------------------------------ eggNOG functional annotation, per gene
fa <- annot %>%
  transmute(protein = query, preferred_name = na_dash(Preferred_name), description = na_dash(Description),
            EC = na_dash(EC), KEGG_ko = na_dash(gsub("ko:", "", KEGG_ko)),
            KEGG_pathway = na_dash(KEGG_Pathway), KEGG_reaction = na_dash(KEGG_Reaction),
            COG_category = na_dash(COG_category), max_annot_lvl, evalue) %>%
  inner_join(g2p, by = "protein")
gene_fa <- fa %>% group_by(gene) %>% summarise(
  proteins_annotated = paste(protein, collapse = ","),
  preferred_name = collapse_set(preferred_name),
  description = description[which.min(evalue)],
  EC = collapse_set(EC), KEGG_ko = collapse_set(KEGG_ko),
  KEGG_pathway = collapse_set(grep("^map", unlist(strsplit(na.omit(KEGG_pathway), ",")), value = TRUE)),
  KEGG_reaction = collapse_set(KEGG_reaction), COG_category = collapse_set(COG_category),
  eggnog_max_annot_lvl = max_annot_lvl[which.min(evalue)], .groups = "drop")

## ------------------------------------------------ D. melanogaster orthologs, per gene
dm <- orth %>% filter(species == "Drosophila melanogaster(7227)") %>%
  separate_rows(orthologs, sep = ",") %>%
  mutate(FBpp = sub("^\\*", "", orthologs), used_for_annotation = grepl("^\\*", orthologs)) %>%
  inner_join(g2p, by = c("query" = "protein"))
fb <- suppressMessages(AnnotationDbi::select(org.Dm.eg.db, keys = unique(dm$FBpp), keytype = "FLYBASEPROT",
                                              columns = c("FLYBASE", "ENTREZID", "SYMBOL")))
cat("FlyBase proteins:", length(unique(dm$FBpp)), "| mapped to Entrez:",
    sum(!is.na(fb$ENTREZID[!duplicated(fb$FLYBASEPROT)])), "\n")
dm <- dm %>% left_join(fb, by = c("FBpp" = "FLYBASEPROT"))
# relationship per gene: the "tightest" one among its isoforms
rank_type <- c(one2one = 1, many2one = 2, one2many = 3, many2many = 4)
gene_dm <- dm %>% group_by(gene) %>% summarise(
  dmel_orth_type = names(rank_type)[min(rank_type[orth_type])],
  dmel_FBgn = collapse_set(FLYBASE), dmel_entrez = collapse_set(ENTREZID),
  dmel_symbol = collapse_set(SYMBOL), .groups = "drop")
# long table gene -> Drosophila gene (Entrez), for step 3
dm_long <- dm %>% filter(!is.na(ENTREZID)) %>%
  distinct(gene, orth_type, dmel_FBgn = FLYBASE, dmel_entrez = ENTREZID, dmel_symbol = SYMBOL)
write.table(dm_long, file.path(out_dir, "gene2dmel_orthologs_long.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

## ------------------------------------------------ gene table
genes <- data.frame(gene = unique(c(g2p$gene, deseq$gene_id))) %>%
  left_join(gene_fa, by = "gene") %>% left_join(gene_dm, by = "gene") %>%
  mutate(in_deseq2 = gene %in% deseq$gene_id, DE = gene %in% de$gene_id)
genes <- genes %>% left_join(deseq %>% dplyr::select(gene = gene_id, log2FC_infected_vs_control,
                                              stat_infected_vs_control, padj, direction), by = "gene")
write.table(genes, file.path(out_dir, "gene_annotation_all_genes.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE, na = "")
de_annot <- genes %>% filter(DE) %>% arrange(desc(abs(log2FC_infected_vs_control)))
write.csv(de_annot, file.path(out_dir, "de_genes_annotated.csv"), row.names = FALSE, na = "")

## ------------------------------------------------ coverage report
cover <- function(d, label) {
  n <- nrow(d)
  data.frame(set = label, n_genes = n,
             eggNOG_hit = sum(!is.na(d$proteins_annotated)),
             preferred_name = sum(!is.na(d$preferred_name)),
             KEGG_KO = sum(!is.na(d$KEGG_ko)),
             KEGG_reaction = sum(!is.na(d$KEGG_reaction)),
             Dmel_ortholog = sum(!is.na(d$dmel_entrez)),
             Dmel_one2one = sum(d$dmel_orth_type %in% "one2one" & !is.na(d$dmel_entrez)),
             KO_or_Dmel = sum(!is.na(d$KEGG_ko) | !is.na(d$dmel_entrez)))
}
cv <- rbind(cover(genes %>% filter(in_deseq2), "all genes tested by DESeq2"),
            cover(de_annot, "DE genes (padj<0.05)"),
            cover(de_annot %>% filter(direction == "up_in_infected"), "DE up in infected"),
            cover(de_annot %>% filter(direction == "down_in_infected"), "DE down in infected"))
pct <- cv; pct[, -(1:2)] <- round(100 * cv[, -(1:2)] / cv$n_genes, 1)
write.table(cv, file.path(out_dir, "annotation_coverage_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pct, file.path(out_dir, "annotation_coverage_percent.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
options(width = 200)
cat("\n== Annotation coverage (counts)\n"); print(cv, row.names = FALSE)
cat("\n== Annotation coverage (% of genes in the set)\n"); print(pct, row.names = FALSE)
cat("\nDrosophila ortholog relationship among DE genes with an ortholog:\n")
print(table(de_annot$dmel_orth_type[!is.na(de_annot$dmel_entrez)]))
cat("\nDistinct Drosophila genes hit by >1 tested Diatraea gene:",
    sum(table(unlist(strsplit(na.omit(genes$dmel_entrez[genes$in_deseq2]), ","))) > 1), "\n")
cat("\nTop DE genes:\n")
print(head(de_annot %>% dplyr::select(gene, direction, log2FC_infected_vs_control, preferred_name, KEGG_ko,
                               dmel_symbol, dmel_orth_type), 25), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
