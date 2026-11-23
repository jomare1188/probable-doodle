#!/usr/bin/env Rscript
# Independent validation of steps 0-3: recomputes key numbers from the raw inputs and compares
# them with the saved outputs. Prints PASS / FAIL / NOTE per test; stops if any test fails.
# Run from the repository root with the integration_diatrea env:  Rscript scripts/step5a_validate_outputs.R
suppressPackageStartupMessages({ library(readxl); library(dplyr) })
res <- list()
test <- function(id, ok, msg) {
  status <- if (isTRUE(ok)) "PASS" else if (identical(ok, NA)) "NOTE" else "FAIL"
  res[[length(res) + 1]] <<- data.frame(test = id, status = status, detail = msg)
  cat(sprintf("[%s] %-8s %s\n", status, id, msg))
}
S <- c("control_1", "control_2", "control_3", "infected_1", "infected_2", "infected_3")

## ============================================================ step 0 - transcriptome
up   <- read.csv("data/transcriptome/up_regulated.csv")      # gene IDs become row names (header is one short)
down <- read.csv("data/transcriptome/down_regulated.csv")
full <- read.csv("data/transcriptome/dea_all_genes.csv")
de   <- read.csv("results/step0/de_genes_infected_vs_control.csv")
test("0.T1", nrow(up) == 331 && nrow(down) == 143, sprintf("rows: up_regulated %d, down_regulated %d", nrow(up), nrow(down)))
test("0.T2", all(up$log2FoldChange > 0) && all(down$log2FoldChange < 0),
    "published LFC > 0 in up_regulated and < 0 in down_regulated (control/infected)")
test("0.T3", all(up$baseMeanA > up$baseMeanB) && all(down$baseMeanA < down$baseMeanB),
    "base means: up_regulated = higher in control (baseMeanA)")
pub <- rbind(data.frame(g = rownames(up), lfc = up$log2FoldChange), data.frame(g = rownames(down), lfc = down$log2FoldChange))
m <- merge(pub, de, by.x = "g", by.y = "gene_id")
test("0.T4", nrow(m) == 474 && all(abs(m$log2FC_infected_vs_control + m$lfc) < 1e-12),
    "harmonised log2FC = -published log2FC for all 474 DE genes")
test("0.T5", sum(de$direction == "up_in_infected") == 143 && sum(de$direction == "down_in_infected") == 331,
    "143 up / 331 down in infected")
test("0.T6", setequal(rownames(full)[full$padj < 0.05], de$gene_id), "DE lists == padj < 0.05 in the full DESeq2 table")
test("0.T7", all(abs(full[de$gene_id, "log2FoldChange"] + de$log2FC_infected_vs_control) < 1e-12),
    "DE values equal the full-table values")

## ============================================================ step 0 - metabolome
x <- as.data.frame(read_excel("data/metabolome/metabolome.xlsx"), check.names = FALSE)
ms <- read.csv("results/step0/metabolome_summary_infected_vs_control.csv", check.names = FALSE)
mets <- setdiff(names(x), c("TREATMENT", "REP"))
ctl <- x$TREATMENT == "non-contaminated"
l2 <- sapply(mets, function(k) if (all(x[ctl, k] == 0)) NA else log2(mean(x[!ctl, k]) / mean(x[ctl, k])))
test("0.M1", sum(ctl) == 6 && sum(!ctl) == 6 && length(mets) == 18, "6 control, 6 infected, 18 metabolites")
test("0.M2", isTRUE(all.equal(unname(l2[ms$metabolite]), ms$log2FC_mean_infected_vs_control)),
    "log2FC (means) recomputed from the xlsx matches")
only <- mets[sapply(mets, function(k) all(x[ctl, k] == 0) && any(x[!ctl, k] > 0))]
test("0.M3", setequal(only, ms$metabolite[ms$detected_only_in_infected]),
    paste("detected only in infected:", paste(sort(only), collapse = ", ")))
test("0.M4", all(ms$direction == ms$published_direction), "own direction == published direction for all 18")
test("0.M5", sum(ms$differential_published) == 15, "15 differential (published p_adj < 0.05)")

## ============================================================ step 0a - cascade
pr <- read.delim(gzfile("results/step0a_microbiome_cascade/per_read_final_assignment.tsv.gz"))
filt <- read.delim("results/step0a_microbiome_cascade/read_set_filtering.tsv")
test("0a.1", !anyDuplicated(pr[, c("sample", "read_id")]), "each read assigned exactly once")
test("0a.2", all(table(pr$sample)[filt$sample] == filt$reads_100_2000bp),
    sprintf("reads per sample == filtered read set (total %d)", nrow(pr)))
ct <- read.delim("results/step0a_microbiome_cascade/abundance_table_genus_cascade.tsv", check.names = FALSE)
micro <- table(pr$sample[pr$category %in% c("Bacteria", "Archaea", "Fungi")])[S]
test("0a.3", all(colSums(ct[, S]) == micro), "genus table column sums == microbial reads per sample")
test("0a.4", all(rowSums(ct[, S]) == ct$total) && !anyDuplicated(ct$genus), "total == row sums; genus names unique")
pez <- c("Nanoscypha", "Wynnea", "Helvella", "Otidea", "Galiella", "Pseudopithyella", "Kompsoscypha")
test("0a.5", !any(ct$genus %in% pez), "no Pezizales artefact genera in the cascade table")
fu <- unlist(ct[ct$genus == "Fusarium", S])
test("0a.6", all(fu[4:6] > fu[1:3]), paste("Fusarium:", paste(fu, collapse = "/")))

## ============================================================ step 0b - DA
da <- read.csv("results/step0b_microbiome_DA/microbiome_DA.csv")
prev <- ct$genus[rowSums(ct[, S] > 0) >= 3]
test("0b.1", setequal(prev, da$genus), sprintf("prevalence filter (>=3 of 6) reproduces the %d genera", length(prev)))
agree <- ifelse(!is.na(da$ancombc2_structural_zero) &
                  ((da$ancombc2_structural_zero == "absent in control" & da$log2_effect_infected_vs_control > 0) |
                   (da$ancombc2_structural_zero == "absent in infected" & da$log2_effect_infected_vs_control < 0)), TRUE,
                !is.na(da$ancombc2_log2fc) & sign(da$ancombc2_log2fc) == sign(da$log2_effect_infected_vs_control))
rule <- (da$BH < 0.1 | (abs(da$aldex2_effect_size) > 1 & agree)) & da$feature_type == "genus"
test("0b.2", all(rule == da$differential), sprintf("rule re-applied gives the same %d differential genera", sum(rule)))
rel <- sweep(as.matrix(ct[, S]) + 0.5, 2, colSums(ct[, S]), "/"); rownames(rel) <- ct$genus
lr <- log2(rowMeans(rel[da$genus, 4:6]) / rowMeans(rel[da$genus, 1:3]))
dd <- da$differential
test("0b.3", mean(sign(lr[dd]) == sign(da$log2_effect_infected_vs_control[dd])) > 0.95,
    sprintf("ALDEx2 direction agrees with raw relative-abundance ratio for %.1f%% of differential genera",
            100 * mean(sign(lr[dd]) == sign(da$log2_effect_infected_vs_control[dd]))))
test("0b.4", all(c("genus", "phylum", "kingdom", "log2_effect_infected_vs_control", "p", "BH", "method", "direction") %in% names(da)),
    "microbiome_DA.csv has the required columns")
t12 <- da$differential & grepl("^(1|2)", da$confidence_tier)
test("0b.5", sum(t12 & da$kingdom == "Bacteria") == 114, sprintf("tier 1-2 bacterial genera: %d", sum(t12 & da$kingdom == "Bacteria")))
test("0b.6", !any(da$differential & da$feature_type != "genus"), "no _uncl. group called differential")

## ============================================================ step 1 - annotation
ga <- read.csv("results/step1_gene_annotation/de_genes_annotated.csv")
test("1.1", setequal(ga$gene, de$gene_id), "annotated DE table contains exactly the 474 DE genes")
emap <- read.delim("data/annotation/proteins.emapper.emapper.annotations.gz",
                   comment.char = "#", quote = "")
g2p <- read.table("data/annotation/gene2protein.txt", col.names = c("gene", "protein"))
ko_raw <- tapply(gsub("ko:", "", emap$KEGG_ko[match(g2p$protein, emap$query)]), g2p$gene, function(v) {
  v <- unique(unlist(strsplit(v[!is.na(v) & v != "-"], ","))); if (length(v)) paste(sort(v), collapse = ",") else NA })
norm <- function(v) ifelse(is.na(v) | v == "", "<none>", as.character(v))   # CSV "" vs R NA
ok <- norm(ga$KEGG_ko) == norm(unname(ko_raw[ga$gene]))
test("1.2", all(ok), sprintf("KOs of all DE genes re-derived from the raw emapper file (%d/%d match)", sum(ok), length(ok)))
test("1.3", round(100 * mean(!is.na(ga$KEGG_ko) & ga$KEGG_ko != ""), 1) == 60.5 &&
          round(100 * mean(!is.na(ga$dmel_entrez) & ga$dmel_entrez != ""), 1) == 61.0,
    "coverage: 60.5% KO, 61.0% Drosophila ortholog")

## ============================================================ step 2 - metabolite IDs
ids <- read.csv("results/step2_metabolite_ids/metabolite_ids_final.csv")
kg <- ids$KEGG[ids$KEGG != ""]
test("2.1", !anyDuplicated(kg) && all(grepl("^C\\d{5}$", kg)), sprintf("%d KEGG IDs, unique and well-formed", length(kg)))
names_kegg <- sapply(kg, function(k) {
  f <- list.files("resources/web_cache/step2_pubchem_kegg", pattern = paste0("kegg_get_", k), full.names = TRUE)
  if (length(f)) sub("^NAME\\s+", "", grep("^NAME", readLines(f[1]), value = TRUE)[1]) else NA })
cat("        KEGG names:", paste(sprintf("%s=%s", ids$metabolite[ids$KEGG != ""], names_kegg), collapse = " | "), "\n")
test("2.2", !any(is.na(names_kegg)), "every accepted KEGG ID has a cached KEGG record (names printed above)")
test("2.3", all(ids$KEGG[ids$metabolite == "kynurenic acid"] == "C01717"), "kynurenic acid: wrong PubChem xref C14645 excluded")

## ============================================================ step 3 - multiGSEA
mg <- read.csv("results/step3_multigsea/multigsea_kegg_dmel.csv")
txr <- read.csv("results/step3_multigsea/transcriptome_rank_dmel_CG.csv")
st <- setNames(full$stat, rownames(full))
test("3.1", all(abs(txr$stat + st[txr$gene]) < 1e-12), "transcriptome ranks == -published Wald stat")
mr <- read.csv("results/step3_multigsea/metabolome_rank_kegg.csv")
test("3.2", all(sign(mr$rank) == ifelse(mr$direction == "up_in_infected", 1, -1)), "metabolite rank signs == direction")
z <- function(p) qnorm(p / 2, lower.tail = FALSE)            # Stouffer on two-sided p-values
st_re <- apply(mg[, c("transcriptome_pval", "metabolome_pval")], 1, function(p) {
  p <- p[!is.na(p)]; if (!length(p)) NA else 2 * pnorm(sum(z(p)) / sqrt(length(p)), lower.tail = FALSE) })
st_re1 <- apply(mg[, c("transcriptome_pval", "metabolome_pval")], 1, function(p) {
  p <- p[!is.na(p)]; if (!length(p)) NA else pnorm(sum(qnorm(p, lower.tail = FALSE)) / sqrt(length(p)), lower.tail = FALSE) })
same2 <- isTRUE(all.equal(st_re, mg$combined_pval_stouffer, tolerance = 1e-6))
same1 <- isTRUE(all.equal(st_re1, mg$combined_pval_stouffer, tolerance = 1e-6))
test("3.3", same1 || same2, sprintf("Stouffer recomputed by hand matches multiGSEA (%s formula)",
                                    if (same1) "one-sided z = qnorm(1-p)" else "two-sided"))
two <- mg[!is.na(mg$metabolome_pval), ]
disc <- two$pathway[sign(two$transcriptome_NES) != sign(two$metabolome_NES)]
test("3.4", if (length(disc)) NA else TRUE,
    paste("pathways where the two layers point in OPPOSITE directions:", paste(disc, collapse = "; ")))
test("3.5", sum(mg$combined_padj < 0.05, na.rm = TRUE) == 6, "6 pathways with combined padj < 0.05")

out <- do.call(rbind, res)
write.csv(out, "results/validation_steps0-3.csv", row.names = FALSE)
cat("\nSummary:", paste(names(table(out$status)), table(out$status), collapse = ", "), "\n")
if (any(out$status == "FAIL")) stop("validation failed: ", paste(out$test[out$status == "FAIL"], collapse = ", "))
