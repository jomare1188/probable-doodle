#!/usr/bin/env Rscript
# Step 0 - QC and harmonization. Every table is put in the direction infected vs control.
# Nothing published is re-computed: DESeq2 statistics and metabolite p-values are carried over as is.
# Run from the repository root:  Rscript scripts/step0_qc_harmonize.R
suppressPackageStartupMessages(library(readxl))

out_dir  <- "results/step0"
dea_up   <- "data/transcriptome/up_regulated.csv"
dea_down <- "data/transcriptome/down_regulated.csv"
dea_full <- "data/transcriptome/dea_all_genes.csv"
met_xlsx <- "data/metabolome/metabolome.xlsx"
met_pub  <- "data/metabolome/statistical_results.csv"
mb_ready <- "data/microbiome/abundance_table_ready.tsv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
outlier_fold <- 5

## ---------------------------------------------------------------- transcriptome
# The gene ID column has no header, so the header has one field less than the rows:
# read.csv then uses the first column as row names, which is what we want.
read_deseq <- function(f, strict = TRUE) {
  x <- read.csv(f, check.names = FALSE)
  stopifnot(all(x$sampleA == "control"), all(x$sampleB == "infected"))
  # published log2FoldChange = log2(control / infected); the base means below verify this.
  # The model is ~ W_1 + group (RUVs), so for genes with LFC ~ 0 the raw base-mean ratio
  # can have the other sign: strict for the DE lists, only counted for the full table.
  disagree <- sign(log2(x$baseMeanA / x$baseMeanB)) != sign(x$log2FoldChange)
  if (strict) stopifnot(!any(disagree)) else
    cat("genes whose base-mean ratio and model LFC disagree in sign:", sum(disagree),
        "| max |LFC| among them:", signif(max(abs(x$log2FoldChange[disagree])), 3), "\n")
  data.frame(gene_id = rownames(x),
             baseMean_control = x$baseMeanA, baseMean_infected = x$baseMeanB,
             baseMean = x$baseMean,
             log2FC_infected_vs_control = -x$log2FoldChange,
             lfcSE = x$lfcSE, stat_infected_vs_control = -x$stat,
             pvalue = x$pvalue, padj = x$padj, row.names = NULL)
}
direction <- function(l) ifelse(l > 0, "up_in_infected", "down_in_infected")

de <- rbind(read_deseq(dea_up), read_deseq(dea_down))
de$direction <- direction(de$log2FC_infected_vs_control)
stopifnot(!anyDuplicated(de$gene_id), all(de$padj < 0.05))
write.csv(de, file.path(out_dir, "de_genes_infected_vs_control.csv"), row.names = FALSE)

cat("== Transcriptome (infected vs control) ==\n")
cat("up_regulated.csv rows:", sum(de$direction == "down_in_infected"),
    "-> down in infected after the flip\n")
cat("down_regulated.csv rows:", sum(de$direction == "up_in_infected"),
    "-> up in infected after the flip\n")
print(table(de$direction))

full <- read_deseq(dea_full, strict = FALSE)
full$direction <- direction(full$log2FC_infected_vs_control)
full$significant <- !is.na(full$padj) & full$padj < 0.05
write.csv(full, file.path(out_dir, "deseq2_all_genes_infected_vs_control.csv"), row.names = FALSE)
cat("Full DESeq2 table:", nrow(full), "genes;", sum(full$significant), "with padj < 0.05;",
    sum(is.na(full$padj)), "with padj NA\n")
cat("DE lists identical to padj<0.05 of the full table:",
    setequal(full$gene_id[full$significant], de$gene_id), "\n")

## ---------------------------------------------------------------- metabolome
met <- as.data.frame(read_excel(met_xlsx), check.names = FALSE)
met$group <- ifelse(met$TREATMENT == "non-contaminated", "control", "infected")
cat("\n== Metabolome ==\n"); print(table(met$TREATMENT, met$group))
mets <- setdiff(names(met), c("TREATMENT", "REP", "group"))
stopifnot(length(mets) == 18, all(table(met$group) == 6))

# published per-metabolite statistics (names there are snake_case versions of the column names)
key <- function(x) sub("^x(?=[0-9])", "", gsub("_+$", "", gsub("[^a-z0-9]+", "_", tolower(x))), perl = TRUE)
pub <- read.csv(met_pub)
pub$key <- key(pub$metabolite)
stopifnot(setequal(pub$key, key(mets)))

long <- do.call(rbind, lapply(mets, function(m)
  data.frame(metabolite = m, group = met$group, rep = met$REP, value = met[[m]])))
long$group_median <- ave(long$value, long$metabolite, long$group, FUN = median)
long$fold_from_median <- ifelse(long$group_median > 0 & long$value > 0,
                                long$value / long$group_median, NA)
long$flag_outlier <- !is.na(long$fold_from_median) &
  (long$fold_from_median > outlier_fold | long$fold_from_median < 1 / outlier_fold)
write.csv(long, file.path(out_dir, "metabolome_long_with_flags.csv"), row.names = FALSE)

summ <- do.call(rbind, lapply(mets, function(m) {
  ctl <- met[[m]][met$group == "control"]; inf <- met[[m]][met$group == "infected"]
  only_inf <- all(ctl == 0) && any(inf > 0)
  p <- pub[pub$key == key(m), ]
  data.frame(metabolite = m,
             mean_control = mean(ctl), mean_infected = mean(inf),
             median_control = median(ctl), median_infected = median(inf),
             n_zero_control = sum(ctl == 0), n_zero_infected = sum(inf == 0),
             detected_only_in_infected = only_inf,
             log2FC_mean_infected_vs_control = if (only_inf) NA else log2(mean(inf) / mean(ctl)),
             log2FC_median_infected_vs_control = if (only_inf) NA else log2(median(inf) / median(ctl)),
             direction = if (mean(inf) > mean(ctl)) "up_in_infected" else "down_in_infected",
             n_flagged_values = sum(long$flag_outlier[long$metabolite == m]),
             published_p = p$p_value, published_p_adj = p$p_adj,
             published_significance = p$significance,
             published_direction = ifelse(p$log2_fc > 0, "up_in_infected", "down_in_infected"))
}))
summ$differential_published <- summ$published_p_adj < 0.05
write.csv(summ, file.path(out_dir, "metabolome_summary_infected_vs_control.csv"), row.names = FALSE)

options(width = 200)
print(format(summ[, c("metabolite", "detected_only_in_infected", "log2FC_mean_infected_vs_control",
                      "log2FC_median_infected_vs_control", "direction", "n_flagged_values",
                      "published_p_adj", "differential_published", "published_direction")],
             digits = 3), row.names = FALSE)
cat("\nValues more than", outlier_fold, "x from their group median (reported, not changed):\n")
print(format(long[long$flag_outlier, c("metabolite", "group", "rep", "value", "group_median",
                                       "fold_from_median")], digits = 4), row.names = FALSE)

## ---------------------------------------------------------------- microbiome (input description only)
ab <- read.delim(mb_ready, check.names = FALSE)
samples <- setdiff(names(ab), c("tax", "total"))
cat("\n== Microbiome ==\nsamples:", samples, "\n")
bad_total <- rowSums(ab[, samples]) != ab$total
cat("rows where 'total' differs from the sum of the six samples:", sum(bad_total), "\n")
print(cbind(ab[bad_total, ], row_sum = rowSums(ab[bad_total, samples])), row.names = FALSE)
missing_tax <- is.na(ab$tax) | ab$tax == ""
tax <- do.call(rbind, strsplit(ab$tax[!missing_tax], ";"))
cat("rows:", nrow(ab), "| missing taxonomy:", sum(missing_tax),
    "(reads:", sum(ab$total[missing_tax]), ")\n")
cat("rank 1:\n"); print(table(tax[, 1]))
cat("rank 2 within Eukaryota:\n"); print(table(tax[tax[, 1] == "Eukaryota", 2]))
cat("duplicated genus names:", sum(duplicated(tax[, 7])), "\n")
cat("library sizes:\n"); print(colSums(ab[, samples]))
cat("smallest non-zero count:", min(ab[, samples][ab[, samples] > 0]), "\n")
