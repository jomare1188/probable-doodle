# PCA plots of the RUVr, RUVg and RUVs corrections tried before choosing RUVs k = 1.
# Same code as the exploration part of the first ruv.r, on the run's own nf-core dds.
#
# Usage: Rscript ruv_exploration.r <star_salmon_dir> <samplesheet.csv> <outdir>
#   defaults: ..  ../../../../raw_reads/samples_paired_relaxed.csv  .
# Env: conda R_popstat_jorge

library("RUVSeq")
library("DESeq2")
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)
star_salmon <- if (length(args) >= 1) args[1] else ".."
samplesheet <- if (length(args) >= 2) args[2] else "../../../../raw_reads/samples_paired_relaxed.csv"
outdir      <- if (length(args) >= 3) args[3] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load(file.path(star_salmon, "deseq2_qc", "deseq2.dds.RData"))
metadata <- read.table(samplesheet, header = TRUE, sep = ",")
metadata <- metadata[match(colnames(dds), metadata$sample), ]
colors <- c("#FF0000", "#00A08A")  # wesanderson Darjeeling1, 2 colours

plot_pca <- function(normalized_counts, file) {
  d <- DESeqDataSetFromMatrix(countData = normalized_counts, colData = metadata, design = ~ group)
  pca_data <- plotPCA(varianceStabilizingTransformation(d), intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    scale_colour_manual(values = colors) +
    theme_bw(base_size = 22)
  ggsave(p, filename = file.path(outdir, file), units = "cm", width = 15*1.3, height = 15, dpi = 320)
}

# ---- RUVr ----
design <- model.matrix(~ dds$Group1)
y <- DGEList(counts = counts(dds), group = dds$Group1)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type = "deviance")
for (k in 1:2) {
  plot_pca(RUVr(y$counts, rownames(y), k = k, res)$normalizedCounts, paste0("k", k, "_RUVr_groups.png"))
}

# ---- RUVg: control genes = 1% of genes with the lowest coefficient of variation ----
counts_mat <- counts(dds, normalized = TRUE)
gene_cv <- apply(counts_mat, 1, sd) / rowMeans(counts_mat)
n_control <- ceiling(0.01 * nrow(counts_mat))
control_genes <- rownames(counts_mat)[order(gene_cv, decreasing = FALSE)[1:n_control]]
plot_pca(RUVg(counts(dds), control_genes, k = 2)$normalizedCounts, "RUVg_groups.png")

# ---- RUVs ----
replicates <- makeGroups(dds$Group1)
for (k in 1:4) {
  plot_pca(RUVs(counts(dds), rownames(dds), replicates, k = k)$normalizedCounts, paste0("k", k, "_RUVs_groups.png"))
}
