# RUVs (k = 1) correction + DESeq2 differential expression, control vs infected.
# Same method as run_paired_samples/star_salmon/deseq2_qc/ruv.r (the RUVs k=1 branch),
# with relative paths so it can be run on any nf-core/rnaseq star_salmon output.
#
# Usage: Rscript ruv_dea.r <star_salmon_dir> <samplesheet.csv> <outdir>
#   defaults: ..  ../../../../raw_reads/samples_paired_relaxed.csv  .
# Env: conda R_popstat_jorge (RUVSeq 1.40, DESeq2 1.46, same as the strict analysis;
#      DESeq2 >= 1.44 changed the greaterAbs test, older versions give 2x larger p-values)

library("RUVSeq")
library("DESeq2")
library("tximport")
library("tidyverse")
library("BiocParallel")

args <- commandArgs(trailingOnly = TRUE)
star_salmon <- if (length(args) >= 1) args[1] else ".."
samplesheet <- if (length(args) >= 2) args[2] else "../../../../raw_reads/samples_paired_relaxed.csv"
outdir      <- if (length(args) >= 3) args[3] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
register(MulticoreParam(8))

# nf-core DESEQ2_QC object (Group1 = control / infected, from the sample names)
load(file.path(star_salmon, "deseq2_qc", "deseq2.dds.RData"))

vst <- vst(dds)
write.table(counts(dds), file = file.path(outdir, "raw_counts.csv"), quote = F, col.names = T, row.names = T, sep = ",")
write.table(assay(vst), file = file.path(outdir, "vst_transform.csv"), quote = F, col.names = T, row.names = T, sep = ",")

metadata <- read.table(samplesheet, header = TRUE, sep = ",")

# ---- RUVs k = 1 ----
replicates <- makeGroups(dds$Group1)
set_RUVs <- RUVs(counts(dds), rownames(dds), replicates, k = 1)

# colData in the same sample order as the count matrix
meta_ruv <- metadata[match(colnames(set_RUVs$normalizedCounts), metadata$sample), ]
dataRUVs <- DESeqDataSetFromMatrix(countData = set_RUVs$normalizedCounts,
                                   colData = meta_ruv,
                                   design = ~ group)
RUVs_vst <- varianceStabilizingTransformation(dataRUVs)
write.table(set_RUVs$W, file = file.path(outdir, "RUVs_k1_W.tsv"), sep = "\t", quote = F)
write.table(assay(dataRUVs), file = file.path(outdir, "RUVs_k1_normalized_counts.tsv"), sep = "\t", quote = F)

# PCA for RUVs
colors <- c("#FF0000", "#00A08A")  # wesanderson Darjeeling1, 2 colours
pca_data_s <- plotPCA(RUVs_vst, intgroup = "group", returnData = TRUE)
percentVar_s <- round(100 * attr(pca_data_s, "percentVar"))

p_s <- ggplot(pca_data_s, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_s[1], "%")) +
  ylab(paste0("PC2: ", percentVar_s[2], "%")) +
  scale_colour_manual(values = colors) +
  theme_bw(base_size = 22)

ggsave(p_s, filename = file.path(outdir, "k1_RUVs_groups.png"), units = "cm", width = 15*1.3, height = 15, dpi = 320)

# ---- Differential expression ----
sample_files <- file.path(star_salmon, pull(metadata, "sample"), "quant.sf")
names(sample_files) <- pull(metadata, "sample")
tx2gene <- read.table(file.path(star_salmon, "tx2gene.tsv"), header = T)

count_data <- tximport(files = sample_files,
                       type = "salmon",
                       tx2gene = tx2gene,
                       ignoreTxVersion = F,
                       ignoreAfterBar = T)

raw <- DESeqDataSetFromTximport(txi = count_data,
                                colData = metadata,
                                design = ~ group)

# keep genes with counts > 0 in more than 3 samples
filter_genes <- rowSums(counts(raw) > 0) > 3
fi <- raw[filter_genes, ]
message("genes before / after filter: ", nrow(raw), " / ", nrow(fi))

W <- set_RUVs$W
stopifnot(rownames(W) == colnames(fi))
colData(fi)$W_1 <- W[, 1]
design(fi) <- ~ W_1 + group

dea <- DESeq(fi, parallel = T)
dea_contrast <- results(dea, lfcThreshold = 1, altHypothesis = "greaterAbs", parallel = T,
                        contrast = c("group", "control", "infected"))
dea_df <- as.data.frame(dea_contrast)

baseMeanA <- rowMeans(counts(dea, normalized = TRUE)[, colData(dea)$group == "control"])
baseMeanB <- rowMeans(counts(dea, normalized = TRUE)[, colData(dea)$group == "infected"])

res <- cbind(baseMeanA, baseMeanB, dea_df)
res <- cbind(sampleA = "control", sampleB = "infected", as.data.frame(res))
res <- res[complete.cases(res), ]
write.table(res, file.path(outdir, "dea_all_genes.csv"), sep = ",", quote = F)

final <- res %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
final <- final[rev(order(final$log2FoldChange)), ]

# up = higher in control, down = higher in infected
up <- final %>% filter(log2FoldChange > 0)
down <- final %>% filter(log2FoldChange < 0)

write.table(up, file.path(outdir, "up_regulated.csv"), sep = ",", quote = F)
write.table(down, file.path(outdir, "down_regulated.csv"), sep = ",", quote = F)
writeLines(rownames(up), file.path(outdir, "genes_up.txt"))
writeLines(rownames(down), file.path(outdir, "genes_down.txt"))

message("up-regulated: ", nrow(up), "   down-regulated: ", nrow(down))
