# Differential expression + GO enrichment for t3_vs_t4
# Same methodology as run_sugarcane_diatrea/star_salmon/deseq2_qc/ruv.r
# RUVs, RUVg and RUVr (k = 1); each method writes its results to its own folder
# T3 = Leaf_120h_Fusarium, T4 = Leaf_120h_Diatrea+Fusarium

library("RUVSeq")
library("DESeq2")
library("tximport")
library("tidyverse")

base_dir <- "/dados04/jorge/rnaseq_diatraea"
run_dir  <- file.path(base_dir, "rnaseq/t3_vs_t4/star_salmon")

groupA <- "Leaf_120h_Fusarium"
groupB <- "Leaf_120h_Diatrea+Fusarium"
k <- 1

metadata <- read.table(file.path(base_dir, "raw_reads/samples_t3_vs_t4.csv"), header = TRUE, sep = ",")

load("deseq2.dds.RData")
stopifnot(colnames(dds) == metadata$sample)

colors = c("#FF0000", "#00A08A") # wesanderson Darjeeling1

plot_ruv_pca <- function(normalizedCounts, filename) {
  dataRUV <- DESeqDataSetFromMatrix(countData = normalizedCounts,
                                    colData = metadata,
                                    design = ~ Group)
  RUV_vst <- varianceStabilizingTransformation(dataRUV)
  pca_data <- plotPCA(RUV_vst, intgroup = "Group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    scale_colour_manual(values = colors) +
    theme_bw(base_size=22)
  ggsave(p, filename = filename, units = "cm", width = 15*1.3, height = 15, dpi = 320)
}

# ---- RUVs ----
replicates <- makeGroups(metadata$Group)
set_RUVs <- RUVs(counts(dds), rownames(dds), replicates, k = k)

# ---- RUVg ----
# ---- Select control genes by CV ----
counts_mat <- counts(dds, normalized = TRUE)  # normalized counts
gene_means <- rowMeans(counts_mat)
gene_sds   <- apply(counts_mat, 1, sd)
gene_cv    <- gene_sds / gene_means  # coefficient of variation
# rank genes by CV (ascending: stable first)
gene_rank <- order(gene_cv, decreasing = FALSE)
# take top 1% most stable genes
n_control <- ceiling(0.01 * nrow(counts_mat))
control_genes <- rownames(counts_mat)[gene_rank[1:n_control]]
length(control_genes)
set_RUVg <- RUVg(counts(dds), control_genes, k = k)

# ---- RUVr ----
design_edger <- model.matrix(~metadata$Group)
y <- DGEList(counts=counts(dds), group=metadata$Group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design_edger)
y <- estimateGLMTagwiseDisp(y, design_edger)
fit <- glmFit(y, design_edger)
res_dev <- residuals(fit, type="deviance")
set_RUVr <- RUVr(y$counts, rownames(y), k = k, res_dev)

ruv_sets <- list(RUVs = set_RUVs, RUVg = set_RUVg, RUVr = set_RUVr)

# Make Differential expression analysis
# load files paths
sample_files = file.path(run_dir, pull(metadata, "sample"), "quant.sf")
# name table columns
names(sample_files) = pull(metadata, "sample")
# relate genes to transcripts
tx2gene = read.table(file.path(run_dir, "salmon.merged.tx2gene.tsv"), header = T)
# GENE MODE
# import count data to tximport
count_data = tximport(files = sample_files,
        type = "salmon",
        tx2gene = tx2gene,
        ignoreTxVersion = F,
        ignoreAfterBar = T)

raw <- DESeqDataSetFromTximport(txi = count_data,
        colData = metadata,
        design = ~ Group)

dim(raw)
temp <- as.data.frame(counts(raw))
logic <- (apply(temp, c(1,2), function(x){x>0}))
filter_genes <- rowSums(logic) > 3
fi <- raw[filter_genes,]
dim(fi)
temp <- NULL

run_dea <- function(W, outdir) {
  stopifnot(rownames(W) == colnames(fi))
  fi_w <- fi
  colData(fi_w)$W_1 <- W[,1]
  design(fi_w) <- ~ W_1 + Group

  ### Differencial expression analyses
  dea <- DESeq(fi_w)
  dea_contrast <- results(dea, lfcThreshold = 1, altHypothesis = "greaterAbs", contrast = c("Group", groupA, groupB))
  dea_df <- as.data.frame(dea_contrast)

  baseMeanA <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == groupA])
  baseMeanB <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$Group == groupB])

  res = cbind(baseMeanA, baseMeanB, dea_df)
  res = cbind(sampleA = groupA, sampleB = groupB, as.data.frame(res))
  res = res[complete.cases(res),]

  final <- res %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
  final <- final[rev(order(final$log2FoldChange)),]

  up <- final %>% filter(log2FoldChange > 0)
  down <- final %>% filter(log2FoldChange < 0)
  cat(outdir, "- Up-regulated:", nrow(up), " Down-regulated:", nrow(down), "\n")

  write.table(up, file.path(outdir, "up_regulated.csv"), sep = ",", quote = F)
  write.table(down, file.path(outdir, "down_regulated.csv"), sep = ",", quote = F)
  list(up = up, down = down)
}


# GO enrichment
library(topGO)

GO <- read.table(file.path(base_dir, "reference_genomes/sugarcane/annotation/go_table_locus_name.txt"), header=FALSE, stringsAsFactors=FALSE, colClasses = c("character", "character"))

colnames(GO) <- c("Gene", "GO_term")
GO$GO_term <- paste0("GO:", GO$GO_term)

# Group by Gene and aggregate GO terms into a single string separated by spaces
formatted_GO <- aggregate(GO_term ~ Gene, GO, function(x) paste(x, collapse=" "))

gene2GO <- strsplit(formatted_GO$GO_term, " ")
names(gene2GO) <- formatted_GO$Gene
geneNames <- names(gene2GO)

go_header <- data.frame(GO.ID=character(), Term=character(), Annotated=integer(), Significant=integer(),
                        Expected=numeric(), Classic=character(), p.adj=numeric())

run_go <- function(deg, prefix) {
  # select set of genes to make overrepresentation test
  MyInterestingGenes <- sub("\\.v2\\.1$", "", rownames(deg))
  geneList <- factor(as.integer(geneNames %in% MyInterestingGenes), levels = c(0, 1))
  names(geneList) <- geneNames
  if (sum(geneList == 1) == 0) {
    message(prefix, ": no DE genes with GO annotation, skipping")
    write.table(go_header, file = paste0(prefix, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")
    return(invisible(NULL))
  }
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO)
  allGO = usedGO(GOdata)

  Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  # Make results table
  table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
  # Filter not significant values for classic algorithm
  table1 <- filter(table, Classic < 0.05)
  # Performing BH correction on our p values FDR
  p.adj <- round(p.adjust(table1$Classic, method="BH"), digits = 4)
  # Create the file with all the statistics from GO analysis
  all_res_final <- cbind(table1, p.adj)
  all_res_final <- all_res_final[order(all_res_final$p.adj),]
  # Get list of significant GO after multiple testing correction
  results.table.bh = all_res_final[which(all_res_final$p.adj <= 0.05),]
  write.table(results.table.bh, file = paste0(prefix, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")
  cat(prefix, ": ", nrow(results.table.bh), " significant GO terms\n", sep = "")
  if (nrow(results.table.bh) == 0) return(invisible(NULL))

  ntop <- 24
  ggdata <- results.table.bh[1:min(ntop, nrow(results.table.bh)),]
  ggdata <- ggdata[complete.cases(ggdata), ]
  ggdata$p.adj <- as.numeric(ggdata$p.adj)
  ggdata <- ggdata[order(ggdata$p.adj),]
  ggdata$Term <- factor(ggdata$Term, levels = rev(unique(ggdata$Term))) # fixes order

  gg1 <- ggplot(ggdata, aes(x = Term, y = -log10(p.adj), size = Significant)) +
    geom_point(colour = "black") +
    scale_size(range = c(2.5, 12.5)) +
    xlab('GO Term') +
    ylab('-log(p)') +
    labs(title = 'GO Biological processes', size = 'Significant') +
    theme_bw(base_size = 24) +
    coord_flip()

  ggsave(paste0(prefix, ".svg"), plot = gg1, device = "svg", width = 40, height = 30, dpi = 300, units = "cm")
  ggsave(paste0(prefix, ".pdf"), plot = gg1, device = "pdf", width = 40, height = 30, dpi = 300, units = "cm")
}


for (method in names(ruv_sets)) {
  set_RUV <- ruv_sets[[method]]
  outdir <- paste0(method, "_k", k)
  dir.create(outdir, showWarnings = FALSE)
  plot_ruv_pca(set_RUV$normalizedCounts, file.path(outdir, paste0("k", k, "_", method, "_groups.png")))
  deg <- run_dea(set_RUV$W, outdir)
  run_go(deg$up, file.path(outdir, "GO_up"))
  run_go(deg$down, file.path(outdir, "GO_down"))
}
