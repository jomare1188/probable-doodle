library("RUVSeq")
library("DESeq2")
library("tximport")
library("tidyverse")
library("wesanderson")

load("deseq2.dds.RData")

vst <- vst(dds)

# write raw counts
raw_counts <- counts(dds)

write.table(raw_counts, file = "raw_counts.csv", quote = F, col.names = T, row.names = T, sep = ",")
write.table(assay(vst), file = "vst_transform.csv", quote = F, col.names = T, row.names = T, sep = ",")



metadata <- read.table("../../../../raw_reads/samples.csv", header = TRUE, sep = ",")

# RUVr
# some magic
design <- model.matrix(~vst$Group1)
y <- DGEList(counts=counts(dds), group=vst$Group1)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# try for different k values
set_RUVr <- RUVr(y$counts,rownames(y), k=2, res)
dataRUVr <- DESeqDataSetFromMatrix(countData = set_RUVr$normalizedCounts,
                                   colData = metadata,
                                   design = ~ group)

RUVr_vst <- varianceStabilizingTransformation(dataRUVr)

# plot PCA for RUVr(normalized) and vst(transformed) data
colors = wes_palette("Darjeeling1",2 , type = "discrete")
pca_data <- plotPCA(RUVr_vst, intgroup = "group", returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))
# plot R color condition
colors = wes_palette("Darjeeling1",2 , type = "discrete")
# Plooot 
p <- ggplot(pca_data, aes(x = PC1, y = PC2,  color = as.factor(group)))
p <- p + geom_point(size=3)
p <- p + labs() 
p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
#p <- p + geom_text_repel(aes(label=name), size=2.5, max.overlaps = 20 )
p <- p + scale_colour_manual(values = colors)
#p <- p + scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8))
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman", size=22),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))
ggsave(p, filename = "k2_RUVr_groups.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

#write.table(assay(dataRUVr), file = "../results/matrices/sorghum_RUVr_k1_raw_counts.tsv", sep = "\t" ,quote = F, col.names = T)
#write.table(assay(RUVr_vst), file = "../results/matrices/sorghum_RUVr_k1_vst_counts.tsv", sep = "\t" ,quote = F, col.names = T)

# ---- RUVg ----
# Define control genes (replace with your actual controls)

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

length(control_genes)  # how many controls were picked
head(control_genes)    # quick preview

set_RUVg <- RUVg(counts(dds), control_genes, k = 2)

dataRUVg <- DESeqDataSetFromMatrix(countData = set_RUVg$normalizedCounts,
                                   colData = metadata,
                                   design = ~ group)

RUVg_vst <- varianceStabilizingTransformation(dataRUVg)

# PCA for RUVg
pca_data_g <- plotPCA(RUVg_vst, intgroup = "group", returnData = TRUE)
percentVar_g <- round(100 * attr(pca_data_g, "percentVar"))

p_g <- ggplot(pca_data_g, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar_g[1], "%")) +
  ylab(paste0("PC2: ", percentVar_g[2], "%")) +
  scale_colour_manual(values = colors) +
  theme_bw(base_size=22)

ggsave(p_g, filename = "RUVg_groups.png", units = "cm", width = 15*1.3, height = 15, dpi = 320)

write.table(assay(dataRUVg), file = "../results/matrices/sorghum_RUVg_k2_raw_counts.tsv", sep = "\t", quote = F)
write.table(assay(RUVg_vst), file = "../results/matrices/sorghum_RUVg_k2_vst_counts.tsv", sep = "\t", quote = F)

# ---- RUVs ----
replicates <- makeGroups(dds$Group1)
set_RUVs <- RUVs(counts(dds), rownames(dds), replicates, k = 1)

dataRUVs <- DESeqDataSetFromMatrix(countData = set_RUVs$normalizedCounts,
                                   colData = metadata,
                                   design = ~ group)

RUVs_vst <- varianceStabilizingTransformation(dataRUVs)

# PCA for RUVs
pca_data_s <- plotPCA(RUVs_vst, intgroup = "group", returnData = TRUE)
percentVar_s <- round(100 * attr(pca_data_s, "percentVar"))

p_s <- ggplot(pca_data_s, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar_s[1], "%")) +
  ylab(paste0("PC2: ", percentVar_s[2], "%")) +
  scale_colour_manual(values = colors) +
  theme_bw(base_size=22)

ggsave(p_s, filename = "k3_RUVs_groups.png", units = "cm", width = 15*1.3, height = 15, dpi = 320)

write.table(assay(dataRUVs), file = "../results/matrices/sorghum_RUVs_k2_raw_counts.tsv", sep = "\t", quote = F)
write.table(assay(RUVs_vst), file = "../results/matrices/sorghum_RUVs_k2_vst_counts.tsv", sep = "\t", quote = F)



# Make Differential expression analysis 
# load files paths

sample_files = paste0("/home/diegoj/rnaseq_diatraea/rnaseq/run_paired_samples/star_salmon/", pull(metadata , "sample"), "/quant.sf")
# name table columns
names(sample_files) = pull(metadata, "sample")
# relate genes to transcripts
tx2gene = read.table("/home/diegoj/rnaseq_diatraea/rnaseq/run_paired_samples/star_salmon/tx2gene.tsv", header = T)
# GENE MODE
# import count data to tximport
count_data = tximport( files = sample_files,
	type = "salmon",
        tx2gene =  tx2gene,
        ignoreTxVersion = F,
        ignoreAfterBar = T)

raw <- DESeqDataSetFromTximport(txi = count_data,
	colData = metadata,
        design = ~ group)


dim(raw)
temp <- as.data.frame(counts(raw))
logic <-(apply(temp,c(1,2), function(x){x>0}))
filter_genes <- rowSums(logic)>3
fi <- raw[filter_genes,]
dim(fi)
temp <- NULL
data <- estimateSizeFactors(fi)
vst <- varianceStabilizingTransformation(data)

###################################################
########## Use RUVs k=1 ###########################
###################################################

W <- set_RUVs$W
stopifnot(rownames(W) == colnames(fi))  # or reorder W appropriately
colData(fi)$W_1 <- W[,1]
design(fi) <- ~ W_1 + group


### Differencial expression analyses
dea <- DESeq(fi, parallel = T)
dea_contrast <- results(dea, lfcThreshold= 1, altHypothesis="greaterAbs", parallel = T, contrast = c("group", "control", "infected"))
#save(dea_contrast, file = paste0("./DEA/by_groups/", group1,"_vs_", group2, ".Rdata", sep = ""))

dea_df <- as.data.frame(dea_contrast)
##  Save files for trinotate

baseMeanA <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$group == "control"])
baseMeanB <- rowMeans(counts(dea, normalized=TRUE)[,colData(dea)$group == "infected"])

res = cbind(baseMeanA, baseMeanB, dea_df)
res = cbind(sampleA="control", sampleB="infected", as.data.frame(res))
res = res[complete.cases(res),]

final <- res %>% filter(abs(log2FoldChange) > 1 & padj < 0.05)
final <- final[rev(order(final$log2FoldChange)),]

up <- final %>% filter(log2FoldChange > 0)
down <- final %>% filter(log2FoldChange < 0)


write.table(up, "up_regulated.csv", sep = ",", quote = F) 
write.table(down, "down_regulated.csv", sep = ",", quote = F)

#write.table(res, file = paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )
#write.table(res, file = paste0(outdir,"Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.DE_results"), sep = "\t", quote = F )
#write.table(cbind(as.character(filter$Group), filter$Sample_file), file= paste0(outdir, group1, "_vs_", group2,"/Trinity.isoform.counts.matrix.", group1, "_vs_", group2, ".DESeq.samples"), col.names = F, row.names = F, sep = "\t", quote = F)



# GO enrichment 
library(topGO)
library(ggplot2)
library(dplyr)
library(parallel)
#library(clusterProfiler)
numCores <- 150

GO <- read.table("/home/diegoj/rnaseq_diatraea/panzzer/annot_01/gene_go_annotations.txt", header=FALSE, stringsAsFactors=FALSE)

# Rename columns if necessary (assuming your data has two columns)
colnames(GO) <- c("Gene", "GO_term")

# Group by Gene and aggregate GO terms into a single string separated by spaces
formatted_GO <- aggregate(GO_term ~ Gene, GO, function(x) paste(x, collapse=" "))
gene2GO <- strsplit(formatted_GO$GO_term , " ")
names(gene2GO) <- formatted_GO$Gene
geneNames <- names(gene2GO)


# select set of genes to make overrepresentation test
MyInterestingGenes <- rownames(up)
geneList <- factor(as.integer(geneNames %in% MyInterestingGenes))
names(geneList) <- geneNames
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)
allGO=usedGO(GOdata)
# MIRAR SI DA VALOR Pcorregido el TOP GO
# classic ingnora la topologia del go
# revisar algoritmos

Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultsWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# Make results  table
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <- cbind(table1,p.adj)
all_res_final <- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]
# Save first top 50 ontolgies sorted by adjusted pvalues
write.table(results.table.bh, file = "GO_up.csv", quote=FALSE, row.names=FALSE, sep = ",")

ntop <- 20

ggdata <- results.table.bh[1:ntop,]
ggdata <- ggdata[complete.cases(ggdata), ]
#ggdata <- all_res_final
#aux <- go2term(ggdata$GO.ID)
#colnames(aux) <- c("GO.ID", "Lterm")

#ggdata <- merge(ggdata, aux, by = "GO.ID")
ggdata$p.adj <- as.numeric(ggdata$p.adj)

ggdata <- ggdata[order(ggdata$p.adj),]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order

#gg1 <- ggplot(ggdata, aes(x = Lterm, y = -log10(p.adj) ))+
#  geom_point(size = 6, colour = "black") +
#  scale_size(range = c(2.5,12.5)) +
#  xlab('GO Term') +
#  ylab('-log(p)') +
#  labs(title = 'GO Biological processes')+
#  theme_bw(base_size = 24) +
#  coord_flip()

gg1 <- ggplot(ggdata, aes(x = Term, y = -log10(p.adj), size = Significant)) +
  geom_point(colour = "black") +
  scale_size(range = c(2.5, 12.5)) +
  xlab('GO Term') +
  ylab('-log(p)') +
  labs(title = 'GO Biological processes', size = 'Significant') +
  theme_bw(base_size = 24) +
  coord_flip()

ggsave("GO_up.svg", device = "svg", width = 40, height = 30, dpi = 300, units = "cm")

