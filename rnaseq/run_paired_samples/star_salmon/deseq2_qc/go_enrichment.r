# topGO over-representation (BP, classic Fisher) for up- and down-regulated genes.
# Same method as run_paired_samples/star_salmon/deseq2_qc/ruv.r (GO section) with the
# full term names of fix_go_terms.r, run for both gene sets, and numeric p-values
# for the p < 0.05 filter (see run_go below).
#
# Usage: Rscript go_enrichment.r <deg_dir> <gene_go_annotations.txt> <outdir>
#   defaults: .  ../../../../panzzer/annot_01/gene_go_annotations.txt  .
# Env: conda R_popstat_jorge (topGO 2.58, same as the strict analysis)

library(topGO)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
deg_dir <- if (length(args) >= 1) args[1] else "."
go_file <- if (length(args) >= 2) args[2] else "../../../../panzzer/annot_01/gene_go_annotations.txt"
outdir  <- if (length(args) >= 3) args[3] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

GO <- read.table(go_file, header = FALSE, stringsAsFactors = FALSE)
colnames(GO) <- c("Gene", "GO_term")
formatted_GO <- aggregate(GO_term ~ Gene, GO, function(x) paste(x, collapse = " "))
gene2GO <- strsplit(formatted_GO$GO_term, " ")
names(gene2GO) <- formatted_GO$Gene
geneNames <- names(gene2GO)

run_go <- function(set, ntop = 20) {
  degs <- read.table(file.path(deg_dir, paste0(set, "_regulated.csv")), header = TRUE, sep = ",")
  geneList <- factor(as.integer(geneNames %in% rownames(degs)))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO)
  allGO <- usedGO(GOdata)
  Classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic', numChar = 1000)
  # GenTable returns p-values as text ("3.3e-10", "< 1e-30"); a text comparison with 0.05
  # silently drops the most significant terms (as happened in the first strict analysis),
  # so use topGO's numeric p-values instead
  table$Classic <- score(Classic)[table$GO.ID]

  table1 <- filter(table, Classic < 0.05)
  # 4 significant digits (rounding to 4 decimals turns very small p.adj into 0 = Inf in the plot)
  p.adj <- signif(p.adjust(table1$Classic, method = "BH"), digits = 4)
  all_res_final <- cbind(table1, p.adj)
  all_res_final <- all_res_final[order(all_res_final$p.adj), ]
  results.table.bh <- all_res_final[which(all_res_final$p.adj <= 0.05), ]
  write.table(results.table.bh, file = file.path(outdir, paste0("GO_", set, ".csv")),
              quote = TRUE, row.names = FALSE, sep = ",")
  message("GO_", set, ": ", nrow(results.table.bh), " terms with p.adj <= 0.05")

  # same plot as fix_go_terms.r
  ggdata <- results.table.bh[1:min(ntop, nrow(results.table.bh)), ]
  ggdata <- ggdata[complete.cases(ggdata), ]
  if (nrow(ggdata) == 0) return(invisible(NULL))
  ggdata$p.adj <- as.numeric(ggdata$p.adj)
  ggdata <- ggdata[order(ggdata$p.adj), ]
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term))

  gg1 <- ggplot(ggdata, aes(x = Term, y = -log10(p.adj), size = Significant)) +
    geom_point(colour = "black") +
    scale_size(range = c(2.5, 12.5)) +
    xlab('GO Term') +
    ylab('-log(p)') +
    labs(title = 'GO Biological processes', size = 'Significant') +
    theme_bw(base_size = 24) +
    coord_flip()

  ggsave(file.path(outdir, paste0("GO_", set, ".svg")), plot = gg1, device = "svg", width = 40, height = 30, dpi = 300, units = "cm")
  ggsave(file.path(outdir, paste0("GO_", set, ".pdf")), plot = gg1, device = "pdf", width = 40, height = 30, dpi = 300, units = "cm")
}

run_go("up")
run_go("down")
