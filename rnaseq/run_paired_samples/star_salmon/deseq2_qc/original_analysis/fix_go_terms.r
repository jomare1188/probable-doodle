# Restore full GO term names in GO_up.csv / GO_down.csv and redraw the GO plots.
# topGO's GenTable truncates terms to 40 chars by default, which made different
# terms share a label (and collapse into one row of the plot). Statistics are
# untouched: only the Term column is replaced using GO.db, looked up by GO.ID.

library(GO.db)
library(ggplot2)

fix_go <- function(prefix, ntop = 20) {
  go <- read.table(paste0(prefix, ".csv"), header = TRUE, sep = ",", quote = "\"",
                   stringsAsFactors = FALSE, comment.char = "")
  full <- unname(Term(GOTERM[go$GO.ID]))

  # sanity check: the old (possibly truncated) label must be a prefix of the full name
  # (commas ignored: the old files had them stripped)
  old <- sub("\\.\\.\\.$", "", go$Term)
  stopifnot(startsWith(gsub(",", "", full), gsub(",", "", old)))
  cat(prefix, ": ", sum(go$Term != full), " of ", nrow(go), " terms restored\n", sep = "")
  go$Term <- full

  write.table(go, paste0(prefix, ".csv"), quote = TRUE, row.names = FALSE, sep = ",")

  # same plot as ruv.r
  ggdata <- go[1:min(ntop, nrow(go)),]
  ggdata <- ggdata[complete.cases(ggdata), ]
  ggdata$p.adj <- as.numeric(ggdata$p.adj)
  ggdata <- ggdata[order(ggdata$p.adj),]
  ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order

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

fix_go("GO_up")
fix_go("GO_down")
