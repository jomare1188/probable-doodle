# =============================================================================
# KEGG enrichment + GO-KEGG gene networks for stem t3_vs_t4 (no RUV)
# Same methodology as t3_vs_t4/star_salmon/deseq2_qc/kegg.r (leaf)
# Run from stem_t3_vs_t4/star_salmon/deseq2_qc after dea_go.r
# =============================================================================

library(readr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)
library(svglite)

base_dir <- "/dados04/jorge/rnaseq_diatraea"
res_dir  <- "."

annot_file <- file.path(base_dir, "reference_genomes/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.annotation_info.txt")

# =============================================================================
# 1. DATA LOADING
# =============================================================================

# gene (locusName) - KO table from the Phytozome annotation
annot_df <- read.delim(annot_file)
kos_genes <- annot_df %>%
  select(gene = locusName, KEGG_ko = KO) %>%
  filter(KEGG_ko != "" & !is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = "\\s+") %>%
  distinct()

# gene - GO table (for network edges)
all_go <- annot_df %>%
  select(gene = locusName, GO = GO) %>%
  filter(GO != "" & !is.na(GO)) %>%
  separate_rows(GO, sep = "\\s+") %>%
  distinct()

universe_kos <- unique(kos_genes$KEGG_ko)

# DEA tables, gene ids without the .v2.1 suffix to match the annotation
read_dea <- function(file) {
  read.table(file, header = TRUE, sep = ",") %>%
    tibble::rownames_to_column("gene") %>%
    mutate(gene = sub("\\.v2\\.1$", "", gene))
}

dea_up   <- read_dea(file.path(res_dir, "up_regulated.csv"))
dea_down <- read_dea(file.path(res_dir, "down_regulated.csv"))

go_up   <- read.table(file.path(res_dir, "GO_up.csv"),   header = TRUE, sep = ",")
go_down <- read.table(file.path(res_dir, "GO_down.csv"), header = TRUE, sep = ",")

# =============================================================================
# 2. KEGG ENRICHMENT
# =============================================================================

perform_kegg_enrichment <- function(gene_kos, universe_kos,
                                    pvalue_cutoff = 0.05,
                                    qvalue_cutoff = 0.05) {
  enrichKEGG(
    gene = gene_kos,
    universe = universe_kos,
    organism = "ko",
    pAdjustMethod = "fdr",
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )
}

run_kegg <- function(dea, prefix, title) {
  gene_kos <- kos_genes %>% filter(gene %in% dea$gene) %>% pull(KEGG_ko)
  cat(prefix, ": ", length(unique(dea$gene[dea$gene %in% kos_genes$gene])), " of ",
      nrow(dea), " DE genes have a KO (", length(unique(gene_kos)), " KOs)\n", sep = "")

  ekegg <- if (length(gene_kos) > 0) perform_kegg_enrichment(gene_kos, universe_kos) else NULL
  sig <- if (is.null(ekegg)) data.frame() else ekegg@result %>% filter(p.adjust < 0.05)
  cat(prefix, ": ", nrow(sig), " significant KEGG pathways\n", sep = "")
  write.table(sig, file.path(res_dir, paste0(prefix, ".csv")), sep = ",", quote = TRUE, row.names = FALSE)

  if (nrow(sig) > 0) {
    p <- dotplot(ekegg, title = title)
    ggsave(file.path(res_dir, paste0(prefix, ".svg")), p, width = 8, height = 6, bg = "white")
    ggsave(file.path(res_dir, paste0(prefix, ".png")), p, width = 8, height = 6, dpi = 300, bg = "white")
    ggsave(file.path(res_dir, paste0(prefix, ".pdf")), p, width = 8, height = 6)
  }
  ekegg
}

kegg_up   <- run_kegg(dea_up,   "kegg_up",   "KEGG Enrichment - Up-regulated")
kegg_down <- run_kegg(dea_down, "kegg_down", "KEGG Enrichment - Down-regulated")

# =============================================================================
# 3. NETWORK PREPARATION
# =============================================================================

prepare_kegg_edges <- function(kegg_result, dea_data, padj_cutoff = 0.05) {
  if (is.null(kegg_result)) return(data.frame(gene = character(), Description = character()))
  kegg_result@result %>%
    filter(p.adjust < padj_cutoff) %>%
    select(Description, geneID) %>%
    mutate(geneID = strsplit(geneID, "/")) %>%
    unnest(geneID) %>%
    mutate(geneID = trimws(geneID)) %>%
    left_join(kos_genes, by = c("geneID" = "KEGG_ko"), relationship = "many-to-many") %>%
    filter(gene %in% dea_data$gene) %>%
    select(gene, Description)
}

prepare_go_edges <- function(go_results, dea_data, top_n = 20) {
  if (nrow(go_results) == 0) return(data.frame(gene = character(), Description = character()))
  go_results %>%
    arrange(p.adj) %>%
    head(top_n) %>%
    left_join(all_go, by = c("GO.ID" = "GO"), relationship = "many-to-many") %>%
    filter(gene %in% dea_data$gene) %>%
    select(gene, Description = Term)
}

combine_edge_lists <- function(go_edges, kegg_edges, dea_data) {
  rbind(go_edges %>% mutate(Class = "GO"),
        kegg_edges %>% mutate(Class = "KEGG")) %>%
    filter(!is.na(gene) & gene != "NA") %>%
    distinct() %>%
    left_join(dea_data %>% select(gene, log2FoldChange), by = "gene")
}

# =============================================================================
# 4. NETWORK VISUALIZATION
# =============================================================================

create_network_graph <- function(edge_list, tf_genes = NULL) {
  g <- graph_from_data_frame(d = edge_list, directed = FALSE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

  V(g)$type <- ifelse(V(g)$name %in% edge_list$Description, "Description", "Gene")
  if (!is.null(tf_genes)) {
    V(g)$type[V(g)$name %in% tf_genes] <- "Transcription Factor"
  }

  desc_class <- edge_list %>% distinct(Description, Class)
  V(g)$Class <- desc_class$Class[match(V(g)$name, desc_class$Description)]

  V(g)$Category <- case_when(
    V(g)$type == "Transcription Factor" ~ "Transcription Factor",
    V(g)$type == "Gene" ~ "Gene",
    V(g)$Class == "GO" ~ "GO term",
    V(g)$Class == "KEGG" ~ "KEGG pathway",
    TRUE ~ NA_character_
  )

  V(g)$label <- ifelse(V(g)$type == "Gene", "", V(g)$name)
  V(g)$label.cex <- rescale(degree(g), to = c(0.6, 0.8))
  g
}

plot_gene_network <- function(graph, title, layout = "tree") {
  color_palette <- c(
    "Gene" = "gray80",
    "GO term" = "#66c2a5",
    "KEGG pathway" = "#fc8d62",
    "Transcription Factor" = "#FF0033"
  )
  colors_to_use <- color_palette[names(color_palette) %in% unique(V(graph)$Category)]

  ggraph(graph, layout = layout) +
    geom_edge_link(alpha = 0.4, colour = "grey70") +
    geom_node_point(aes(color = Category), size = 4, show.legend = TRUE) +
    geom_node_text(aes(label = label), repel = TRUE,
                   angle = 60, size = 3, color = "black") +
    scale_color_manual(name = "Node Type", values = colors_to_use) +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
}

save_network_outputs <- function(plot, graph, prefix,
                                 width = 16, height = 9, dpi = 300) {
  ggsave(paste0(prefix, ".svg"), plot, width = width, height = height, bg = "white")
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height)

  edges <- igraph::as_data_frame(graph, what = "edges") %>%
    mutate(
      source_Class = "gene",
      target_Class = V(graph)$Class[match(to, V(graph)$name)]
    )
  write.table(edges, paste0(prefix, "_edges.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

build_network <- function(kegg_result, go_results, dea_data, direction, title) {
  combined <- combine_edge_lists(prepare_go_edges(go_results, dea_data),
                                 prepare_kegg_edges(kegg_result, dea_data),
                                 dea_data)
  cat(direction, " network: ", nrow(combined), " gene-term edges\n", sep = "")
  if (nrow(combined) == 0) return(invisible(NULL))
  g <- create_network_graph(combined)
  save_network_outputs(plot_gene_network(g, title), g,
                       file.path(res_dir, paste0("gene_network_", direction)))
}

build_network(kegg_up,   go_up,   dea_up,   "up",   "Up-Regulated Gene–Functional Term Network")
build_network(kegg_down, go_down, dea_down, "down", "Down-Regulated Gene–Functional Term Network")
