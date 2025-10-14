
library(readr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(clusterProfiler)

emapper_df <- read_tsv("/home/j/proteins.emapper.emapper.annotations", comment = "#", 
                       guess_max = 100000)  
                            
gene2protein <- read.table("/home/j/gene2protein.txt", col.names = c("gene", "protein"))

emapper_with_genes <- emapper_df %>%
  left_join(gene2protein, by = c("query" = "protein"))


# Create a table with IDs and ko:Knumber from input file 1
KOs <- emapper_with_genes %>% select(gene,KEGG_ko)

KOs_long <- KOs %>%  filter(KEGG_ko != "-" & KEGG_ko != "" & !is.na(KEGG_ko)) %>% separate_rows(KEGG_ko, sep = ",")
colnames(KOs_long) <- c("gene", "KEGG_ko")


KOs_long_formatted <- KOs_long %>%
  mutate(KEGG_ko = sub("^ko:", "", KEGG_ko))

DOWN <- rownames(read.table("/home/j/down_regulated.csv", header = T, sep = ","))
UP   <- rownames(read.table("/home/j/up_regulated.csv", header = T, sep = ","))

library(dplyr)

####
# Get KOs for UP genes
up_gene_kos <- gene2protein %>%
  filter(gene %in% UP) %>%
  left_join(KOs_long, by = "gene") %>%
  distinct(gene, KEGG_ko) %>%
  filter(!is.na(KEGG_ko))

down_gene_kos <- gene2protein %>%
  filter(gene %in% DOWN) %>%
  left_join(KOs_long, by = "gene") %>%
  distinct(gene, KEGG_ko) %>%
  filter(!is.na(KEGG_ko))

up_kos <- up_gene_kos$KEGG_ko
up_kos <- sub("^ko:", "", up_kos)

down_kos <- down_gene_kos$KEGG_ko
down_kos <- sub("^ko:", "", down_kos)

universe_kos <- KOs_long$KEGG_ko
universe_kos <- sub("^ko:", "", universe_kos)

universe_genes <- KOs_long$gene


ekegg_up <- enrichKEGG(
  gene = up_kos,
  universe = universe_kos ,
  organism = "ko",
  pAdjustMethod = "fdr",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)


ekegg_down <- enrichKEGG(
  gene = down_kos,
  universe = universe_kos,
  organism = "ko",
  pAdjustMethod = "fdr",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

dotplot(ekegg_down)
dotplot(ekegg_up)




kegg_up_df <- ekegg_up@result
kos_up_enriched <- kegg_up_df %>% filter(p.adjust < 0.05) %>% select(Description, geneID)

# Split the geneID column by "/"
edge_list <- kos_up_enriched %>%
  mutate(geneID = strsplit(geneID, "/")) %>%  # split into list
  unnest(geneID) %>%                          # expand into rows
  mutate(geneID = trimws(geneID))  

# reformat KOS 
KOs_long$KEGG_ko  <- sub("^ko:", "", KOs_long$KEGG_ko)
 
# Join edge list with gene-to-KO table
edge_list_genes <- edge_list %>%  left_join(KOs_long, by = c("geneID" = "KEGG_ko"))

# Optional: reorder columns for clarity
edge_list_genes <- edge_list_genes %>%  select(gene, Description, geneID)

# Inspect result
head(edge_list_genes)

# add logfold2change

dea_up <- read.table("/home/j/up_regulated.csv", sep = ",", header = T)
dea_up <- dea_up %>%  tibble::rownames_to_column("gene")

KEGG_edge_list_with_lfc <- edge_list_genes %>% left_join(dea_up %>% select(gene, log2FoldChange), by = "gene")

KEGG_edge_list_with_lfc <- KEGG_edge_list_with_lfc[!is.na(KEGG_edge_list_with_lfc$log2FoldChange), ]




### add TOP GO RESULTS 

all_GO <- read.table("/home/j/gene_go_annotations.txt", col.names = c("gene", "GO"))

up_GO <- read.table("/home/j/GO_up.csv", header = T, sep = ",")

up_GO_genes <- up_GO %>%
  left_join(all_GO, by = c("GO.ID" = "GO")) %>%          # bring the gene column
  group_by(GO.ID, Term, Annotated, Significant, Expected, Classic, p.adj) %>% 
  summarise(gene = paste(unique(gene), collapse = "/"), .groups = "drop")  # collapse genes

up_GO_genes <- head(up_GO_genes[order(up_GO_genes$p.adj ,decreasing = F),], 20)

go_edge_list <- up_GO_genes %>% separate_rows(gene, sep = "/") %>% select("gene", "Term") 

# add logfold changue

dea_up <- read.table("/home/j/up_regulated.csv", sep = ",", header = T)
dea_up <- dea_up %>%  tibble::rownames_to_column("gene")

go_edge_list_with_lfc <- go_edge_list %>% left_join(dea_up %>% select(gene, log2FoldChange), by = "gene")

go_edge_list_with_lfc <- go_edge_list_with_lfc[!is.na(go_edge_list_with_lfc$log2FoldChange), ]



# make big combined edge list

go_edge_list <- go_edge_list_with_lfc %>% select(gene, Term) %>% mutate(Class = "GO")

edge_list_KEGG <- KEGG_edge_list_with_lfc %>% select(gene, Description) %>% mutate(Class = "KEGG")



colnames(go_edge_list) <- colnames(edge_list_KEGG)
combined_edge_list  <- rbind(go_edge_list, edge_list_KEGG)

combined_edge_list <- combined_edge_list %>%  filter(!is.na(gene) & gene != "NA") 

combined_edge_list_with_fc <- combined_edge_list %>% left_join(dea_up %>% select(gene, log2FoldChange), by = "gene")## Import to igraph

library(igraph)

# ---- 1. Build the graph ----
g <- graph_from_data_frame(d = combined_edge_list_with_fc, directed = FALSE, )
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# ---- 2. Set node attributes ----
# Get a vector indicating which nodes are GO/KEGG terms
V(g)$type <- ifelse(V(g)$name %in% combined_edge_list$Description, "Description", "Gene")

# Assign the Class (GO/KEGG) only to description nodes
desc_class <- combined_edge_list %>% distinct(Description, Class)

# Match the class info to vertex names
V(g)$Class <- desc_class$Class[match(V(g)$name, desc_class$Description)]

# ---- 3. Define colors ----
# We'll color nodes by Class only for Description nodes
V(g)$color <- ifelse(
  V(g)$type == "Gene", "gray80",              # genes: gray
  ifelse(V(g)$Class == "GO", "#66c2a5", "#fc8d62")  # GO=greenish, KEGG=orangeish
)

# ---- Plot ----

library(scales)

###########
# Calculate degrees
degrees <- degree(g)

# Set label sizes proportional to degree
V(g)$label.cex <- scales::rescale(degrees, to = c(0.6, 0.8))  # Adjust range as needed

# Set labels: empty string for Genes, keep name for GO and KEGG
V(g)$label <- ifelse(V(g)$Class == "Gene", "", V(g)$name)

###############################################3

library(ggraph)
library(ggrepel)
library(dplyr)

# Ensure proper class labels for all vertices
V(g)$Category <- ifelse(
  V(g)$type == "Gene", "Gene",
  ifelse(V(g)$Class == "GO", "GO term", "KEGG pathway")
)

library(ggraph)
library(ggrepel)
library(dplyr)

# Make sure the category column exists and is a factor
V(g)$Category <- ifelse(
  V(g)$type == "Gene", "Gene",
  ifelse(V(g)$Class == "GO", "GO term", "KEGG pathway")
)

# Convert to factor to force ggplot to create a legend
V(g)$Category <- factor(V(g)$Category, levels = c("Gene", "GO term", "KEGG pathway"))

# Build the plot
net_plot <- ggraph(g, layout = 'tree') +
  geom_edge_link(alpha = 0.4, colour = "grey70") +
  geom_node_point(aes(color = Category), size = 4, show.legend = TRUE) +
  geom_node_text(aes(label = label), repel = T, angle = 60, size = 3, color = "black") +
  scale_color_manual(
    name = "Node Type",
    values = c("Gene" = "gray80", "GO term" = "#66c2a5", "KEGG pathway" = "#fc8d62")
  ) +
  theme_void() +
  ggtitle("Up Regulated Gene–Functional Term Network") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # or "bottom" or "top"
  )


# Set a standard size (adjust for your network complexity)
plot_width <- 8*2# in inches
plot_height <- 6*1.5   # in inches
plot_dpi <- 300    # for PNG

# SVG
ggsave("gene_network_up.svg", net_plot, width = plot_width, height = plot_height, bg = "white")
# PNG (high-res for publications)
ggsave("gene_network_up.png", net_plot, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
# PDF
ggsave("gene_network_up.pdf", net_plot, width = plot_width, height = plot_height)

#####################################################################################################


edges <- igraph::as_data_frame(g, what = "edges")

# Add Class of the target node
edges <- edges %>%
  mutate(
    source_Class = "gene",
    target_Class = V(g)$Class[match(to, V(g)$name)]
  )

# Save to TSV
write.table(edges, "up_network_edges_with_class.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)




############ Now with Down regulated genes ########################3
# get kos significatly enriched 

kegg_down_df <- ekegg_down@result

kos_down_enriched <- kegg_down_df %>% filter(p.adjust < 0.05) %>% select(Description, geneID)


# Split the geneID column by "/"
edge_list <- kos_down_enriched %>%
  mutate(geneID = strsplit(geneID, "/")) %>%  # split into list
  unnest(geneID) %>%                          # expand into rows
  mutate(geneID = trimws(geneID))  

# reformat KOS 
KOs_long$KEGG_ko  <- sub("^ko:", "", KOs_long$KEGG_ko)

# Join edge list with gene-to-KO table
edge_list_genes <- edge_list %>%  left_join(KOs_long, by = c("geneID" = "KEGG_ko"))

# Optional: reorder columns for clarity
edge_list_genes <- edge_list_genes %>%  select(gene, Description, geneID)
###


###
# Inspect result
head(edge_list_genes)

# add logfold2change

dea_down <- read.table("/home/j/down_regulated.csv", sep = ",", header = T)
dea_down <- dea_down %>%  tibble::rownames_to_column("gene")

KEGG_edge_list_with_lfc <- edge_list_genes %>% left_join(dea_down %>% select(gene, log2FoldChange), by = "gene")

KEGG_edge_list_with_lfc <- KEGG_edge_list_with_lfc[!is.na(KEGG_edge_list_with_lfc$log2FoldChange), ]


### add TOP GO RESULTS 

all_GO <- read.table("/home/j/gene_go_annotations.txt", col.names = c("gene", "GO"))

down_GO <- read.table("/home/j/GO_down.csv", header = T, sep = ",")

down_GO_genes <- down_GO %>%
  left_join(all_GO, by = c("GO.ID" = "GO")) %>%          # bring the gene column
  group_by(GO.ID, Term, Annotated, Significant, Expected, Classic, p.adj) %>% 
  summarise(gene = paste(unique(gene), collapse = "/"), .groups = "drop")  # collapse genes

down_GO_genes <- head(down_GO_genes[order(down_GO_genes$p.adj ,decreasing = F),], 20)

go_edge_list <- down_GO_genes %>% separate_rows(gene, sep = "/") %>% select("gene", "Term") 

# add logfold changue

dea_down <- read.table("/home/j/down_regulated.csv", sep = ",", header = T)
dea_down <- dea_down %>%  tibble::rownames_to_column("gene")

go_edge_list_with_lfc <- go_edge_list %>% left_join(dea_down %>% select(gene, log2FoldChange), by = "gene")

go_edge_list_with_lfc <- go_edge_list_with_lfc[!is.na(go_edge_list_with_lfc$log2FoldChange), ]


# make big combined edge list

go_edge_list <- go_edge_list_with_lfc %>% select(gene, Term) %>% mutate(Class = "GO")

edge_list_KEGG <- KEGG_edge_list_with_lfc %>% select(gene, Description) %>% mutate(Class = "KEGG")



colnames(go_edge_list) <- colnames(edge_list_KEGG)
combined_edge_list  <- rbind(go_edge_list, edge_list_KEGG)

combined_edge_list <- combined_edge_list %>%  filter(!is.na(gene) & gene != "NA") 


##############

# ---- 1. Build the graph ----
g <- graph_from_data_frame(d = combined_edge_list, directed = FALSE, )
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# ---- 2. Set node attributes ----
# Get a vector indicating which nodes are GO/KEGG terms
V(g)$type <- ifelse(V(g)$name %in% combined_edge_list$Description, "Description", "Gene")

# Assign the Class (GO/KEGG) only to description nodes
desc_class <- combined_edge_list %>% distinct(Description, Class)

# Match the class info to vertex names
V(g)$Class <- desc_class$Class[match(V(g)$name, desc_class$Description)]

# ---- 3. Define colors ----
# We'll color nodes by Class only for Description nodes
V(g)$color <- ifelse(
  V(g)$type == "Gene", "gray80",              # genes: gray
  ifelse(V(g)$Class == "GO", "#66c2a5", "#fc8d62")  # GO=greenish, KEGG=orangeish
)

# ---- Plot ----

library(scales)

###########
# Calculate degrees
degrees <- degree(g)

# Set label sizes proportional to degree
V(g)$label.cex <- scales::rescale(degrees, to = c(0.6, 0.8))  # Adjust range as needed

# Set labels: empty string for Genes, keep name for GO and KEGG
V(g)$label <- ifelse(V(g)$Class == "Gene", "", V(g)$name)

###############################################3

library(ggraph)
library(ggrepel)
library(dplyr)

# Ensure proper class labels for all vertices
V(g)$Category <- ifelse(
  V(g)$type == "Gene", "Gene",
  ifelse(V(g)$Class == "GO", "GO term", "KEGG pathway")
)

# Convert to factor to force ggplot to create a legend
V(g)$Category <- factor(V(g)$Category, levels = c("Gene", "GO term", "KEGG pathway"))

# Build the plot
net_plot <- ggraph(g, layout = 'tree') +
  geom_edge_link(alpha = 0.4, colour = "grey70") +
  geom_node_point(aes(color = Category), size = 4, show.legend = TRUE) +
  geom_node_text(aes(label = label), repel = T, angle = 60, size = 3, color = "black") +
  scale_color_manual(
    name = "Node Type",
    values = c("Gene" = "gray80", "GO term" = "#66c2a5", "KEGG pathway" = "#fc8d62")
  ) +
  theme_void() +
  ggtitle("Up Regulated Gene–Functional Term Network") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"  # or "bottom" or "top"
  )


# Set a standard size (adjust for your network complexity)
plot_width <- 8*2# in inches
plot_height <- 6*1.5   # in inches
plot_dpi <- 300    # for PNG

# SVG
ggsave("gene_network_down.svg", net_plot, width = plot_width, height = plot_height, bg = "white")
# PNG (high-res for publications)
ggsave("gene_network_down.png", net_plot, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
# PDF
ggsave("gene_network_down.pdf", net_plot, width = plot_width, height = plot_height)

#####################################################################################################


edges <- igraph::as_data_frame(g, what = "edges")

# Add Class of the target node
edges <- edges %>%
  mutate(
    source_Class = "gene",
    target_Class = V(g)$Class[match(to, V(g)$name)]
  )

# Save to TSV
write.table(edges, "down_network_edges_with_class.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)




#############

























## Import to igraph

library(igraph)

# ---- 1. Build the graph ----
g <- graph_from_data_frame(d = combined_edge_list, directed = FALSE, )
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

# ---- 2. Set node attributes ----
# Get a vector indicating which nodes are GO/KEGG terms
V(g)$type <- ifelse(V(g)$name %in% combined_edge_list$Description, "Description", "Gene")

# Assign the Class (GO/KEGG) only to description nodes
desc_class <- combined_edge_list %>%
  distinct(Description, Class)

# Match the class info to vertex names
V(g)$Class <- desc_class$Class[match(V(g)$name, desc_class$Description)]

# ---- 3. Define colors ----
# We'll color nodes by Class only for Description nodes
V(g)$color <- ifelse(
  V(g)$type == "Gene", "gray80",              # genes: gray
  ifelse(V(g)$Class == "GO", "#66c2a5", "#fc8d62")  # GO=greenish, KEGG=orangeish
)

# ---- 4. Define sizes ----
V(g)$size <- ifelse(V(g)$type == "Gene", 3, 6)

# ---- Layout options ----

# 2. Or alternatively:
layout_bip <- layout_with_graphopt(g, niter = 800)

# ---- Plot ----


