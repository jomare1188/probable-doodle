# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(topGO)
library(dplyr)
library(ggplot2)
library(scales)
library(openxlsx)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
go_raw <- read.table("insert corresponding file", header=FALSE, stringsAsFactors=FALSE)

# Insert file!!!
up_regulated <- read.csv("insert corresponding file", header=TRUE, stringsAsFactors=FALSE)

# Insert file!!!
down_regulated <- read.csv("insert corresponding file", header=TRUE, stringsAsFactors=FALSE)

colnames(go_raw) <- c("Gene", "GO_term")

# ==============================================================================
# Step3_Data_Preprocessing
# ==============================================================================
formatted_go <- aggregate(GO_term ~ Gene, go_raw, function(x) paste(x, collapse=" "))
gene2go <- strsplit(formatted_go$GO_term , " ")
names(gene2go) <- formatted_go$Gene
gene_names <- names(gene2go)

# ==============================================================================
# Step4_Enrichment_Analysis_Function
# ==============================================================================
process_go_enrichment <- function(gene_set, direction = "UP") {
  col_1_ids <- as.character(trimws(gene_set[,1]))
  row_ids <- as.character(trimws(rownames(gene_set)))
  
  if (sum(col_1_ids %in% gene_names) > 0) {
    my_interesting_genes <- col_1_ids
  } else if (sum(row_ids %in% gene_names) > 0) {
    my_interesting_genes <- row_ids
  } else {
    stop(paste("Error: IDs do not match for", direction))
  }
  
  gene_presence <- as.integer(gene_names %in% my_interesting_genes)
  gene_list <- factor(gene_presence, levels = c(0, 1))
  names(gene_list) <- gene_names
  
  go_data <- new("topGOdata", ontology = "BP", allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = gene2go)
  
  all_go <- usedGO(go_data)
  classic_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  
  res_table <- GenTable(go_data, Classic = classic_test, topNodes = length(all_go), orderBy = 'Classic')
  res_table$Classic <- as.numeric(gsub("< ", "", res_table$Classic))
  
  res_table$p.adj <- p.adjust(res_table$Classic, method = "BH")
  
  filtered_table <- res_table %>% dplyr::filter(p.adj <= 0.05)
  results_table_bh <- filtered_table %>% arrange(p.adj)
  
  if(nrow(results_table_bh) > 0) results_table_bh$Direction <- direction
  
  return(results_table_bh)
}

# ==============================================================================
# Step5_Execute_Enrichment_Analysis
# ==============================================================================
results_up <- process_go_enrichment(up_regulated, "UP")
results_down <- process_go_enrichment(down_regulated, "DOWN")

all_results <- rbind(results_up, results_down)
write.csv(all_results, "all_go_terms.csv", row.names = FALSE)

# ==============================================================================
# Step6_Load_Curated_Terms
# ==============================================================================
# Insert file!!!
cleaned_data <- read.csv("insert corresponding file", header=TRUE, stringsAsFactors=FALSE)

# ==============================================================================
# Step7_Curated_Scale_Setup
# ==============================================================================
ggdata_up <- cleaned_data %>% filter(Direction == "UP")
ggdata_down <- cleaned_data %>% filter(Direction == "DOWN")

min_p_limit <- 1e-5 

if(nrow(ggdata_up) > 0) ggdata_up$p.adj[ggdata_up$p.adj < min_p_limit | is.na(ggdata_up$p.adj)] <- min_p_limit
if(nrow(ggdata_down) > 0) ggdata_down$p.adj[ggdata_down$p.adj < min_p_limit | is.na(ggdata_down$p.adj)] <- min_p_limit

if(nrow(ggdata_up) > 0) ggdata_up$log_p_scaled <- -log10(ggdata_up$p.adj)
if(nrow(ggdata_down) > 0) ggdata_down$log_p_scaled <- log10(ggdata_down$p.adj) 

ggdata_combined <- rbind(ggdata_up, ggdata_down)
ggdata_combined$Term_unique <- paste0(ggdata_combined$Term, " (", ggdata_combined$Direction, ")")
ggdata_combined <- ggdata_combined[order(ggdata_combined$log_p_scaled), ]
ggdata_combined$Term_unique <- factor(ggdata_combined$Term_unique, levels = ggdata_combined$Term_unique)

max_limit <- ceiling(max(abs(ggdata_combined$log_p_scaled), na.rm = TRUE))

# ==============================================================================
# Step8_Render_Curated_Plot
# ==============================================================================
color_modern <- c("UP" = "#573C75", "DOWN" = "#00A087")

gg_curated <- ggplot(ggdata_combined, aes(x = log_p_scaled, y = Term_unique, size = Significant)) +
  geom_point(aes(fill = Direction), shape = 21, colour = "white", stroke = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.6) +
  scale_size_continuous(range = c(4, 12), name = "Gene count", breaks = pretty_breaks(n = 4)) +
  scale_fill_manual(values = color_modern, name = "Regulation") +
  scale_x_continuous(labels = abs, limits = c(-max_limit, max_limit), breaks = scales::pretty_breaks(n = 8)) +
  scale_y_discrete(expand = expansion(add = c(0.8, 1.1))) +
  labs(title = "GO Enrichment Analysis (Curated)", 
       subtitle = "Biological Processes: DOWN-regulated (left) vs UP-regulated (right)",
       x = expression(paste("-log"["10"], "(", italic("P")["adj"], ")")), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(colour = "#e7e8e8", linewidth = 0.5, linetype = "dotted"),
    panel.grid.major.x = element_line(colour = "#dddddd", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(face = "bold", size = 16, color = "#2C3E50", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "#7F8C8D", margin = margin(b = 15)),
    axis.text.y = element_text(color = "#34495E", size = 9, margin = margin(r = 5)),
    axis.text.x = element_text(color = "#34495E", size = 10, face = "bold"),
    axis.title.x = element_text(color = "#2C3E50", size = 12, face = "bold", margin = margin(t = 12)),
    axis.ticks.x = element_line(colour = "#BDC3C7", linewidth = 0.5),
    axis.ticks.length.x = unit(0.15, "cm"),
    axis.line.x = element_line(colour = "#BDC3C7", linewidth = 0.5), 
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2C3E50", size = 11),
    legend.text = element_text(color = "#34495E", size = 10),
    legend.key = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6, colour = "white")),
         size = guide_legend(override.aes = list(colour = "#7F8C8D")))

print(gg_curated)
ggsave("GO_enrichment_Curated_Final.pdf", gg_curated, width = 14.5, height = 10, dpi = 600)
ggsave("GO_enrichment_Curated_Final.png", gg_curated, width = 14.5, height = 10, dpi = 900)

# ==============================================================================
# Step9_Global_Enrichment_Function
# ==============================================================================
# Insert file!!!
go_raw_global <- read.table("insert corresponding file", header=FALSE, stringsAsFactors=FALSE)

# Insert file!!!
up_global <- read.csv("insert corresponding file", header=TRUE, stringsAsFactors=FALSE)

# Insert file!!!
down_global <- read.csv("insert corresponding file", header=TRUE, stringsAsFactors=FALSE)

colnames(go_raw_global) <- c("Gene", "GO_term")

formatted_go_global <- aggregate(GO_term ~ Gene, go_raw_global, function(x) paste(x, collapse=" "))
gene2go_global <- strsplit(formatted_go_global$GO_term , " ")
names(gene2go_global) <- formatted_go_global$Gene
gene_names_global <- names(gene2go_global)

process_go_enrichment_global <- function(gene_set, direction = "UP") {
  col_1_ids <- as.character(trimws(gene_set[,1]))
  row_ids <- as.character(trimws(rownames(gene_set)))
  
  if (sum(col_1_ids %in= gene_names_global) > 0) {
    my_interesting_genes <- col_1_ids
  } else if (sum(row_ids %in% gene_names_global) > 0) {
    my_interesting_genes <- row_ids
  } else {
    stop(paste("Error: IDs do not match for", direction))
  }
  
  gene_presence <- as.integer(gene_names_global %in% my_interesting_genes)
  gene_list <- factor(gene_presence, levels = c(0, 1))
  names(gene_list) <- gene_names_global
  
  go_data <- new("topGOdata", ontology = "BP", allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = gene2go_global)
  
  all_go <- usedGO(go_data)
  classic_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
  
  res_table <- GenTable(go_data, Classic = classic_test, topNodes = length(all_go), orderBy = 'Classic')
  res_table$Classic <- as.numeric(gsub("< ", "", res_table$Classic))
  
  res_table$p.adj <- p.adjust(res_table$Classic, method = "BH")
  
  filtered_table <- res_table %>% dplyr::filter(Classic <= 0.05)
  results_table_bh <- filtered_table %>% arrange(Classic)
  
  if(nrow(results_table_bh) > 0) results_table_bh$Direction <- direction
  
  return(results_table_bh)
}

# ==============================================================================
# Step10_Execute_Global_Enrichment
# ==============================================================================
results_up_global <- process_go_enrichment_global(up_global, "UP")
results_down_global <- process_go_enrichment_global(down_global, "DOWN")

all_results_global <- rbind(results_up_global, results_down_global)
write.csv(all_results_global, "all_go_terms.csv", row.names = FALSE)

# ==============================================================================
# Step11_Global_Scale_Setup
# ==============================================================================
ggdata_up_global <- all_results_global %>% filter(Direction == "UP")
ggdata_down_global <- all_results_global %>% filter(Direction == "DOWN")

min_p_limit <- 1e-5 

if(nrow(ggdata_up_global) > 0) ggdata_up_global$Classic[ggdata_up_global$Classic < min_p_limit | is.na(ggdata_up_global$Classic)] <- min_p_limit
if(nrow(ggdata_down_global) > 0) ggdata_down_global$Classic[ggdata_down_global$Classic < min_p_limit | is.na(ggdata_down_global$Classic)] <- min_p_limit

if(nrow(ggdata_up_global) > 0) ggdata_up_global$log_p_scaled <- -log10(ggdata_up_global$Classic)
if(nrow(ggdata_down_global) > 0) ggdata_down_global$log_p_scaled <- log10(ggdata_down_global$Classic) 

ggdata_combined_global <- rbind(ggdata_up_global, ggdata_down_global)
ggdata_combined_global$Term_unique <- paste0(ggdata_combined_global$Term, " [", ggdata_combined_global$GO.ID, "] (", ggdata_combined_global$Direction, ")")

ggdata_combined_global <- ggdata_combined_global[order(ggdata_combined_global$log_p_scaled), ]
ggdata_combined_global$Term_unique <- make.unique(as.character(ggdata_combined_global$Term_unique))
ggdata_combined_global$Term_unique <- factor(ggdata_combined_global$Term_unique, levels = unique(ggdata_combined_global$Term_unique))

max_limit_global <- ceiling(max(abs(ggdata_combined_global$log_p_scaled), na.rm = TRUE))

# ==============================================================================
# Step12_Render_Global_Plot
# ==============================================================================
gg_all <- ggplot(ggdata_combined_global, aes(x = log_p_scaled, y = Term_unique, size = Significant)) +
  geom_point(aes(fill = Direction), shape = 21, colour = "white", stroke = 0.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7F8C8D", linewidth = 0.6) +
  scale_size_continuous(range = c(2, 9), name = "Gene count", breaks = pretty_breaks(n = 4)) +
  scale_fill_manual(values = color_modern, name = "Regulation") +
  scale_x_continuous(labels = abs, limits = c(-max_limit_global, max_limit_global), breaks = scales::pretty_breaks(n = 8)) +
  scale_y_discrete(expand = expansion(add = c(0.8, 1.1))) +
  labs(title = "GO Enrichment Analysis (Global View)", 
       subtitle = "Biological Processes by Nominal P-value: DOWN-regulated (left) vs UP-regulated (right)",
       x = expression(paste("-log"["10"], "(", italic("P"), "-value)")), y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(colour = "#e7e8e8", linewidth = 0.5, linetype = "dotted"),
    panel.grid.major.x = element_line(colour = "#dddddd", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(face = "bold", size = 16, color = "#2C3E50", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "#7F8C8D", margin = margin(b = 15)),
    axis.text.y = element_text(color = "#34495E", size = if(nrow(ggdata_combined_global) > 60) 6 else 9, margin = margin(r = 5)),
    axis.text.x = element_text(color = "#34495E", size = 10, face = "bold"),
    axis.title.x = element_text(color = "#2C3E50", size = 12, face = "bold", margin = margin(t = 12)),
    axis.ticks.x = element_line(colour = "#BDC3C7", linewidth = 0.5),
    axis.ticks.length.x = unit(0.15, "cm"),
    axis.line.x = element_line(colour = "#BDC3C7", linewidth = 0.5), 
    legend.position = "right",
    legend.title = element_text(face = "bold", color = "#2C3E50", size = 11),
    legend.text = element_text(color = "#34495E", size = 10),
    legend.key = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5, colour = "white")),
         size = guide_legend(override.aes = list(colour = "#7F8C8D")))

dynamic_height <- max(10, nrow(ggdata_combined_global) * 0.15)

print(gg_all)
ggsave("GO_enrichment_Global_View.pdf", gg_all, width = 14, height = dynamic_height, dpi = 600, limitsize = FALSE)
ggsave("GO_enrichment_Global_View.png", gg_all, width = 16, height = dynamic_height, dpi = 600, limitsize = FALSE)

# ==============================================================================
# Step13_Supplementary_Data_Export
# ==============================================================================
supplementary_go <- all_results_global %>%
  arrange(Direction, p.adj) %>% 
  dplyr::select(
    Direction,
    `GO_ID` = GO.ID,
    `Biological_Process` = Term,
    `Total_Annotated_Genes` = Annotated,
    `Significant_Genes_in_Dataset` = Significant,
    `Expected_Genes` = Expected,
    `Nominal_p_value` = Classic,
    `FDR_Adjusted_p_value` = p.adj
  )

wb_go <- createWorkbook()
addWorksheet(wb_go, "GO_Enrichment_Full")

header_style <- createStyle(
  textDecoration = "bold", 
  halign = "center", 
  fgFill = "#EFEFEF", 
  border = "TopBottomLeftRight"
)

writeData(wb_go, "GO_Enrichment_Full", supplementary_go, headerStyle = header_style)
setColWidths(wb_go, "GO_Enrichment_Full", cols = 1:ncol(supplementary_go), widths = "auto")
saveWorkbook(wb_go, "Supplementary_Data_GO_Enrichment.xlsx", overwrite = TRUE)