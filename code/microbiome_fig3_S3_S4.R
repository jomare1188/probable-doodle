# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(ggtree)
library(ALDEx2)
library(patchwork)
library(ggrepel) 
library(ggnewscale)
library(Polychrome)
library(RColorBrewer)
library(scales)
library(broom)

# ==============================================================================
# Step2_Hierarchical_Merge_Functions
# ==============================================================================
load_abundance <- function(path, suffix) {
  files <- list.files(path, pattern = "abundance.*\\.tsv$", full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) stop(paste("Error: File not found in", path))
  read_tsv(files[1], show_col_types = FALSE) %>%
    rename(tax = 1) %>%
    dplyr::select(-matches("total")) %>%
    rename_with(~paste0(., "_", suffix), -tax)
}

df1_silva    <- load_abundance("resultados_silva_138_2", "silva")
df2_16s      <- load_abundance("resultados_ncbi_16s_18s_28s_ITS", "n16s")
df3_plus8    <- load_abundance("resultados_ncbi_PlusPF-8", "p8")
df4_plus2026 <- load_abundance("resultados_ncbi_pluspf2026", "p2026")
df5_gg2      <- load_abundance("resultados_greengenes", "gg2")

df_merged <- df1_silva %>%
  full_join(df2_16s, by = "tax") %>% 
  full_join(df3_plus8, by = "tax") %>%
  full_join(df4_plus2026, by = "tax") %>% 
  full_join(df5_gg2, by = "tax")

samples <- c("control_1", "control_2", "control_3", "infected_1", "infected_2", "infected_3")

df_hierarchical <- df_merged %>%
  dplyr::select(tax) %>%
  bind_cols(map_dfc(samples, function(am) {
    tibble(!!am := coalesce(df_merged[[paste0(am, "_silva")]], df_merged[[paste0(am, "_n16s")]],
                            df_merged[[paste0(am, "_p8")]], df_merged[[paste0(am, "_p2026")]],
                            df_merged[[paste0(am, "_gg2")]]))
  })) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0)))

# ==============================================================================
# Step3_Biological_Curative_Filtering
# ==============================================================================
exclusion_terms <- c(
  "Homo", "sapiens", "Cutibacterium", "Leptailurus", "Eupleridae", "Timema", 
  "Babiana", "Ouratea", "Abiotrophia", "Akkermansia", "Alloprevotella", "Bacteroides", 
  "Bifidobacterium", "Blautia", "Dialister", "Eggerthella", "Faecalibacterium", 
  "Gardnerella", "Methanobrevibacter", "Neisseria", "Prevotella", "Ruminococcus", 
  "Rothia", "Veillonella", "Dehalococcoides", "Anaerococcus", "Malassezia", 
  "Candida", "Streptococcus"
)

df_separated <- df_hierarchical %>%
  filter(!str_detect(tax, "^(?i)(unclassified|unknown);(?i)(unknown|unclassified)")) %>%
  filter(!str_detect(tax, paste(exclusion_terms, collapse = "|"))) %>%
  dplyr::select(tax, all_of(samples)) %>%
  separate(tax, into = c("L1", "L2", "L3", "L4", "L5", "L6", "Genus"), sep = ";", fill = "right")

df_microbiome <- df_separated %>%
  filter(L1 %in% c("Bacteria", "Archaea") | 
           (L1 == "Eukaryota" & L3 %in% c("Ascomycota", "Basidiomycota", "Mucoromycota", "Zoopagomycota", "Chytridiomycota"))) %>%
  filter(L1 != "Viruses")

tax_fixed <- df_microbiome %>%
  dplyr::select(L1:Genus) %>%
  mutate(across(everything(), ~ ifelse(str_detect(., ".*(?i)(none|unknown|unclassified|no_rank|no rank).*") | is.na(.), "Unclassified", .)))

for (i in 2:ncol(tax_fixed)) {
  tax_fixed[[i]] <- ifelse(tax_fixed[[i]] == "Unclassified", paste0(tax_fixed[[i-1]], "_uncl."), tax_fixed[[i]])
}

df_clean_genus <- df_microbiome %>%
  dplyr::select(all_of(samples)) %>%
  bind_cols(tax_fixed) %>%
  mutate(Genus = str_replace_all(Genus, "(_[A-Za-z]|_[0-9]+|\\.[0-9]+)+$", "")) %>%
  mutate(Genus = str_replace_all(Genus, "(_uncl\\.)+", "_uncl.")) %>%
  mutate(total = rowSums(dplyr::select(., all_of(samples))))

tax_consensus <- df_clean_genus %>%
  group_by(Genus) %>%
  arrange(desc(total)) %>%
  slice(1) %>% 
  dplyr::select(L1, L2, L3, L4, L5, L6, Genus) %>%
  ungroup()

df_summed <- df_clean_genus %>%
  group_by(Genus) %>%
  summarise(across(all_of(samples), sum))

df_final_consolidated <- tax_consensus %>%
  left_join(df_summed, by = "Genus") %>%
  unite("tax", L1:Genus, sep = ";") %>%
  mutate(total = rowSums(dplyr::select(., all_of(samples)))) %>%
  arrange(desc(total))

write_tsv(df_final_consolidated, "abundance_table_genus_CONSOLIDATED.tsv")

# ==============================================================================
# Step4_Phyloseq_Object_Construction
# ==============================================================================
# Insert file!!!
raw_data <- read_tsv("insert corresponding file", show_col_types = FALSE) %>%
  dplyr::select(-total)

tax_matrix <- raw_data %>%
  dplyr::select(tax) %>%
  separate(tax, into = c("Kingdom", "L2", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", fill = "right") %>%
  as.matrix()

tax_matrix <- apply(tax_matrix, 2, str_trim)
tax_matrix[is.na(tax_matrix)] <- "Unclassified"

feature_ids <- make.unique(ifelse(str_detect(tax_matrix[,"Genus"], "_uncl.") | tax_matrix[,"Genus"] == "Unclassified", 
                                  tax_matrix[,"Family"], 
                                  tax_matrix[,"Genus"]))

rownames(tax_matrix) <- feature_ids
abundance_data <- as.matrix(raw_data[,-1])
rownames(abundance_data) <- feature_ids

metadata <- data.frame(
  treatment = factor(rep(c("control", "infected"), each = 3), levels = c("control", "infected")),
  row.names = colnames(abundance_data)
)

tax_dist <- vegan::taxa2dist(as.data.frame(tax_matrix), varstep = FALSE)
rooted_tree <- ape::as.phylo(hclust(tax_dist))

ps <- phyloseq(otu_table(abundance_data, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(tax_matrix),
               phy_tree(rooted_tree))

# ==============================================================================
# Step5_Differential_Expression_ALDEx2
# ==============================================================================
conds <- as.character(sample_data(ps)$treatment)
aldex_res <- aldex(otu_table(ps), conds, mc.samples = 128, test = "t", effect = TRUE, denom = "all")

# ==============================================================================
# Step6_Ordination_Differential_Taxa_PCoA
# ==============================================================================
significant_taxa_with_tax <- aldex_res %>%
  rownames_to_column(var = "Taxa_ID") %>%
  filter(we.ep < 0.05 & abs(effect) > 1.0)

thesis_colors <- c("control" = "#407126", "infected" = "#A23232")

if (exists("significant_taxa_with_tax") && nrow(significant_taxa_with_tax) > 0) {
  significant_genus_names <- significant_taxa_with_tax$Taxa_ID
  significant_genus_names <- significant_genus_names[!is.na(significant_genus_names)]
  
  if (length(significant_genus_names) > 0) {
    ps_significant <- prune_taxa(significant_genus_names, ps)
    
    if (ntaxa(ps_significant) > 0 && nsamples(ps_significant) > 0) {
      ps_clr <- ps_significant 
      otu_mat <- as(otu_table(ps_clr), "matrix")
      otu_mat_pseudo <- otu_mat + 1 
      gm_mean_per_sample <- apply(otu_mat_pseudo, 2, function(x) exp(mean(log(x))))
      otu_clr_mat <- sweep(log(otu_mat_pseudo), 2, log(gm_mean_per_sample), FUN = "-")
      otu_table(ps_clr) <- otu_table(as.matrix(otu_clr_mat), taxa_are_rows = TRUE)
      
      dist_aitchison_significant <- phyloseq::distance(ps_clr, method = "euclidean")
      pcoa_aitchison_significant <- ordinate(ps_clr, method = "PCoA", distance = dist_aitchison_significant)
      
      plot_df_pcoa <- as.data.frame(pcoa_aitchison_significant$vectors[, 1:2])
      colnames(plot_df_pcoa) <- c("Axis.1", "Axis.2")
      plot_df_pcoa$treatment <- sample_data(ps_clr)$treatment 
      plot_df_pcoa$Sample <- rownames(plot_df_pcoa)
      
      get_ellipse <- function(x, y, conf = 0.95) {
        cov_mat <- cov(cbind(x, y))
        center <- c(mean(x), mean(y))
        eig <- eigen(cov_mat)
        angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
        radius <- sqrt(qchisq(conf, df = 2))
        a <- radius * sqrt(eig$values[1])
        b <- radius * sqrt(eig$values[2])
        theta <- seq(0, 2 * pi, length.out = 100)
        x_ell <- center[1] + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
        y_ell <- center[2] + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)
        return(data.frame(Axis.1 = x_ell, Axis.2 = y_ell))
      }
      
      ellipse_list <- by(plot_df_pcoa, plot_df_pcoa$treatment, function(df) {
        ell <- get_ellipse(df$Axis.1, df$Axis.2, conf = 0.95)
        ell$treatment <- df$treatment[1] 
        return(ell)
      })
      ellipse_data <- do.call(rbind, ellipse_list)
      
      plot_pcoa_significant <- ggplot(plot_df_pcoa, aes(x = Axis.1, y = Axis.2, color = treatment)) +
        geom_polygon(data = ellipse_data, aes(x = Axis.1, y = Axis.2, fill = treatment), alpha = 0.2, color = NA, show.legend = FALSE) +
        geom_point(size = 3, alpha = 0.8) + 
        scale_color_manual(values = thesis_colors, labels = c("Control", "Infected")) +
        scale_fill_manual(values = thesis_colors) +
        labs(
          title = "PCoA (Aitchison) - Differential Genera",
          subtitle = paste0("Composition of ", ntaxa(ps_significant), " Significant Genera"),
          x = paste0("PC1 (", round(pcoa_aitchison_significant$values$Relative_eig[1] * 100, 1), "%)"),
          y = paste0("PC2 (", round(pcoa_aitchison_significant$values$Relative_eig[2] * 100, 1), "%)")
        ) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "bottom"
        )
      
      if(!dir.exists("results_merged_R")) dir.create("results_merged_R")
      ggsave("results_merged_R/pcoa_aitchison_significant_genus_clr.png", plot_pcoa_significant, width = 10, height = 8, dpi = 600, bg = "white")
    }
  }
}

# ==============================================================================
# Step7_Alpha_Diversity_Metrics
# ==============================================================================
alpha_diversity <- estimate_richness(ps, measures = "Shannon")
alpha_diversity$Sample <- rownames(alpha_diversity)

sample_info <- as(sample_data(ps), "data.frame")
sample_info$Sample <- rownames(sample_info)
alpha_div_with_treatment <- merge(alpha_diversity, sample_info, by = "Sample")

alpha_data <- alpha_div_with_treatment %>%
  dplyr::select(Sample, treatment, Shannon) %>%
  dplyr::mutate(treatment = factor(treatment, levels = c("control", "infected"), labels = c("Control", "Infected")))

shannon_stats <- alpha_data %>%
  rstatix::t_test(Shannon ~ treatment, var.equal = TRUE) %>%
  dplyr::mutate(
    significance = case_when(
      p < 0.001 ~ "***",
      p < 0.01 ~ "**", 
      p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write_tsv(shannon_stats, "results_merged_R/SUPP_Alpha_Diversity_T_Test.tsv")

if (shannon_stats$p < 0.05) {
  alpha_data$letters <- ifelse(alpha_data$treatment == "Control", "a", "b")
} else {
  alpha_data$letters <- "a"
}

alpha_summary <- alpha_data %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(max_value = max(Shannon, na.rm = TRUE), letter = unique(letters)[1], .groups = "drop") %>%
  dplyr::mutate(y_position = max_value + max_value * 0.05)

alpha_data_with_positions <- alpha_data %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(
    group_index = row_number(),
    n_samples = n(),
    x_offset = (group_index - (n_samples + 1)/2) * 0.1,  
    x_position = as.numeric(treatment) + x_offset
  ) %>%
  dplyr::ungroup()

alpha_plot <- ggplot(alpha_data, aes(x = treatment, y = Shannon, fill = treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(data = alpha_data_with_positions, aes(x = x_position, y = Shannon), size = 2, alpha = 0.7, color = "black", inherit.aes = FALSE) +
  geom_text(data = alpha_summary, aes(x = treatment, y = y_position, label = letter), inherit.aes = FALSE, size = 10/2.83465, family = "sans", fontface = "bold") +
  scale_fill_manual(values = thesis_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 9, family = "sans", color = "black", face="bold"),
    axis.title.y = element_text(size = 9, family = "sans", color = "black", face="bold"),
    axis.text.x = element_text(size = 8, family = "sans", color = "black"),
    axis.text.y = element_text(size = 8, family = "sans", color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm")
  ) +
  labs(
    x = "Treatment",
    y = "Shannon Diversity Index",
    title = "Alpha Diversity Comparison",
    subtitle = paste0("Student's t-test: p = ", round(shannon_stats$p, 4), " (", shannon_stats$significance, ")")
  )

ggsave("results_merged_R/alpha_diversity_shannon.png", plot = alpha_plot, width = 7, height = 6, dpi = 600)

# ==============================================================================
# Step8_Beta_Diversity_Distance_Matrices
# ==============================================================================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
otu_matrix <- as(otu_table(ps_rel), "matrix")
if (taxa_are_rows(ps_rel)) {
  otu_matrix <- t(otu_matrix)
}

bray_dist <- vegdist(otu_matrix, method = "bray")
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_variance <- pcoa_result$eig[1:2] / sum(pcoa_result$eig) * 100

pcoa_data <- data.frame(
  Sample = rownames(pcoa_result$points),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2]
)

pcoa_data <- merge(pcoa_data, sample_info[, c("Sample", "treatment")], by = "Sample")
pcoa_data$treatment <- factor(pcoa_data$treatment, levels = c("control", "infected"), labels = c("Control", "Infected"))
sample_info_ordered <- sample_info[match(pcoa_data$Sample, sample_info$Sample), ]

num_permutations <- min(999, factorial(length(pcoa_data$Sample)) / (factorial(table(pcoa_data$treatment)[1]) * factorial(table(pcoa_data$treatment)[2])))
if (is.infinite(num_permutations) || num_permutations < 10) { 
  num_permutations <- 19 
}

permanova_result <- adonis2(bray_dist ~ treatment, data = sample_info_ordered, permutations = num_permutations, method = "bray")
write_tsv(broom::tidy(permanova_result), "results_merged_R/SUPP_Beta_Diversity_PERMANOVA.tsv")

r2_value <- round(permanova_result$R2[1], 3)
p_value <- permanova_result$`Pr(>F)`[1]
p_text <- if (p_value < 0.001) "p < 0.001" else paste("p =", round(p_value, 3))

create_sd_ellipse <- function(data_subset, treatment_name, sd_factor = 1.5) {
  if (nrow(data_subset) < 2) {
    return(data.frame(PC1 = numeric(0), PC2 = numeric(0), treatment = character(0)))
  }
  center_x <- mean(data_subset$PC1, na.rm = TRUE)
  center_y <- mean(data_subset$PC2, na.rm = TRUE)
  sd_x <- sd(data_subset$PC1, na.rm = TRUE)
  sd_y <- sd(data_subset$PC2, na.rm = TRUE)
  
  if (is.na(sd_x) || sd_x == 0) sd_x <- 0.02
  if (is.na(sd_y) || sd_y == 0) sd_y <- 0.02
  
  angles <- seq(0, 2*pi, length.out = 100)
  ellipse_data <- data.frame(
    PC1 = center_x + sd_factor * sd_x * cos(angles),
    PC2 = center_y + sd_factor * sd_y * sin(angles),
    treatment = rep(treatment_name, 100)
  )
  return(ellipse_data)
}

control_data <- pcoa_data[pcoa_data$treatment == "Control", ]
infected_data <- pcoa_data[pcoa_data$treatment == "Infected", ]

ellipse_sd_control <- create_sd_ellipse(control_data, "Control", sd_factor = 1.5)
ellipse_sd_infected <- create_sd_ellipse(infected_data, "Infected", sd_factor = 1.5)
ellipse_sd_data <- rbind(ellipse_sd_control, ellipse_sd_infected)

beta_plot <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_polygon(data = ellipse_sd_data, aes(x = PC1, y = PC2, fill = treatment), alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = thesis_colors) +
  scale_fill_manual(values = thesis_colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 9, family = "sans", color = "black", face="bold"),
    axis.title.y = element_text(size = 9, family = "sans", color = "black", face="bold"),
    axis.text.x = element_text(size = 8, family = "sans", color = "black"),
    axis.text.y = element_text(size = 8, family = "sans", color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm")
  ) +
  labs(
    title = "Beta Diversity (Bray-Curtis PCoA)",
    subtitle = paste0("PERMANOVA: R² = ", r2_value, ", ", p_text),
    x = paste0("PC1 (", round(pcoa_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pcoa_variance[2], 1), "%)")
  )

ggsave("results_merged_R/beta_diversity_bray.png", plot = beta_plot, width = 7, height = 6, dpi = 600)

# ==============================================================================
# Step9_Phylogenetic_Tree_Ring_Setup
# ==============================================================================
tree <- phy_tree(ps)
tree_transformed <- tree
tree_transformed$edge.length <- sqrt(tree$edge.length)

tax_table_df <- as.data.frame(tax_table(ps)) %>% tibble::rownames_to_column(var = "OTU_ID")
tree_tips <- tree$tip.label
tax_table_filtered <- tax_table_df %>% dplyr::filter(OTU_ID %in% tree_tips)

class_data <- tax_table_filtered %>%
  dplyr::mutate(Class_Display = ifelse(!is.na(Class) & Class != "Unclassified", Class, "Unclassified")) %>%
  dplyr::select(OTU_ID, Class_Display)

otu_table_df <- as.data.frame(otu_table(ps)) %>%
  tibble::rownames_to_column(var = "OTU_ID") %>%
  dplyr::filter(OTU_ID %in= tree_tips)

sample_data_df <- as_tibble(as.data.frame(sample_data(ps)), rownames = "Sample_ID")

presence_data <- otu_table_df %>%
  tidyr::pivot_longer(cols = -OTU_ID, names_to = "Sample_ID", values_to = "Abundance") %>%
  dplyr::left_join(sample_data_df %>% dplyr::select(Sample_ID, treatment), by = "Sample_ID") %>%
  dplyr::filter(Abundance > 0) %>%
  dplyr::group_by(OTU_ID, treatment) %>%
  dplyr::summarise(Count = n(), .groups = "drop") %>%
  dplyr::group_by(OTU_ID) %>%
  dplyr::summarise(
    Present_in_Control = any(treatment == "control"),
    Present_in_Infected = any(treatment == "infected"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Presence_Status = case_when(
      Present_in_Control & Present_in_Infected ~ "Both",
      Present_in_Control ~ "Control",
      Present_in_Infected ~ "Infected",
      TRUE ~ "None"
    )
  ) %>%
  dplyr::select(OTU_ID, Presence_Status)

significance_data <- aldex_res %>%
  tibble::rownames_to_column(var = "OTU_ID") %>%
  dplyr::filter(OTU_ID %in% tree_tips) %>%
  dplyr::mutate(Significant = case_when(we.ep < 0.05 & abs(effect) > 1.0 ~ "Significant", TRUE ~ "Non-significant")) %>%
  dplyr::select(OTU_ID, Significant)

n_classes <- length(unique(class_data$Class_Display))
if (n_classes <= 12) {
  library(RColorBrewer)
  class_colors <- brewer.pal(min(n_classes, 12), "Set3")
  if (n_classes < 12) class_colors <- class_colors[1:n_classes]
} else {
  set.seed(123) 
  class_colors <- createPalette(N = n_classes, seedcolors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"), range = c(30, 300), target = "normal", M = 5000)
}
names(class_colors) <- unique(class_data$Class_Display)

presence_colors <- c("Control" = "#407126", "Infected" = "#A23232", "Both" = "#674a42", "None" = "#BDC3C7")
significance_colors <- c("Significant" = "#9c1a15", "Non-significant" = "#D7D6D2", "Not tested" = "#ECF0F1")

class_data_for_ring <- class_data %>% dplyr::mutate(value = 1) %>% dplyr::rename(taxa = OTU_ID)
presence_data_for_ring <- presence_data %>% dplyr::mutate(value = 1) %>% dplyr::rename(taxa = OTU_ID)
significance_data_for_ring <- significance_data %>% dplyr::mutate(value = 1) %>% dplyr::rename(taxa = OTU_ID)

p_tree <- ggtree(tree_transformed, layout = "circular", linewidth = 0.2) +
  geom_fruit(data = class_data_for_ring, geom = geom_tile, mapping = aes(y = taxa, x = value, fill = Class_Display), width = 2.0, offset = 0.10, pwidth = 0.02, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = class_colors, name = "Class") +
  new_scale_fill() +
  geom_fruit(data = presence_data_for_ring, geom = geom_tile, mapping = aes(y = taxa, x = value, fill = Presence_Status), width = 2.0, offset = 0.13, pwidth = 0.02, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = presence_colors, name = "Presence") +
  new_scale_fill() +
  geom_fruit(data = significance_data_for_ring, geom = geom_tile, mapping = aes(y = taxa, x = value, fill = Significant), width = 2.0, offset = 0.16, pwidth = 0.02, color = "white", linewidth = 0.1) +
  scale_fill_manual(values = significance_colors, name = "ALDEx2 Significance") +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 8), legend.title = element_text(size = 10, face = "bold"), legend.key.size = unit(0.4, "cm"))

ggsave("results_merged_R/phylogenetic_tree_rings.png", plot = p_tree, width = 12, height = 10, dpi = 900, bg = "white")

# ==============================================================================
# Step10_Phylogenetic_Tree_Heatmap_Full
# ==============================================================================
# Insert file!!!
raw_data_heat <- read_tsv("insert corresponding file", show_col_types = FALSE) %>%
  dplyr::select(-total)

tax_matrix_heat <- raw_data_heat %>%
  dplyr::select(tax) %>%
  separate(tax, into = c("Kingdom", "L2", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", fill = "right") %>%
  as.matrix()

tax_matrix_heat <- apply(tax_matrix_heat, 2, str_trim)
tax_matrix_heat[is.na(tax_matrix_heat)] <- "Unclassified"

feature_ids_heat <- make.unique(ifelse(tax_matrix_heat[,"Genus"] == "Unclassified" | str_detect(tax_matrix_heat[,"Genus"], "_uncl"), 
                                       tax_matrix_heat[,"Family"], 
                                       tax_matrix_heat[,"Genus"]))

abundance_data_heat <- as.matrix(raw_data_heat[,-1])
rownames(abundance_data_heat) <- feature_ids_heat
rownames(tax_matrix_heat) <- feature_ids_heat

metadata_heat <- data.frame(treatment = factor(rep(c("control", "infected"), each = 3)), row.names = colnames(abundance_data_heat))
ps_heat <- phyloseq(otu_table(abundance_data_heat, taxa_are_rows = TRUE), tax_table(tax_matrix_heat), sample_data(metadata_heat))

tax_dist_heat <- dist.gene(as.data.frame(tax_table(ps_heat)))
ps_heat <- merge_phyloseq(ps_heat, as.phylo(hclust(tax_dist_heat)))

conds_heat <- as.character(sample_data(ps_heat)$treatment)
aldex_res_heat <- aldex(otu_table(ps_heat), conds_heat, mc.samples = 128, test = "t", effect = TRUE, denom = "all")

sig_results <- aldex_res_heat %>%
  rownames_to_column(var = "Genus_Match") %>%
  mutate(Genus_Match = str_trim(Genus_Match)) %>%
  filter(we.ep < 0.05 & abs(effect) > 2.5)

ps_rel_all <- transform_sample_counts(ps_heat, function(x) x / sum(x))
ps_sig <- prune_taxa(sig_results$Genus_Match, ps_rel_all)

if (ntaxa(ps_sig) > 0) {
  eff_vec <- sig_results$effect
  names(eff_vec) <- sig_results$Genus_Match
  
  tree_ordered <- reorder(phy_tree(ps_sig), "postorder")
  phy_tree(ps_sig) <- tree_ordered
  
  tax_df_sig <- as.data.frame(tax_table(ps_sig)) %>%
    rownames_to_column(var = "label") %>%
    mutate(label_display = ifelse(Genus != "Unclassified", Genus, Family))
  
  abund_df <- as.data.frame(otu_table(ps_sig)) %>%
    rownames_to_column(var = "label") %>%
    pivot_longer(-label, names_to = "Sample", values_to = "Value") %>%
    mutate(Treatment = ifelse(str_detect(Sample, "control"), "Control", "Infected")) %>%
    group_by(label, Treatment) %>%
    summarise(MeanAbund = log10(mean(Value) + 1e-6), .groups = "drop")
  
  effect_df <- sig_results %>%
    rename(label = Genus_Match) %>%
    filter(label %in% taxa_names(ps_sig))
  
  tree_plot <- ggtree(phy_tree(ps_sig), layout = "rectangular", linewidth = 0.5) %<+% tax_df_sig +
    geom_tiplab(aes(label = label_display), size = 6, offset = 0.02, fontface = "italic") +
    hexpand(0.3) + theme_tree2()
  
  tip_order <- tree_plot$data %>% filter(isTip) %>% arrange(y) %>% pull(label)
  mid_abund <- median(abund_df$MeanAbund, na.rm = TRUE) 
  
  heatmap_abund <- ggplot(abund_df, aes(x = Treatment, y = label, fill = MeanAbund)) + 
    geom_tile(color = NA) + 
    scale_y_discrete(limits = tip_order) + 
    scale_fill_gradient2(
      low = "#0000FF", mid = "#F0FFFF", high = "#DC143C", midpoint = -3.5, 
      limits = c(-6, -1.0), oob = scales::squish, name = "Log10\nAbund."
    ) +
    theme_void() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  
  heatmap_eff <- ggplot(effect_df, aes(x = "Effect", y = label, fill = effect)) + 
    geom_tile(color = NA) + 
    scale_y_discrete(limits = tip_order) + 
    scale_fill_gradient2(low = "#116968", mid = "white", high = "#340F5A", limits = c(-6, 6), oob = scales::squish, name = "Effect") +
    theme_void() + theme(axis.text.x = element_text(size = 6))
  
  final_fig <- tree_plot + heatmap_abund + heatmap_eff + plot_layout(widths = c(2.2, 0.9, 0.4))
  ggsave("results_merged_R/tree_heatmap_full.png", final_fig, width = 15, height = 25, dpi = 900, bg = "white")
}

# ==============================================================================
# Step11_Phylogenetic_Tree_Heatmap_Top_Taxa
# ==============================================================================
vip_taxa <- c("Fusarium", "Streptomyces", "Pantoea", "Microbacterium", "Arthrobacter", "Curtobacterium", "Leifsonia", "Paenibacillus")

aldex_formatted_511 <- aldex_res %>%
  rownames_to_column(var = "Genus_Match") %>%
  mutate(Genus_Match = str_trim(Genus_Match))

sig_vip_511 <- aldex_formatted_511 %>% filter(Genus_Match %in% vip_taxa)
sig_top50_511 <- aldex_formatted_511 %>% filter(!Genus_Match %in% vip_taxa) %>% filter(we.ep < 0.05) %>% arrange(desc(abs(effect))) %>% slice_head(n = 30)                         
sig_results_511 <- bind_rows(sig_vip_511, sig_top50_511) %>% distinct(Genus_Match, .keep_all = TRUE)    

ps_rel_all_511 <- transform_sample_counts(ps, function(x) x / sum(x))
ps_sig_511 <- prune_taxa(sig_results_511$Genus_Match, ps_rel_all_511)

if (ntaxa(ps_sig_511) > 0) {
  eff_vec_511 <- sig_results_511$effect
  names(eff_vec_511) <- sig_results_511$Genus_Match
  
  tree_ordered_511 <- reorder(phy_tree(ps_sig_511), "postorder")
  phy_tree(ps_sig_511) <- tree_ordered_511
  
  tax_df_sig_511 <- as.data.frame(tax_table(ps_sig_511)) %>%
    rownames_to_column(var = "label") %>%
    mutate(label_display = ifelse(Genus != "Unclassified", Genus, Family))
  
  abund_df_511 <- as.data.frame(otu_table(ps_sig_511)) %>%
    rownames_to_column(var = "label") %>%
    pivot_longer(-label, names_to = "Sample", values_to = "Value") %>%
    mutate(Treatment = ifelse(str_detect(Sample, "control"), "Control", "Infected")) %>%
    group_by(label, Treatment) %>%
    summarise(MeanAbund = log10(mean(Value) + 1e-6), .groups = "drop")
  
  effect_df_511 <- sig_results_511 %>% rename(label = Genus_Match) %>% filter(label %in% taxa_names(ps_sig_511))
  
  tree_plot_511 <- ggtree(phy_tree(ps_sig_511), layout = "rectangular", linewidth = 0.5) %<+% tax_df_sig_511 +
    geom_tiplab(aes(label = label_display), size = 6, offset = 0.02, fontface = "italic") +
    hexpand(0.3) + theme_tree2()
  
  tip_order_511 <- tree_plot_511$data %>% filter(isTip) %>% arrange(y) %>% pull(label)
  mid_abund_511 <- median(abund_df_511$MeanAbund, na.rm = TRUE) 
  
  heatmap_abund_511 <- ggplot(abund_df_511, aes(x = Treatment, y = label, fill = MeanAbund)) + 
    geom_tile(color = NA) + 
    scale_y_discrete(limits = tip_order_511) + 
    scale_fill_gradient2(low = "#0000FF", mid = "#F0FFFF", high = "#DC143C", midpoint = -3.5, limits = c(-6, -1.0), oob = scales::squish, name = "Log10\nAbund.") +
    theme_void() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  heatmap_eff_511 <- ggplot(effect_df_511, aes(x = "Effect", y = label, fill = effect)) + 
    geom_tile(color = NA) + 
    scale_y_discrete(limits = tip_order_511) + 
    scale_fill_gradient2(low = "#116968", mid = "white", high = "#340F5A", limits = c(-6, 6), oob = scales::squish, name = "Effect") +
    theme_void() + theme(axis.text.x = element_text(size = 8))
  
  final_fig_511 <- tree_plot_511 + heatmap_abund_511 + heatmap_eff_511 + plot_layout(widths = c(2.2, 0.9, 0.4))
  ggsave("results_merged_R/tree_heatmap_Nature_Top30.png", final_fig_511, width = 12, height = 18, dpi = 900, bg = "white")
}

# ==============================================================================
# Step12_Master_Data_Export
# ==============================================================================
if(!dir.exists("results_merged_R")) dir.create("results_merged_R")

tax_full <- as.data.frame(tax_table(ps)) %>% rownames_to_column(var = "Genus_ID")

counts_raw <- as.data.frame(otu_table(ps)) %>%
  rownames_to_column(var = "Genus_ID") %>%
  rename_with(~paste0(.x, "_RawCounts"), -Genus_ID)

ps_rel_export <- transform_sample_counts(ps, function(x) (x / sum(x)) * 100)
counts_rel <- as.data.frame(otu_table(ps_rel_export)) %>%
  rownames_to_column(var = "Genus_ID") %>%
  rename_with(~paste0(.x, "_RelAbund_pct"), -Genus_ID)

aldex_full <- aldex_res %>%
  rownames_to_column(var = "Genus_ID") %>%
  mutate(
    Significancia_ALDEx2 = case_when(
      we.ep < 0.05 & abs(effect) > 1.0 ~ "Significant (Strong)",
      we.ep < 0.05 & abs(effect) <= 1.0 ~ "Significant (Weak)",
      TRUE ~ "Non Significant"
    )
  )

master_taxa_table <- tax_full %>%
  left_join(counts_raw, by = "Genus_ID") %>%
  left_join(counts_rel, by = "Genus_ID") %>%
  left_join(aldex_full, by = "Genus_ID") %>%
  arrange(we.ep)

write_tsv(master_taxa_table, "results_merged_R/MASTER_TABLE_Taxa_Full_Stats.tsv")

if(exists("alpha_data") & exists("pcoa_data")) {
  master_sample_table <- alpha_data %>%
    dplyr::select(Sample, treatment, Shannon) %>%
    left_join(pcoa_data %>% dplyr::select(Sample, PC1, PC2), by = "Sample") %>%
    rename(Tratamento = treatment, Diversidade_Shannon = Shannon, PCoA_Bray_PC1 = PC1, PCoA_Bray_PC2 = PC2)
  
  write_tsv(master_sample_table, "results_merged_R/MASTER_TABLE_Samples_Diversity.tsv")
}

if(exists("abund_df") & exists("effect_df") & exists("tax_df_sig")) {
  abund_wide <- abund_df %>%
    pivot_wider(names_from = Treatment, values_from = MeanAbund) %>%
    rename(Log10_Control = Control, Log10_Infected = Infected)
  
  df_fig_510 <- abund_wide %>%
    left_join(effect_df, by = "label") %>%
    left_join(tax_df_sig, by = "label") %>%
    dplyr::select(
      ID_Taxon = label, Nome_Exibicao = label_display, Log10_Control, Log10_Infected,
      ALDEx2_Effect = effect, p_valor_Welch = we.ep, p_valor_ajustado = we.eBH,
      Kingdom, Phylum, Class, Order, Family, Genus
    ) %>%
    arrange(desc(ALDEx2_Effect))
  
  write_tsv(df_fig_510, "results_merged_R/SUMMARY_Figura_5_10_Tree.tsv")
}