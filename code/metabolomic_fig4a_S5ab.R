# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(RColorBrewer)
library(scales)
library(patchwork)
library(janitor)
library(ggrepel)
library(openxlsx)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
data_raw <- read_excel("insert corresponding file", sheet = 1)
data_raw <- data_raw %>% clean_names()

# ==============================================================================
# Step3_Data_Transformation
# ==============================================================================
expected_cols <- c("treatment", "rep")
missing_cols <- setdiff(expected_cols, names(data_raw))

if(length(missing_cols) > 0) {
  stop("Expected columns not found: ", paste(missing_cols, collapse = ", "))
}

metabolite_cols <- names(data_raw)[!names(data_raw) %in% c("treatment", "rep")]

data_clean <- data_raw %>%
  mutate(
    treatment = factor(treatment, levels = c("non-contaminated", "F. verticillioides")),
    rep = as.numeric(rep),
    sample_id = paste0(treatment, "_", rep)
  ) %>%
  mutate(across(all_of(metabolite_cols), ~ {
    if(is.character(.x)) {
      .x <- gsub(",", ".", .x)
      .x <- gsub("E\\+", "e+", .x)
      .x <- gsub("E-", "e-", .x)
      as.numeric(.x)
    } else {
      as.numeric(.x)
    }
  })) %>%
  pivot_longer(
    cols = all_of(metabolite_cols),
    names_to = "metabolite",
    values_to = "concentration"
  ) %>%
  filter(!is.na(concentration)) %>%
  mutate(
    log_concentration = log10(concentration + 1),
    is_zero = concentration == 0,
    log2_concentration = log2(concentration + 1)
  )

# ==============================================================================
# Step4_Statistical_Testing_Function
# ==============================================================================
perform_statistical_test <- function(data) {
  n_groups <- n_distinct(data$treatment)
  min_samples <- min(table(data$treatment))
  
  if(n_groups < 2 || min_samples < 2) {
    return(tibble(
      p_value = NA_real_,
      mean_control = NA_real_,
      mean_treatment = NA_real_,
      fold_change = NA_real_,
      method = "Insufficient data"
    ))
  }
  
  if(var(data$log_concentration, na.rm = TRUE) == 0) {
    means <- data %>%
      group_by(treatment) %>%
      summarise(mean_log = mean(log_concentration, na.rm = TRUE), .groups = "drop")
    
    mean_control_val <- means$mean_log[means$treatment == "non-contaminated"]
    mean_treatment_val <- means$mean_log[means$treatment == "F. verticillioides"]
    
    return(tibble(
      p_value = 1,
      mean_control = mean_control_val,
      mean_treatment = mean_treatment_val,
      fold_change = 0,
      method = "No variation"
    ))
  }
  
  tryCatch({
    test_result <- t.test(log_concentration ~ treatment, data = data)
    
    means <- data %>%
      group_by(treatment) %>%
      summarise(mean_log = mean(log_concentration, na.rm = TRUE), .groups = "drop")
    
    mean_control_val <- means$mean_log[means$treatment == "non-contaminated"]
    mean_treatment_val <- means$mean_log[means$treatment == "F. verticillioides"]
    
    tibble(
      p_value = test_result$p.value,
      mean_control = mean_control_val,
      mean_treatment = mean_treatment_val,
      fold_change = mean_treatment_val - mean_control_val,
      method = "t-test"
    )
  }, error = function(e) {
    tibble(
      p_value = NA_real_,
      mean_control = NA_real_,
      mean_treatment = NA_real_,
      fold_change = NA_real_,
      method = paste("Error:", e$message)
    )
  })
}

# ==============================================================================
# Step5_Execute_Differential_Analysis
# ==============================================================================
stats_results <- data_clean %>%
  group_by(metabolite) %>%
  do(perform_statistical_test(.)) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    significance = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**", 
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    log2_fc = fold_change / log10(2),
    abs_log2_fc = abs(log2_fc)
  )

write_csv(stats_results, "statistical_results.csv")

# ==============================================================================
# Step6_Conditional_Violin_Boxplot_Function
# ==============================================================================
create_individual_boxplot <- function(metabolite_name, save_plot = TRUE) {
  plot_data <- data_clean %>%
    filter(metabolite == metabolite_name)
  
  stat_info <- stats_results %>%
    filter(metabolite == metabolite_name)
  
  summary_stats <- plot_data %>%
    group_by(treatment) %>%
    summarise(
      n = n(),
      mean_val = mean(log_concentration, na.rm = TRUE),
      median_val = median(log_concentration, na.rm = TRUE),
      sd_val = sd(log_concentration, na.rm = TRUE),
      se_val = sd_val / sqrt(n),
      n_zeros = sum(log_concentration == 0),
      all_zero = all(log_concentration == 0),
      .groups = "drop"
    )
  
  clean_title <- str_replace_all(metabolite_name, "_", " ") %>% 
    str_to_title() %>%
    str_replace("C(\\d+)", "C\\1:") %>%
    str_replace("X7", "7,") %>%
    str_replace("D2", " D2") %>%
    str_replace("E$", " E")
  
  control_all_zero <- summary_stats %>%
    filter(treatment == "non-contaminated") %>%
    pull(all_zero)
  
  treatment_all_zero <- summary_stats %>%
    filter(treatment == "F. verticillioides") %>%
    pull(all_zero)
  
  all_values_zero <- all(plot_data$log_concentration == 0)
  
  control_data <- plot_data %>% filter(treatment == "non-contaminated")
  treatment_data <- plot_data %>% filter(treatment == "F. verticillioides")
  
  min_y <- min(plot_data$log_concentration, na.rm = TRUE)
  max_y <- max(plot_data$log_concentration, na.rm = TRUE)
  range_y <- max_y - min_y
  
  if(all_values_zero) {
    y_min_plot <- -0.05
    y_max_plot <- 0.15
    y_annotation <- 0.1
    
  } else if(control_all_zero && !treatment_all_zero) {
    if(max_y <= 1) {
      y_min_plot <- -0.05
      y_max_plot <- max_y + (max_y * 0.15)
      y_annotation <- max_y + (max_y * 0.08)
    } else if(max_y <= 5) {
      y_min_plot <- -0.1
      y_max_plot <- max_y + (max_y * 0.12)
      y_annotation <- max_y + (max_y * 0.06)
    } else {
      y_min_plot <- -0.15
      y_max_plot <- max_y + (max_y * 0.1)
      y_annotation <- max_y + (max_y * 0.05)
    }
    
  } else if(!control_all_zero && treatment_all_zero) {
    if(max_y <= 1) {
      y_min_plot <- -0.05
      y_max_plot <- max_y + (max_y * 0.15)
      y_annotation <- max_y + (max_y * 0.08)
    } else {
      y_min_plot <- -0.1
      y_max_plot <- max_y + (max_y * 0.12)
      y_annotation <- max_y + (max_y * 0.06)
    }
    
  } else {
    if(range_y < 0.5) {
      y_buffer <- 0.15
    } else {
      y_buffer <- range_y * 0.1
    }
    
    y_min_plot <- max(0, min_y - y_buffer)
    y_max_plot <- max_y + y_buffer
    y_annotation <- max_y + (y_buffer * 0.5)
  }
  
  p_label <- if(!is.na(stat_info$p_adj)) {
    if(stat_info$p_adj < 0.001) {
      "p < 0.001"
    } else {
      paste0("p = ", format(round(stat_info$p_adj, 3), nsmall = 3))
    }
  } else {
    "p = NA"
  }
  
  p <- ggplot(plot_data, aes(x = treatment, y = log_concentration, fill = treatment))
  
  if(!control_all_zero) {
    p <- p + geom_violin(data = control_data, alpha = 0.6, trim = FALSE, scale = "width")
  }
  
  if(!treatment_all_zero) {
    p <- p + geom_violin(data = treatment_data, alpha = 0.6, trim = FALSE, scale = "width")
  }
  
  p <- p +
    geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, color = "black", linewidth = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.8, size = 2.0, color = "black", stroke = 0.3) +
    scale_fill_manual(
      values = c("non-contaminated" = "#407126", "F. verticillioides" = "#957e2a"),
      labels = c("Control", "F. verticillioides")
    ) +
    coord_cartesian(ylim = c(y_min_plot, y_max_plot)) +
    labs(
      title = clean_title,
      subtitle = paste0("n = ", summary_stats$n[1], " vs ", summary_stats$n[2], " | ", p_label),
      x = "Treatment",
      y = expression(Log[10]~"(Concentration + 1)"),
      fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.75),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    {if(!is.na(stat_info$p_adj) && stat_info$p_adj < 0.05) {
      annotate("segment", x = 1, xend = 2, y = y_annotation, yend = y_annotation, color = "black", linewidth = 0.5)
    }} +
    {if(!is.na(stat_info$p_adj) && stat_info$p_adj < 0.05) {
      annotate("text", x = 1.5, y = y_annotation + (y_max_plot - y_min_plot) * 0.02, 
               label = p_label, size = 4, hjust = 0.5, fontface = "bold", color = "red")
    }}
  
  if((control_all_zero || treatment_all_zero) && y_min_plot < 0 && !all_values_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", color = "gray60", alpha = 0.6, linewidth = 0.3)
  }
  
  if(save_plot) {
    dir.create("individual_boxplots", showWarnings = FALSE)
    filename_safe <- str_replace_all(metabolite_name, "[^A-Za-z0-9_]", "_")
    ggsave(paste0("individual_boxplots/", filename_safe, "_boxplot.png"), plot = p, width = 4, height = 8, dpi = 600, bg = "white")
    ggsave(paste0("individual_boxplots/", filename_safe, "_boxplot.pdf"), plot = p, width = 6, height = 6)
  }
  
  return(p)
}

# ==============================================================================
# Step7_Compact_Violin_Boxplot_Function
# ==============================================================================
create_compact_boxplot <- function(metabolite_name) {
  plot_data <- data_clean %>%
    filter(metabolite == metabolite_name)
  
  stat_info <- stats_results %>%
    filter(metabolite == metabolite_name)
  
  clean_title <- str_replace_all(metabolite_name, "_", " ") %>% 
    str_to_title() %>%
    str_trunc(18)
  
  summary_stats <- plot_data %>%
    group_by(treatment) %>%
    summarise(
      n = n(),
      all_zero = all(log_concentration == 0),
      .groups = "drop"
    )
  
  control_all_zero <- summary_stats %>%
    filter(treatment == "non-contaminated") %>%
    pull(all_zero)
  
  treatment_all_zero <- summary_stats %>%
    filter(treatment == "F. verticillioides") %>%
    pull(all_zero)
  
  all_values_zero <- all(plot_data$log_concentration == 0)
  
  control_data <- plot_data %>% filter(treatment == "non-contaminated")
  treatment_data <- plot_data %>% filter(treatment == "F. verticillioides")
  
  min_y <- min(plot_data$log_concentration, na.rm = TRUE)
  max_y <- max(plot_data$log_concentration, na.rm = TRUE)
  
  if(all_values_zero) {
    y_min_plot <- -0.03
    y_max_plot <- 0.08
    y_annotation <- 0.05
  } else if(control_all_zero && !treatment_all_zero) {
    if(max_y <= 1) {
      y_min_plot <- -0.03
      y_max_plot <- max_y + (max_y * 0.1)
      y_annotation <- max_y + (max_y * 0.05)
    } else {
      y_min_plot <- -0.05
      y_max_plot <- max_y + (max_y * 0.08)
      y_annotation <- max_y + (max_y * 0.04)
    }
  } else if(!control_all_zero && treatment_all_zero) {
    if(max_y <= 1) {
      y_min_plot <- -0.03
      y_max_plot <- max_y + (max_y * 0.1)
      y_annotation <- max_y + (max_y * 0.05)
    } else {
      y_min_plot <- -0.05
      y_max_plot <- max_y + (max_y * 0.08)
      y_annotation <- max_y + (max_y * 0.04)
    }
  } else {
    range_y <- max_y - min_y
    y_buffer <- max(0.1, range_y * 0.08)
    y_min_plot <- max(0, min_y - y_buffer)
    y_max_plot <- max_y + y_buffer
    y_annotation <- max_y + (y_buffer * 0.4)
  }
  
  p <- ggplot(plot_data, aes(x = treatment, y = log_concentration, fill = treatment))
  
  if(!control_all_zero) {
    p <- p + geom_violin(data = control_data, alpha = 0.5, trim = FALSE, scale = "width")
  }
  
  if(!treatment_all_zero) {
    p <- p + geom_violin(data = treatment_data, alpha = 0.5, trim = FALSE, scale = "width")
  }
  
  p <- p +
    geom_boxplot(width = 0.25, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.7, size = 1.2, color = "black") +
    scale_fill_manual(values = c("non-contaminated" = "#407126", "F. verticillioides" = "#957e2a")) +
    coord_cartesian(ylim = c(y_min_plot, y_max_plot)) +
    labs(
      title = clean_title,
      x = NULL,
      y = expression(Log[10]~"(Conc. + 1)")
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  if(!is.na(stat_info$p_adj)) {
    sig_label <- if(stat_info$p_adj < 0.001) {
      "p<0.001"
    } else {
      paste0("p=", format(round(stat_info$p_adj, 3), nsmall = 3))
    }
    
    p <- p + annotate("text", x = 1.5, y = y_annotation, label = sig_label, 
                      size = 2, hjust = 0.5, fontface = "bold",
                      color = ifelse(stat_info$p_adj < 0.05, "red", "gray50"))
  }
  
  return(p)
}

# ==============================================================================
# Step8_Batch_Plot_Generation
# ==============================================================================
metabolite_names <- unique(data_clean$metabolite)
dir.create("individual_boxplots", showWarnings = FALSE)
individual_plots <- list()

for(i in seq_along(metabolite_names)) {
  metabolite <- metabolite_names[i]
  tryCatch({
    individual_plot <- create_individual_boxplot(metabolite, save_plot = TRUE)
    individual_plots[[metabolite]] <- individual_plot
  }, error = function(e) {
    next
  })
}

compact_plots <- map(metabolite_names, create_compact_boxplot)
names(compact_plots) <- metabolite_names
p1_combined <- wrap_plots(compact_plots, ncol = 4)

# ==============================================================================
# Step9_Volcano_Plot_Rendering
# ==============================================================================
p2 <- stats_results %>%
  filter(!is.na(p_adj) & !is.na(log2_fc)) %>%
  mutate(
    color_group = case_when(
      p_adj < 0.05 & abs_log2_fc > 1 ~ "Highly significant",
      p_adj < 0.05 ~ "Significant (p < 0.05)",
      abs_log2_fc > 1 ~ "Large change (|FC| > 2)",
      TRUE ~ "Not significant"
    ),
    metabolite_clean = str_replace_all(metabolite, "_", " ") %>% str_to_title()
  ) %>%
  ggplot(aes(x = log2_fc, y = -log10(p_adj))) +
  geom_point(aes(color = color_group), size = 3, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.7) +
  geom_text_repel(
    data = . %>% filter(p_adj < 0.05 & abs_log2_fc > 1),
    aes(label = metabolite_clean),
    size = 3, max.overlaps = 15,
    box.padding = 0.5, point.padding = 0.5,
    segment.color = 'grey50'
  ) +
  scale_color_manual(values = c(
    "Highly significant" = "#E31A1C",
    "Significant (p < 0.05)" = "#FF7F00", 
    "Large change (|FC| > 2)" = "#1F78B4",
    "Not significant" = "gray60"
  )) +
  labs(
    title = "Volcano Plot - Differential Metabolite Analysis",
    subtitle = "F. verticillioides vs Control (log10(concentration + 1))",
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"(Adjusted p-value)"),
    color = "Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

dir.create("metabolomics_plots", showWarnings = FALSE)
ggsave("metabolomics_plots/01_boxplots_combined_conditional_violin.png", p1_combined, width = 16, height = 20, dpi = 300, bg = "white")
ggsave("metabolomics_plots/02_volcano_plot.png", p2, width = 12, height = 10, dpi = 300, bg = "white")

# ==============================================================================
# Step10_PCA_Data_Preparation
# ==============================================================================
# Insert file!!!
data_raw_pca <- read_excel("insert corresponding file", sheet = 1) %>% clean_names()

data_pca_prep <- data_raw_pca %>%
  mutate(
    treatment = factor(treatment, levels = c("non-contaminated", "F. verticillioides")),
    rep = as.numeric(rep),
    sample_id = paste0(treatment, "_", rep)
  ) %>%
  mutate(across(-c(treatment, rep, sample_id), ~ {
    if(is.character(.x)) {
      .x <- gsub(",", ".", .x)
      .x <- gsub("E\\+", "e+", .x)
      .x <- gsub("E-", "e-", .x)
      as.numeric(.x)
    } else {
      as.numeric(.x)
    }
  }))

metabolite_cols_pca <- names(data_pca_prep)[!names(data_pca_prep) %in= c("treatment", "rep", "sample_id")]
data_pca_transformed <- data_pca_prep %>% mutate(across(all_of(metabolite_cols_pca), ~ log10(.x + 1)))

pca_matrix <- data_pca_transformed %>% select(all_of(metabolite_cols_pca)) %>% as.matrix()
rownames(pca_matrix) <- data_pca_transformed$sample_id

var_cols <- apply(pca_matrix, 2, var, na.rm = TRUE)
pca_matrix_clean <- pca_matrix[, var_cols > 0, drop = FALSE]

# ==============================================================================
# Step11_PCA_Execution_and_Outlier_Detection
# ==============================================================================
pca_result <- PCA(pca_matrix_clean, graph = FALSE, scale.unit = TRUE)
group_info <- data_pca_transformed %>%
  select(sample_id, treatment, rep) %>%
  mutate(treatment = factor(treatment, levels = c("non-contaminated", "F. verticillioides"), labels = c("Control", "FV")))

pca_coords <- as.data.frame(pca_result$ind$coord[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  left_join(group_info, by = "sample_id")

detect_outliers <- function(data, threshold = 2.0) {
  group_centers <- data %>%
    group_by(treatment) %>%
    summarise(center_pc1 = mean(Dim.1, na.rm = TRUE), center_pc2 = mean(Dim.2, na.rm = TRUE), .groups = "drop")
  
  data_with_dist <- data %>%
    left_join(group_centers, by = "treatment") %>%
    mutate(dist_to_center = sqrt((Dim.1 - center_pc1)^2 + (Dim.2 - center_pc2)^2))
  
  outlier_thresholds <- data_with_dist %>%
    group_by(treatment) %>%
    summarise(mean_dist = mean(dist_to_center, na.rm = TRUE), sd_dist = sd(dist_to_center, na.rm = TRUE), threshold = mean_dist + threshold * sd_dist, .groups = "drop")
  
  outliers <- data_with_dist %>%
    left_join(outlier_thresholds, by = "treatment") %>%
    mutate(is_outlier = dist_to_center > threshold, outlier_type = case_when(is_outlier ~ "Outlier", TRUE ~ "Normal"))
  
  return(outliers)
}

outlier_analysis <- detect_outliers(pca_coords, threshold = 2.0)
outliers_detected <- outlier_analysis %>% filter(is_outlier)

if(nrow(outliers_detected) > 0) {
  pca_data_final <- outlier_analysis %>% filter(!is_outlier) %>% rename(PC1 = Dim.1, PC2 = Dim.2)
  outlier_samples <- outliers_detected$sample_id
  pca_matrix_corrected <- pca_matrix_clean[!rownames(pca_matrix_clean) %in% outlier_samples, ]
  pca_result_final <- PCA(pca_matrix_corrected, graph = FALSE, scale.unit = TRUE)
  
  pca_data_final <- as.data.frame(pca_result_final$ind$coord[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(group_info, by = "sample_id") %>%
    rename(PC1 = Dim.1, PC2 = Dim.2)
  
  pca_variance <- pca_result_final$eig[1:2, 2]
} else {
  pca_data_final <- pca_coords %>% rename(PC1 = Dim.1, PC2 = Dim.2)
  pca_variance <- pca_result$eig[1:2, 2]
}

# ==============================================================================
# Step12_PCA_Visualization_Rendering
# ==============================================================================
pca_plot <- ggplot(pca_data_final, aes(x = PC1, y = PC2, color = treatment, fill = treatment)) +
  stat_ellipse(aes(fill = treatment, group = treatment), geom = "polygon", type = "norm", level = 0.95, alpha = 0.2, color = NA) +
  stat_ellipse(aes(color = treatment, group = treatment), type = "norm", level = 0.95, linewidth = 0.5, linetype = "solid", alpha = 0.6) +
  geom_point(size = 3, alpha = 0.5, stroke = 0.5, color = "white") +
  geom_point(size = 2.0, alpha = 0.5) +
  stat_summary(aes(group = treatment, color = treatment), fun = mean, geom = "point", size = 1, shape = 16, alpha = 0.6) +
  scale_color_manual(values = c("Control" = "#407126", "FV" = "#957e2a"), name = "Treatment") +
  scale_fill_manual(values = c("Control" = "#407126", "FV" = "#957e2a"), name = "Treatment") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 7, family = "sans", color = "black"),
    axis.title.y = element_text(size = 7, family = "sans", color = "black"),
    axis.text.x = element_text(size = 6, family = "sans", color = "black"),
    axis.text.y = element_text(size = 6, family = "sans", color = "black"),
    legend.title = element_text(size = 7, family = "sans", face = "bold", color = "black"),
    legend.text = element_text(size = 6, family = "sans", color = "black"),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = "NA"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm")
  ) +
  labs(
    x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 1), "%)"),
    title = "Principal Component Analysis (PCA) - Metabolomics",
    subtitle = if(nrow(outliers_detected) > 0) {
      paste0("Outliers removed: ", paste(outliers_detected$sample_id, collapse = ", "), " | Ellipses: 95% confidence intervals")
    } else {
      "Ellipses: 95% confidence intervals"
    }
  )

dir.create("pca_plots", showWarnings = FALSE)
ggsave("pca_plots/pca_microbiome_filled_ellipses.png", pca_plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("pca_plots/pca_microbiome_filled_ellipses.pdf", pca_plot, width = 10, height = 8)

write_csv(pca_data_final, "pca_plots/pca_final_coordinates.csv")
if(nrow(outliers_detected) > 0) {
  write_csv(outliers_detected, "pca_plots/outliers_removed.csv")
}

# ==============================================================================
# Step13_Supplementary_Tables_Export
# ==============================================================================
supplementary_metabolomics_table <- stats_results %>%
  arrange(p_adj) %>%
  dplyr::select(
    Metabolite = metabolite,
    `Mean_Control_Log10` = mean_control,
    `Mean_Infected_Log10` = mean_treatment,
    Log2_FoldChange = log2_fc,
    Raw_p_value = p_value,
    FDR_Adjusted_p = p_adj,
    Significance = significance,
    Statistical_Method = method
  ) %>%
  mutate(
    Metabolite = str_replace_all(Metabolite, "_", " ") %>% 
      str_to_title() %>%
      str_replace("C(\\d+)", "C\\1:") %>%
      str_replace("X7", "7,") %>%
      str_replace("D2", " D2") %>%
      str_replace("E$", " E")
  )

dir.create("metabolomics_plots", showWarnings = FALSE)

wb_metab <- createWorkbook()
addWorksheet(wb_metab, "Metabolomics_Stats")

header_style <- createStyle(
  textDecoration = "bold", 
  halign = "center", 
  fgFill = "#EFEFEF", 
  border = "TopBottomLeftRight"
)

writeData(wb_metab, "Metabolomics_Stats", supplementary_metabolomics_table, headerStyle = header_style)
setColWidths(wb_metab, "Metabolomics_Stats", cols = 1:ncol(supplementary_metabolomics_table), widths = "auto")
saveWorkbook(wb_metab, "metabolomics_plots/Supplementary_Metabolomics_Stats.xlsx", overwrite = TRUE)