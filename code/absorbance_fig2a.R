# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(readxl)
library(ggpubr)
library(reshape2)
library(car)
library(emmeans)
library(rstatix)
library(writexl)
library(broom)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
raw_data <- read_excel("insert corresponding file", col_names = TRUE)
names(raw_data)[1:3] <- c("Code", "Description", "Replicate")

# ==============================================================================
# Step3_Data_Transformation
# ==============================================================================
long_data <- melt(
  raw_data,
  id.vars = c("Code", "Description", "Replicate"),
  variable.name = "Time",
  value.name = "Absorbance"
)

long_data$Time <- factor(
  long_data$Time,
  levels = c("0h", "24h", "48h", "72h", "96h", "120h")
)

# ==============================================================================
# Step4_Baseline_Correction
# ==============================================================================
od0 <- long_data %>% 
  filter(Time == "0h") %>% 
  dplyr::select(Code, Replicate, Absorbance_0h = Absorbance)

corrected_data <- long_data %>% 
  left_join(od0, by = c("Code", "Replicate")) %>% 
  mutate(Abs_corrected = as.numeric(Absorbance) - as.numeric(Absorbance_0h))

# ==============================================================================
# Step5_Group_Segmentation
# ==============================================================================
fv <- corrected_data %>% filter(str_detect(Description, "^Fv -"))
ctrl <- corrected_data %>% filter(str_detect(Description, "Controle"))

# ==============================================================================
# Step6_Visualization_Function
# ==============================================================================
plot_group_data <- function(df, plot_title, color_palette) {
  summary_df <- df %>% 
    group_by(Description, Time) %>% 
    summarise(
      mean_val = mean(Abs_corrected, na.rm = TRUE),
      error_val = sd(Abs_corrected, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    mutate(Shape_Type = case_when(
      str_detect(Description, "saliva") ~ "triangle",
      str_detect(Description, "sacarose") ~ "square",
      TRUE ~ "circle"
    ))
  
  ggplot(summary_df, aes(x = Time, y = mean_val, group = Description, color = Description, shape = Shape_Type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_val - error_val, ymax = mean_val + error_val), width = 0.15) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = c("circle" = 16, "square" = 15, "triangle" = 17)) +
    labs(
      title = plot_title,
      x = "Incubation Time (hours)",
      y = "Optical Density at 600 nm",
      color = "Treatment"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    ylim(0, NA)
}

# ==============================================================================
# Step7_Generate_Palettes
# ==============================================================================
colors_fv <- scales::seq_gradient_pal("#F5B041", "#957e2a", "Lab")(seq(0, 1, length.out = length(unique(fv$Description))))
colors_ctrl <- scales::seq_gradient_pal("#BDC3C7", "#2C3E50", "Lab")(seq(0, 1, length.out = length(unique(ctrl$Description))))

# ==============================================================================
# Step8_Render_Individual_Plots
# ==============================================================================
g1 <- plot_group_data(fv, "F. verticillioides", colors_fv)
g2 <- plot_group_data(ctrl, "Control", colors_ctrl)

ymax <- max(c(layer_scales(g1)$y$range$range[2],
              layer_scales(g2)$y$range$range[2]))

g1 <- g1 + ylim(0, ymax)
g2 <- g2 + ylim(0, ymax)

# ==============================================================================
# Step9_Combine_Plots
# ==============================================================================
panel_output <- ggarrange(
  g1, g2,
  ncol = 2, nrow = 1,
  labels = c("A", "B"),
  font.label = list(size = 12, face = "bold"),
  common.legend = FALSE
)

print(panel_output)

# ==============================================================================
# Step10_Descriptive_Statistics
# ==============================================================================
descriptive_table <- corrected_data %>%
  group_by(Time, Description) %>%
  summarise(
    N = n(),
    Mean = mean(Abs_corrected, na.rm = TRUE),
    SE = sd(Abs_corrected, na.rm = TRUE) / sqrt(N),
    Mean_SE = sprintf("%.3f ± %.3f", Mean, SE),
    .groups = "drop"
  ) %>%
  rename(
    Treatment_Description = Description
  )

# ==============================================================================
# Step11_Inferential_Analysis
# ==============================================================================
stat_data <- corrected_data %>% filter(Time != "0h")
time_points <- unique(stat_data$Time)

tests_list <- list()
posthoc_list <- list()

for (t in time_points) {
  df_temp <- stat_data %>% filter(Time == t)
  
  lm_model <- lm(Abs_corrected ~ Description, data = df_temp)
  shapiro_p <- shapiro.test(residuals(lm_model))$p.value
  
  levene_p <- car::leveneTest(Abs_corrected ~ Description, data = df_temp)$`Pr(>F)`[1]
  
  if (shapiro_p > 0.05 && levene_p > 0.05) {
    aov_model <- aov(Abs_corrected ~ Description, data = df_temp)
    
    test_res <- tidy(aov_model) %>% 
      filter(term == "Description") %>%
      mutate(
        Time = t, 
        Method = "One-way ANOVA", 
        Shapiro_p = shapiro_p, 
        Levene_p = levene_p
      ) %>%
      rename(Term = term, Statistic = statistic)
    
    tukey_res <- emmeans(aov_model, pairwise ~ Description, adjust = "tukey")
    ph_res <- as.data.frame(tukey_res$contrasts) %>%
      mutate(
        Time = t, 
        PostHoc_Method = "Tukey HSD"
      ) %>%
      rename(Contrast = contrast, Estimate = estimate, Statistic = t.ratio)
    
  } else {
    kw_res <- kruskal.test(Abs_corrected ~ Description, data = df_temp)
    
    test_res <- data.frame(
      Time = t, 
      Term = "Treatment_Description", 
      df = kw_res$parameter, 
      Statistic = kw_res$statistic, 
      p.value = kw_res$p.value,
      Method = "Kruskal-Wallis", 
      Shapiro_p = shapiro_p, 
      Levene_p = levene_p
    )
    
    dunn_res <- rstatix::dunn_test(df_temp, Abs_corrected ~ Description, p.adjust.method = "bonferroni")
    ph_res <- dunn_res %>%
      mutate(
        Contrast = paste(group1, "-", group2),
        Time = t, 
        PostHoc_Method = "Dunn"
      ) %>%
      rename(p.value = p.adj, Statistic = statistic) %>%
      dplyr::select(Time, Contrast, Statistic, p.value, PostHoc_Method)
  }
  
  tests_list[[as.character(t)]] <- test_res
  posthoc_list[[as.character(t)]] <- ph_res
}

# ==============================================================================
# Step12_Compile_Results
# ==============================================================================
final_tests_table <- bind_rows(tests_list) %>%
  dplyr::select(Time, Method, Shapiro_p, Levene_p, Term, df, Statistic, p.value)

final_posthoc_table <- bind_rows(posthoc_list) %>%
  dplyr::select(Time, PostHoc_Method, Contrast, Statistic, p.value, everything()) %>%
  dplyr::select(-any_of(c("SE", "df", "group1", "group2", "n1", "n2", "p", "p.adj.signif")))

# ==============================================================================
# Step13_Export_Data
# ==============================================================================
export_sheets <- list(
  "1_Descriptive_Stats" = descriptive_table,
  "2_Assumptions_and_Tests" = final_tests_table,
  "3_PostHoc_Comparisons" = final_posthoc_table
)

# Insert file: [Statistical_Analysis_Results.xlsx]
write_xlsx(export_sheets, "Statistical_Analysis_Results.xlsx")