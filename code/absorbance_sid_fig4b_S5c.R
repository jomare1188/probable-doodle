# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(readxl)
library(reshape2)
library(broom)
library(emmeans)
library(multcomp)
library(openxlsx)
library(rstatix)

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
  levels = c("0h", "1 day", "2 day", "3 day", "4 day", "5 day", "6 day", "7 day")
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
# Step5_Feature_Engineering
# ==============================================================================
time_numeric_map <- c("0h" = 0, "1 day" = 24, "2 day" = 48, "3 day" = 72, 
                      "4 day" = 96, "5 day" = 120, "6 day" = 144, "7 day" = 168)

corrected_data$Time_numeric <- as.numeric(time_numeric_map[as.character(corrected_data$Time)])

# ==============================================================================
# Step6_AUC_Calculation_Function
# ==============================================================================
calc_auc <- function(time_vec, abs_vec) {
  sorting_order <- order(time_vec)
  time_sorted <- time_vec[sorting_order]
  abs_sorted <- abs_vec[sorting_order]
  auc_val <- sum(diff(time_sorted) * (head(abs_sorted, -1) + tail(abs_sorted, -1)) / 2)
  return(auc_val)
}

# ==============================================================================
# Step7_Compute_Metrics
# ==============================================================================
auc_df <- corrected_data %>%
  group_by(Description, Replicate) %>%
  summarise(
    AUC = calc_auc(Time_numeric, Abs_corrected),
    .groups = 'drop'
  ) %>%
  rename(Treatment = Description)

auc_summary <- auc_df %>%
  group_by(Treatment) %>%
  summarise(
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    AUC_se = sd(AUC, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = 'drop'
  ) %>%
  arrange(desc(AUC_mean))

plot_data <- corrected_data %>%
  group_by(Description, Time) %>%
  summarise(
    mean_val = mean(Abs_corrected, na.rm = TRUE),
    error_val = sd(Abs_corrected, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  rename(Treatment = Description)

# ==============================================================================
# Step8_AUC_Statistical_Modeling
# ==============================================================================
anova_auc <- aov(AUC ~ Treatment, data = auc_df)
emm_auc <- emmeans(anova_auc, ~ Treatment)
cld_auc <- cld(emm_auc, Letters = letters, adjust = 'tukey')

auc_summary_letters <- as.data.frame(cld_auc) %>%
  dplyr::select(Treatment, emmean, SE, .group) %>%
  rename(AUC_mean = emmean, AUC_se = SE, Letters = .group) %>%
  mutate(Letters = str_trim(Letters)) %>%
  arrange(desc(AUC_mean))

# ==============================================================================
# Step9_Color_Palette_Config
# ==============================================================================
vibrant_colors <- c(
  "OS filtrated" = "#00E676",             
  "OS filt. + EDTA" = "#2196F3",         
  "OS filt. + EDTA + Fe" = "#FFEA99",     
  "OS filt. + Sid." = "#FF9199",            
  "OS filt. + Sid. + Fe" = "#000080",      
  "OS filt. + Sid. + EDTA + Fe" = "#FFA366", 
  "Sid. + EDTA + Fe" = "#00BCD4",          
  "OS completed 20%" = "#FF1744"            
)

# ==============================================================================
# Step10_Plot_Growth_Curves
# ==============================================================================
p1 <- ggplot(plot_data, aes(x = Time, y = mean_val, color = Treatment, group = Treatment)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.7, fill = "white", stroke = 1) +
  geom_errorbar(aes(ymin = mean_val - error_val, ymax = mean_val + error_val), 
                width = 0.1, linewidth = 0.2, color = "#565656") +
  scale_color_manual(values = vibrant_colors) +
  labs(
    x = "Incubation time (h)",
    y = expression(Delta*"OD"["600"]),
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.margin = margin(5, 5, 5, 5, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

png("Figure_1_Growth_Curves.png", width = 180, height = 100, units = "mm", res = 600)
print(p1)
dev.off()

# ==============================================================================
# Step11_Plot_AUC_Bar_Chart
# ==============================================================================
p2 <- ggplot(auc_summary_letters, aes(x = reorder(Treatment, AUC_mean), y = AUC_mean, fill = Treatment)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.7) +
  geom_errorbar(aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
                width = 0.2, linewidth = 0.4, color = "#565656") +
  geom_text(aes(y = AUC_mean + AUC_se + (max(AUC_mean + AUC_se) * 0.05), label = Letters),
            size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(values = vibrant_colors) +
  labs(
    x = NULL,
    y = expression("Area under the curve ("*Delta*"OD"["600"]*" × h)"),
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.margin = margin(5, 5, 5, 5, "mm"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.title = element_text(color = "black", size = 10, face = "bold"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "#E0E0E0", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

png("Figure_2_AUC_Treatments.png", width = 180, height = 100, units = "mm", res = 600)
print(p2)
dev.off()

# ==============================================================================
# Step12_Supplementary_Data_Compilation
# ==============================================================================
sup_descriptive <- auc_df %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    Mean = mean(AUC, na.rm = TRUE),
    SE = sd(AUC, na.rm = TRUE) / sqrt(n),
    `Mean_SE_Format` = sprintf("%.3f ± %.3f", Mean, SE),
    .groups = "drop"
  )

lm_model_auc <- lm(AUC ~ Treatment, data = auc_df)
shapiro_p <- shapiro.test(residuals(lm_model_auc))$p.value
levene_p <- car::leveneTest(AUC ~ Treatment, data = auc_df)$`Pr(>F)`[1]

sup_assumptions <- data.frame(
  Test = c("Shapiro-Wilk (Normality)", "Levene (Homoscedasticity)"),
  p_value = c(shapiro_p, levene_p),
  Assumption_Met = c(shapiro_p > 0.05, levene_p > 0.05)
)

if (shapiro_p > 0.05 && levene_p > 0.05) {
  sup_global_test <- tidy(anova_auc) %>% mutate(Method = "One-way ANOVA")
  
  sup_posthoc <- as.data.frame(pairs(emm_auc)) %>%
    rename(Contrast = contrast, Estimate = estimate, Statistic = t.ratio, p_value = p.value) %>%
    mutate(Method = "Tukey HSD")
} else {
  kw_res <- kruskal.test(AUC ~ Treatment, data = auc_df)
  sup_global_test <- data.frame(
    term = "Treatment", statistic = kw_res$statistic, df = kw_res$parameter,
    p.value = kw_res$p.value, Method = "Kruskal-Wallis"
  )
  
  dunn_res <- rstatix::dunn_test(auc_df, AUC ~ Treatment, p.adjust.method = "bonferroni")
  sup_posthoc <- dunn_res %>%
    mutate(Contrast = paste(group1, "-", group2), Method = "Dunn Test") %>%
    rename(p_value = p.adj, Statistic = statistic) %>%
    dplyr::select(Contrast, Statistic, p_value, Method)
}

sup_cld <- auc_summary_letters %>%
  rename(Mean = AUC_mean, SE = AUC_se) %>%
  dplyr::select(Treatment, Mean, SE, Letters)

# ==============================================================================
# Step13_Export_Data
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "Siderophore_Kinetics_Stats")
title_style <- createStyle(textDecoration = "bold", fontSize = 12)
current_row <- 1

writeData(wb, "Siderophore_Kinetics_Stats", "1. Descriptive Statistics of AUC (Mean ± SE)", startRow = current_row)
addStyle(wb, "Siderophore_Kinetics_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Siderophore_Kinetics_Stats", sup_descriptive, startRow = current_row + 1)
current_row <- current_row + nrow(sup_descriptive) + 3 

writeData(wb, "Siderophore_Kinetics_Stats", "2. Assumptions Verification", startRow = current_row)
addStyle(wb, "Siderophore_Kinetics_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Siderophore_Kinetics_Stats", sup_assumptions, startRow = current_row + 1)
current_row <- current_row + nrow(sup_assumptions) + 3

writeData(wb, "Siderophore_Kinetics_Stats", "3. Global Statistical Test", startRow = current_row)
addStyle(wb, "Siderophore_Kinetics_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Siderophore_Kinetics_Stats", sup_global_test, startRow = current_row + 1)
current_row <- current_row + nrow(sup_global_test) + 3

writeData(wb, "Siderophore_Kinetics_Stats", "4. Post-Hoc Comparisons", startRow = current_row)
addStyle(wb, "Siderophore_Kinetics_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Siderophore_Kinetics_Stats", sup_posthoc, startRow = current_row + 1)
current_row <- current_row + nrow(sup_posthoc) + 3

writeData(wb, "Siderophore_Kinetics_Stats", "5. Compact Letter Display (CLD)", startRow = current_row)
addStyle(wb, "Siderophore_Kinetics_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Siderophore_Kinetics_Stats", sup_cld, startRow = current_row + 1)

setColWidths(wb, "Siderophore_Kinetics_Stats", cols = 1:8, widths = "auto")

saveWorkbook(wb, "Supplementary_Siderophore_Kinetics_Stats.xlsx", overwrite = TRUE)