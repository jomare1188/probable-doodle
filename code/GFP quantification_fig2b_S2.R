# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(ggpubr)
library(readxl)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)
library(rstatix)
library(openxlsx)
library(broom)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
raw_data <- read_excel("insert corresponding file", col_names = TRUE)

# ==============================================================================
# Step3_Data_Preprocessing
# ==============================================================================
raw_data$Treatment <- as.factor(raw_data$Treatment)
raw_data$Rep <- as.factor(raw_data$Rep)

fv_levels <- c("fv_mm", "fv_mm_sucrose", "fv_mm_saliva")
raw_data$Treatment <- factor(raw_data$Treatment, levels = fv_levels)

fv_data <- raw_data %>% filter(Treatment %in% fv_levels)
ymax <- max(fv_data$GFP, na.rm = TRUE)

# ==============================================================================
# Step4_Initial_Statistical_Modeling
# ==============================================================================
anova_fv <- aov(GFP ~ Treatment, data = fv_data)
emmeans_fv <- emmeans(anova_fv, pairwise ~ Treatment, adjust = "tukey")

comparison_letters <- cld(emmeans_fv$emmeans, Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = str_trim(.group))

# ==============================================================================
# Step5_Data_Visualization
# ==============================================================================
fv_plot <- ggplot(fv_data, aes(x = Treatment, y = GFP, fill = Treatment)) +
  geom_violin(
    position = position_dodge(width = 0.8),
    color = "gray30",
    trim = TRUE,
    alpha = 0.5,
    scale = "width",
    width = 0.7,
    linewidth = 0.2
  ) +
  geom_point(
    aes(fill = Treatment),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    shape = 21,
    color = "black",
    stroke = 0.2,
    size = 0.8,
    alpha = 0.75
  ) +
  scale_fill_manual(
    values = c(
      "fv_mm" = "#FCE698",
      "fv_mm_sucrose" = "#F39C12",
      "fv_mm_saliva" = "#933900"
    ),
    labels = c("MM", "MM + sucrose", "MM + OS"),
    name = "Medium"
  ) +
  scale_x_discrete(
    labels = c("MM", "MM + sucrose", "MM + OS")
  ) +
  labs(
    title = "F. verticillioides - Fungal Biomass (GFP signal)",
    x = "Treatment",
    y = "Integrated Density (a.u.)"
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 6, face = "bold")
  ) +
  geom_text(
    data = comparison_letters,
    aes(x = Treatment, y = ymax * 1.08, label = .group),
    position = position_dodge(width = 0.8),
    size = 2,
    color = "black",
    inherit.aes = FALSE
  ) +
  ylim(0, ymax * 1.15)

print(fv_plot)

# ==============================================================================
# Step6_Supplementary_Data_Compilation
# ==============================================================================
sup_descriptive <- fv_data %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    Mean = mean(GFP, na.rm = TRUE),
    SE = sd(GFP, na.rm = TRUE) / sqrt(n),
    `Mean_SE_Format` = sprintf("%.3f ± %.3f", Mean, SE),
    .groups = "drop"
  )

lm_model_fv <- lm(GFP ~ Treatment, data = fv_data)
shapiro_p <- shapiro.test(residuals(lm_model_fv))$p.value
levene_p <- car::leveneTest(GFP ~ Treatment, data = fv_data)$`Pr(>F)`[1]

sup_assumptions <- data.frame(
  Test = c("Shapiro-Wilk (Normality)", "Levene (Homoscedasticity)"),
  p_value = c(shapiro_p, levene_p),
  Assumption_Met = c(shapiro_p > 0.05, levene_p > 0.05)
)

if (shapiro_p > 0.05 && levene_p > 0.05) {
  anova_res <- aov(GFP ~ Treatment, data = fv_data)
  
  sup_global_test <- tidy(anova_res) %>% 
    mutate(Method = "One-way ANOVA")
  
  tukey_res <- emmeans(anova_res, pairwise ~ Treatment, adjust = "tukey")
  sup_posthoc <- as.data.frame(tukey_res$contrasts) %>%
    rename(Contrast = contrast, Estimate = estimate, Statistic = t.ratio, p_value = p.value) %>%
    mutate(Method = "Tukey HSD")
  
} else {
  kw_res <- kruskal.test(GFP ~ Treatment, data = fv_data)
  
  sup_global_test <- data.frame(
    term = "Treatment",
    statistic = kw_res$statistic,
    df = kw_res$parameter,
    p.value = kw_res$p.value,
    Method = "Kruskal-Wallis"
  )
  
  dunn_res <- rstatix::dunn_test(fv_data, GFP ~ Treatment, p.adjust.method = "bonferroni")
  sup_posthoc <- dunn_res %>%
    mutate(Contrast = paste(group1, "-", group2), Method = "Dunn Test") %>%
    rename(p_value = p.adj, Statistic = statistic) %>%
    dplyr::select(Contrast, Statistic, p_value, Method)
}

anova_fv <- aov(GFP ~ Treatment, data = fv_data)
emmeans_fv <- emmeans(anova_fv, pairwise ~ Treatment, adjust = "tukey")
sup_cld <- cld(emmeans_fv$emmeans, Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = str_trim(.group)) %>%
  rename(Letters = .group, Mean = emmean) %>%
  dplyr::select(Treatment, Mean, SE, df, lower.CL, upper.CL, Letters)

# ==============================================================================
# Step7_Export_Data
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "Supplementary_Stats")

title_style <- createStyle(textDecoration = "bold", fontSize = 12)
current_row <- 1

writeData(wb, "Supplementary_Stats", "1. Descriptive Statistics (Mean ± SE)", startRow = current_row)
addStyle(wb, "Supplementary_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Supplementary_Stats", sup_descriptive, startRow = current_row + 1)
current_row <- current_row + nrow(sup_descriptive) + 3 

writeData(wb, "Supplementary_Stats", "2. Assumptions Verification", startRow = current_row)
addStyle(wb, "Supplementary_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Supplementary_Stats", sup_assumptions, startRow = current_row + 1)
current_row <- current_row + nrow(sup_assumptions) + 3

writeData(wb, "Supplementary_Stats", "3. Global Statistical Test", startRow = current_row)
addStyle(wb, "Supplementary_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Supplementary_Stats", sup_global_test, startRow = current_row + 1)
current_row <- current_row + nrow(sup_global_test) + 3

writeData(wb, "Supplementary_Stats", "4. Post-Hoc Comparisons", startRow = current_row)
addStyle(wb, "Supplementary_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Supplementary_Stats", sup_posthoc, startRow = current_row + 1)
current_row <- current_row + nrow(sup_posthoc) + 3

writeData(wb, "Supplementary_Stats", "5. Tukey HSD - Compact Letter Display (CLD)", startRow = current_row)
addStyle(wb, "Supplementary_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Supplementary_Stats", sup_cld, startRow = current_row + 1)

setColWidths(wb, "Supplementary_Stats", cols = 1:8, widths = "auto")

saveWorkbook(wb, "Supplementary_Statistics_Fig2b.xlsx", overwrite = TRUE)