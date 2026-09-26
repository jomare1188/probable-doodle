# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(car)
library(emmeans)
library(broom)
library(openxlsx)
library(readxl)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
raw_data <- read_excel("insert corresponding file")

# ==============================================================================
# Step3_Data_Preprocessing
# ==============================================================================
processed_data <- raw_data %>%
  mutate(
    Treatment = factor(Treatment, levels = c("FV", "FV+OS")),
    Region = factor(Region, levels = c("Up 40-45", "Up 15-20", "Up 0-5", "Down 0-5", "Down 15-20", "Down 40-45")),
    LogAmount = log10(Amount + 1) 
  )

# ==============================================================================
# Step4_Descriptive_Statistics
# ==============================================================================
descriptive_stats <- processed_data %>%
  group_by(Region, Treatment) %>%
  summarise(
    n = n(),
    Mean_Raw = mean(Amount, na.rm = TRUE),
    SE_Raw = sd(Amount, na.rm = TRUE)/sqrt(n),
    Mean_Log = mean(LogAmount, na.rm = TRUE),
    SE_Log = sd(LogAmount, na.rm = TRUE)/sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    `Raw_Mean_±_SE` = sprintf("%.2e ± %.2e", Mean_Raw, SE_Raw),
    `Log_Mean_±_SE` = sprintf("%.3f ± %.3f", Mean_Log, SE_Log)
  )

# ==============================================================================
# Step5_Two_Way_ANOVA_Modeling
# ==============================================================================
anova_model <- aov(LogAmount ~ Treatment * Region, data = processed_data)
anova_results <- tidy(anova_model) %>% 
  mutate(Method = "Two-way ANOVA (Log10 transformed)")

# ==============================================================================
# Step6_Assumption_Testing
# ==============================================================================
shapiro_p <- shapiro.test(residuals(anova_model))$p.value
levene_p <- car::leveneTest(LogAmount ~ Treatment * Region, data = processed_data)$`Pr(>F)`[1]

assumptions_df <- data.frame(
  Test = c("Shapiro-Wilk (Normality of Residuals)", "Levene's Test (Homoscedasticity)"),
  p_value = c(shapiro_p, levene_p),
  Assumption_Met = c(shapiro_p > 0.05, levene_p > 0.05)
)

# ==============================================================================
# Step7_Pairwise_Comparisons
# ==============================================================================
emmeans_res <- emmeans(anova_model, ~ Treatment | Region)
pairwise_comps <- as.data.frame(contrast(emmeans_res, method = "pairwise", adjust = "tukey")) %>%
  mutate(Method = "Tukey HSD")

# ==============================================================================
# Step8_Estimated_Marginal_Means
# ==============================================================================
emmeans_df <- as.data.frame(emmeans_res)
names(emmeans_df)[3] <- "Estimated_Log_Mean"
names(emmeans_df)[4] <- "SE"
names(emmeans_df)[6] <- "Lower_95_CI"
names(emmeans_df)[7] <- "Upper_95_CI"

emmeans_df <- emmeans_df %>%
  dplyr::select(Region, Treatment, Estimated_Log_Mean, SE, Lower_95_CI, Upper_95_CI)

# ==============================================================================
# Step9_Export_Data
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "Spatial_Distribution_Stats")
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

current_row <- 1
writeData(wb, "Spatial_Distribution_Stats", "1. Descriptive Statistics (Raw and Log10 transformed)", startRow = current_row)
addStyle(wb, "Spatial_Distribution_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Spatial_Distribution_Stats", descriptive_stats, startRow = current_row + 1)
current_row <- current_row + nrow(descriptive_stats) + 3

writeData(wb, "Spatial_Distribution_Stats", "2. Assumptions Verification (on Log10 data)", startRow = current_row)
addStyle(wb, "Spatial_Distribution_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Spatial_Distribution_Stats", assumptions_df, startRow = current_row + 1)
current_row <- current_row + nrow(assumptions_df) + 3

writeData(wb, "Spatial_Distribution_Stats", "3. Two-Way ANOVA Results", startRow = current_row)
addStyle(wb, "Spatial_Distribution_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Spatial_Distribution_Stats", anova_results, startRow = current_row + 1)
current_row <- current_row + nrow(anova_results) + 3

writeData(wb, "Spatial_Distribution_Stats", "4. Estimated Marginal Means (95% CI)", startRow = current_row)
addStyle(wb, "Spatial_Distribution_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Spatial_Distribution_Stats", emmeans_df, startRow = current_row + 1)
current_row <- current_row + nrow(emmeans_df) + 3

writeData(wb, "Spatial_Distribution_Stats", "5. Pairwise Comparisons (FV vs FV+OS per Region)", startRow = current_row)
addStyle(wb, "Spatial_Distribution_Stats", title_style, rows = current_row, cols = 1)
writeData(wb, "Spatial_Distribution_Stats", pairwise_comps, startRow = current_row + 1)

setColWidths(wb, "Spatial_Distribution_Stats", cols = 1:10, widths = "auto")
saveWorkbook(wb, "Supplementary_Spatial_Distribution_Stats.xlsx", overwrite = TRUE)