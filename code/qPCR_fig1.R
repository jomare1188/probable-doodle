# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(readxl)
library(car)
library(ggpubr)
library(emmeans)
library(gridExtra)
library(multcomp)

# ==============================================================================
# Step2_Data_Import
# ==============================================================================
# Insert file!!!
raw_data <- read_excel("insert corresponding file", col_names = TRUE)

# ==============================================================================
# Step3_Data_Preprocessing
# ==============================================================================
raw_data$Treat <- as.factor(raw_data$Treat)
raw_data$Time <- factor(raw_data$Time, levels = c("0h", "24h", "72h", "120h"))

filtered_data <- raw_data %>%
  filter(Treat %in% c("1", "2"))

treat_order <- c("1", "2")
filtered_data$Treat <- factor(filtered_data$Treat, levels = treat_order)

table(filtered_data$Treat, filtered_data$Time)

# ==============================================================================
# Step4_Assumption_Testing
# ==============================================================================
shapiro_test <- shapiro.test(filtered_data$Value)
print(shapiro_test)

levene_test <- leveneTest(Value ~ Treat * Time, data = filtered_data)
print(levene_test)

# ==============================================================================
# Step5_Statistical_Analysis
# ==============================================================================
anova_result <- aov(Value ~ Treat * Time, data = filtered_data)
summary(anova_result)

emmeans_result <- emmeans(anova_result, pairwise ~ Treat | Time)
print(emmeans_result)

comparisons_result <- emmeans(anova_result, pairwise ~ Treat | Time, adjust = "bonferroni")
comparisons_df <- as.data.frame(comparisons_result$contrasts)

comparison_letters <- cld(comparisons_result, Letters = letters)

# ==============================================================================
# Step6_Data_Visualization
# ==============================================================================
combined_plot <- ggplot(filtered_data, aes(x = Time, y = Value, fill = Treat)) +
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
    aes(fill = Treat),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    shape = 21,
    color = "black",
    stroke = 0.2,
    size = 0.8,
    alpha = 0.75
  ) +
  scale_fill_manual(
    values = c(
      "#d59500",
      "#63072c"
    ),
    name = "Treatment",
    labels = c("Treatment 2", "Treatment 4")
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    legend.title = element_text(size = 4),
    legend.text = element_text(size = 4),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  stat_compare_means(
    aes(group = Treat), 
    method = "t.test",
    label = "p.signif",
    label.y = max(filtered_data$Value, na.rm = TRUE) * 1.15,
    size = 2
  ) +
  geom_text(
    data = comparison_letters, 
    aes(x = Time, y = max(filtered_data$Value, na.rm = TRUE) * 1.08, label = .group), 
    position = position_dodge(width = 0.8),
    size = 2,
    color = "black"
  )

print(combined_plot)

# ==============================================================================
# Step7_Post_Hoc_T_Tests
# ==============================================================================
for(current_time in levels(filtered_data$Time)) {
  time_data <- filtered_data[filtered_data$Time == current_time, ]
  if(nrow(time_data) > 3) {
    test_result <- t.test(Value ~ Treat, data = time_data)
    cat(sprintf("Time %s: p-value = %.4f\n", current_time, test_result$p.value))
  }
}

# ==============================================================================
# Step8_Export_Graphics
# ==============================================================================
ggsave("qpcr_graph_treatment2_4.tiff", plot = combined_plot, dpi = 600, width = 10, height = 4, units = "cm", compression = "lzw")

# ==============================================================================
# Step9_Descriptive_Statistics
# ==============================================================================
descriptive_stats <- filtered_data %>%
  group_by(Treat, Time) %>%
  summarise(
    n = n(),
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    median_val = median(Value, na.rm = TRUE),
    .groups = 'drop'
  )

print(descriptive_stats)
