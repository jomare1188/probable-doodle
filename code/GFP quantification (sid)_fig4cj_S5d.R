# ==============================================================================
# Step1_Load_Packages
# ==============================================================================
library(tidyverse)
library(ggpubr)
library(readxl)
library(car)
library(emmeans)
library(multcomp)

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

sid_levels <- c(
  "Negative control", "OS filtrated+EDTA+Fe", "OS filtrated+EDTA", "OS filtrated+EDTA+Fe+Desferrioxamine E",
  "OS filtrated+Desferrioxamine E", "OS filtrated", "OS filtrated+Desferrioxamine E+Fe", "OS"
)

raw_data$Treatment <- factor(raw_data$Treatment, levels = sid_levels)

# ==============================================================================
# Step4_Global_Scale_Setup
# ==============================================================================
ymax <- max(raw_data$GFP, na.rm = TRUE)

# ==============================================================================
# Step5_Group_Filtering
# ==============================================================================
sid_data <- raw_data %>% filter(Treatment %in% sid_levels)
sid_data$Treatment <- factor(sid_data$Treatment, levels = sid_levels)

# ==============================================================================
# Step6_Data_Visualization
# ==============================================================================
sid_plot <- ggplot(sid_data, aes(x = Treatment, y = GFP, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, color = "black", width = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.25) +
  scale_fill_manual(
    values = c(
      "Negative control" = "#00BCD4",
      "OS filtrated+EDTA+Fe" = "#FFEA99",
      "OS filtrated+EDTA" = "#2196F3",
      "OS filtrated+EDTA+Fe+Desferrioxamine E" = "#FFA366",
      "OS filtrated+Desferrioxamine E" = "#FF9199",
      "OS filtrated" = "#00E676",
      "OS filtrated+Desferrioxamine E+Fe" = "#000080",
      "OS" = "#FF1744"
    ),
    labels = c("C-", "OS filt.+EDTA+Fe", "OS filt.+EDTA", "OS filt.+EDTA+Fe+sid", "OS filt.+sid", "OS filt.", "OS filt.+sid+Fe", "OS"),
    name = "Treat"
  ) +
  scale_x_discrete(
    labels = c("C-", "OS filt.+EDTA+Fe", "OS filt.+EDTA", "OS filt.+EDTA+Fe+sid", "OS filt.+sid", "OS filt.", "OS filt.+sid+Fe", "OS")
  ) +
  labs(
    title = "F. verticillioides - Fungal Biomass (GFP signal)",
    x = "Treatment",
    y = "GFP Intensity"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  ylim(0, ymax)

print(sid_plot)

# ==============================================================================
# Step7_Statistical_Analysis
# ==============================================================================
fv_anova <- aov(GFP ~ Treatment, data = sid_data)
summary(fv_anova)

fv_emmeans <- emmeans(fv_anova, pairwise ~ Treatment, adjust = "tukey")
cld(fv_emmeans, Letters = letters)