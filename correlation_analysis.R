# ==============================================================================
# Spearman Correlation Analysis
# Purpose: Correlation between P. dorei abundance and clinical markers
# ==============================================================================

library(tidyverse)
library(ggplot2)

# 1. Load Data
# ------------------------------------------------------------------------------
# Input data should contain the target microbe and clinical variable columns
data <- read.csv("data/clinical_microbiome_merged.csv")

# 2. Statistical Test (Spearman)
# ------------------------------------------------------------------------------
# Replace 'p_dorei' and 'PLATELET' with your actual column names
cor_result <- cor.test(data$p_dorei, data$PLATELET, method = "spearman")

# Create a label for the plot
label_text <- paste0("Spearman rho = ", round(cor_result$estimate, 2), 
                     "\nP-value = ", signif(cor_result$p.value, 3))

print("--- Correlation Result ---")
print(label_text)

# 3. Visualization (Standard Scatter Plot)
# ------------------------------------------------------------------------------
plot_corr <- ggplot(data, aes(x = p_dorei, y = PLATELET)) +
  geom_point(alpha = 0.6) +                  # Basic points
  geom_smooth(method = "lm", se = TRUE) +    # Linear regression line with CI
  annotate("text", x = Inf, y = Inf, label = label_text, 
           hjust = 1.1, vjust = 1.5, size = 5) + # Add stats on plot
  labs(
    x = "Relative Abundance (Log10)",
    y = "Clinical Marker Level"
  ) +
  theme_bw()                                 # Standard clean theme

# Save plot
ggsave("results/correlation_plot_basic.pdf", plot_corr, width = 5, height = 5)