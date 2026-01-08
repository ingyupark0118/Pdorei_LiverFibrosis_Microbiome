# ==============================================================================
# Logistic Regression & ROC Analysis
# Purpose: Evaluate the diagnostic value of P. dorei combined with clinical factors.
# ==============================================================================

library(tidyverse)
library(pROC)
library(broom)
library(ggplot2)

# 1. Data Loading
# ------------------------------------------------------------------------------
# Input should contain: Group, Age, Sex, BMI, and Target_Microbe abundance
data <- read.csv("data/clinical_microbiome_merged.csv", row.names = 1)

# Define Disease Status (0 = Normal, 1 = Disease)
model_data <- data %>%
  mutate(DiseaseStatus = ifelse(Group == "Normal", 0, 1)) %>%
  na.omit()

# 2. Model Construction
# ------------------------------------------------------------------------------
# Model 1: Baseline Clinical Model
model_base <- glm(
  DiseaseStatus ~ Age + Sex + BMI, 
  data = model_data, 
  family = "binomial"
)

# Model 2: Baseline + P. dorei
# Replace 'Abundance_Pdorei' with your actual column name
model_target <- glm(
  DiseaseStatus ~ Age + Sex + BMI + Abundance_Pdorei, 
  data = model_data, 
  family = "binomial"
)

# 3. Odds Ratio Analysis
# ------------------------------------------------------------------------------
or_results <- tidy(model_target, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  select(term, estimate, conf.low, conf.high, p.value)

print("--- Odds Ratios (95% CI) ---")
print(or_results)

# Simple Visualization: Odds Ratios
plot_or <- ggplot(or_results, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Odds Ratio (log scale)", y = "Variables") +
  theme_bw()

ggsave("results/odds_ratio_plot.pdf", plot_or, width = 6, height = 4)

# 4. ROC & DeLong's Test
# ------------------------------------------------------------------------------
roc_base   <- roc(model_data$DiseaseStatus, fitted(model_base))
roc_target <- roc(model_data$DiseaseStatus, fitted(model_target))

# Statistical comparison (DeLong's test)
roc_test <- roc.test(roc_base, roc_target)

print("--- AUROC Comparison ---")
print(paste("Base Model AUC:", round(auc(roc_base), 3)))
print(paste("Target Model AUC:", round(auc(roc_target), 3)))
print(paste("P-value (DeLong):", roc_test$p.value))

# Simple Visualization: ROC Curves
# Using ggroc for basic plotting without custom styling
plot_roc <- ggroc(list(Base = roc_base, Target = roc_target)) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey") +
  labs(title = "ROC Curve Comparison", color = "Model") +
  theme_minimal()

ggsave("results/roc_curve_comparison.pdf", plot_roc, width = 5, height = 5)