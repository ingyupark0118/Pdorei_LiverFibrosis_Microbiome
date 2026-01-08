# ==============================================================================
# Multivariate Association with Linear Models (MaAsLin2)
# Purpose: Identification of specific microbial biomarkers (e.g., P. dorei)
#          associated with liver disease, adjusting for confounding factors.
# ==============================================================================

library(Maaslin2)
library(tidyverse)

# 1. Configuration
# ------------------------------------------------------------------------------
# Set your file paths here
INPUT_METADATA  <- "data/metadata_anonymized.tsv" 
INPUT_ABUNDANCE <- "data/species_abundance_table.tsv"
OUTPUT_DIR      <- "results/maaslin2_output"

# Define variables
# Fixed effects: Target variable + Confounders
FIXED_EFFECTS   <- c("Group", "Age", "Sex", "BMI") 
# Random/Batch effects (if applicable, e.g., "Cohort", "SequencingBatch")
RANDOM_EFFECTS  <- c("Cohort", "Probiotic_Use") 
REFERENCE       <- "Group,Normal"

# 2. Data Loading & Preprocessing
# ------------------------------------------------------------------------------
metadata <- read.delim(INPUT_METADATA, sep = "\t", header = TRUE, row.names = 1) %>%
  mutate(
    Group = factor(Group),
    Sex = factor(Sex),
    Cohort = factor(Cohort)
  )

abundance <- read.delim(INPUT_ABUNDANCE, sep = "\t", header = TRUE, row.names = 1)

# 3. Run MaAsLin2
# ------------------------------------------------------------------------------
fit_results <- Maaslin2(
    input_data = abundance,
    input_metadata = metadata,
    output = OUTPUT_DIR,
    fixed_effects = FIXED_EFFECTS,
    random_effects = RANDOM_EFFECTS, # Optional
    reference = REFERENCE,
    normalization = "NONE",          # Assuming data is already normalized/CLR
    transform = "LOG",
    analysis_method = "LM",
    min_prevalence = 0.1,
    standardize = TRUE,
    plot_heatmap = FALSE,            # Basic plots disabled for batch run
    plot_scatter = FALSE,
    cores = 8
)

# 4. Extract Significant Results (e.g., P. dorei)
# ------------------------------------------------------------------------------
target_species <- "Phocaeicola.*dorei" # Regex for target species

significant_results <- fit_results$results %>%
  filter(grepl(target_species, feature, ignore.case = TRUE)) %>%
  arrange(qval)

print("--- MaAsLin2 Analysis Completed ---")
print(head(significant_results))

write.csv(significant_results, file.path(OUTPUT_DIR, "target_species_stats.csv"), row.names = FALSE)