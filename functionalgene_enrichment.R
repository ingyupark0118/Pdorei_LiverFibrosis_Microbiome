#!/usr/bin/env Rscript

# ==============================================================================
# Functional Enrichment Analysis (KEGG Modules)
# Purpose: Map PanPhlAn gene families to KEGG Modules and perform Fisher's exact test
#          to identify differentially present functional modules between groups.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

rm(list = ls())

# ==============================================================================
# 1. Configuration & Data Loading
# ==============================================================================
# [USER SETTING] Define paths and column names here
BASE_DIR      <- "."                       # Current directory
INPUT_PA      <- "data/gene_presence_absence.tsv"
INPUT_META    <- "data/metadata.tsv"
INPUT_ANNOT   <- "ref/eggnog_annotations.tsv"
INPUT_MAP     <- "ref/ko_to_module.tsv"

OUTPUT_DIR    <- "results/functional_analysis"
GROUP_COL     <- "Group"                   # Column name for grouping (e.g., Disease)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1) Load PanPhlAn Matrix (Gene x Sample)
#    Assumes the first column is the gene ID (UniRef90)
pa_mat <- fread(file.path(BASE_DIR, INPUT_PA), header = TRUE) %>%
  column_to_rownames(names(.)[1])
  
# Clean up sample names (optional, customize as needed)
colnames(pa_mat) <- gsub("^panphlan_result_|\\.tsv$", "", colnames(pa_mat))

# 2) Load Metadata
meta_df <- fread(file.path(BASE_DIR, INPUT_META)) %>%
  column_to_rownames("sample_id") %>%
  filter(row.names(.) %in% colnames(pa_mat))

# Align matrix columns with metadata rows
pa_mat <- pa_mat[, rownames(meta_df), drop = FALSE]

print(paste("Loaded:", ncol(pa_mat), "samples and", nrow(pa_mat), "genes."))

# ==============================================================================
# 2. Functional Mapping: Gene -> KO -> KEGG Module
# ==============================================================================
print("Mapping Genes to KEGG Modules...")

# 1) Load Annotation (UniRef90 -> KO)
#    Parses eggNOG-mapper output
gene_to_ko <- fread(file.path(BASE_DIR, INPUT_ANNOT), skip = 4, fill = TRUE) %>%
  select(query, KEGG_ko) %>%
  mutate(
    uniref = sub("::.*$", "", query),
    ko     = gsub("ko:", "", KEGG_ko)
  ) %>%
  filter(ko != "" & uniref %in% rownames(pa_mat)) %>%
  separate_rows(ko, sep = ",") %>%
  distinct(uniref, ko)

# 2) Load Mapping File (KO -> Module)
ko_to_mod <- fread(file.path(BASE_DIR, INPUT_MAP)) %>%
  mutate(
    ko     = gsub("ko:", "", ko), 
    module = gsub("md:", "", module)
  )

# 3) Create Final Mapping Table (Gene -> Module)
#    Aggregates genes into functional modules
gene_to_mod_map <- gene_to_ko %>%
  inner_join(ko_to_mod, by = "ko") %>%
  select(uniref, module) %>%
  distinct()

# 4) Convert Gene Matrix to Module Matrix (Aggregation)
#    Logic: If any gene in a module is present (1), the module is considered present.
module_pa_mat <- pa_mat %>%
  rownames_to_column("uniref") %>%
  inner_join(gene_to_mod_map, by = "uniref") %>%
  group_by(module) %>%
  summarise(across(all_of(rownames(meta_df)), max)) %>% # Max(0,1) = Presence/Absence
  column_to_rownames("module")

print(paste("Aggregated into", nrow(module_pa_mat), "KEGG Modules."))

# ==============================================================================
# 3. Statistical Analysis (Pairwise Fisher's Exact Test)
# ==============================================================================
print("Running Pairwise Fisher Tests...")

# Function: Perform Fisher's test for a specific pair of groups
run_pairwise_fisher <- function(matrix, meta, group_col, g1, g2) {
  
  # Subset data for the two groups
  target_samples <- rownames(meta)[meta[[group_col]] %in% c(g1, g2)]
  sub_mat  <- matrix[, target_samples, drop=FALSE]
  sub_meta <- meta[target_samples, group_col]
  
  # Apply Fisher test for each module
  results <- apply(sub_mat, 1, function(row) {
    # Create 2x2 contingency table (Presence/Absence vs Group1/Group2)
    tbl <- table(factor(row, levels = c(0, 1)), 
                 factor(sub_meta, levels = c(g1, g2)))
    fisher.test(tbl)$p.value
  })
  
  data.frame(
    Module  = names(results),
    Group1  = g1,
    Group2  = g2,
    P_value = as.numeric(results)
  )
}

# Automatically detect all unique groups and generate combinations
unique_groups <- unique(meta_df[[GROUP_COL]])
comparisons   <- combn(unique_groups, 2, simplify = FALSE)

# Execute tests for all combinations
all_stats <- map_df(comparisons, ~run_pairwise_fisher(
  module_pa_mat, meta_df, GROUP_COL, .x[1], .x[2]
))

# ==============================================================================
# 4. Result Filtering & Saving
# ==============================================================================

# Filter significant results (FDR < 0.05)
sig_results <- all_stats %>%
  group_by(Group1, Group2) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  filter(P_value < 0.05) %>%  # Filter by raw P-value or FDR as needed
  arrange(P_value)

# Display & Save
print(head(sig_results))

write.csv(sig_results, 
          file.path(OUTPUT_DIR, "significant_kegg_modules.csv"), 
          row.names = FALSE)

print("Analysis Completed. Results saved.")