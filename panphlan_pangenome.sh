#!/bin/bash
# Step 4: PanPhlAn Analysis for Phocaeicola dorei
# Target: P. dorei pangenome

SAMPLE_ID="Your_Sample_Name"

# 4-1. Map Sample to Pangenome
panphlan_map.py \
    --input results/01_kneaddata/${SAMPLE_ID}/${SAMPLE_ID}.fastq.gz \
    --indexes databases/panphlan_db/Bacteroides_dorei \
    --pangenome databases/panphlan_db/Bacteroides_dorei_pangenome.tsv \
    --output results/04_panphlan/map/${SAMPLE_ID}.tsv \
    --nproc 8

# 4-2. Profiling (Run this ONCE after mapping all samples)
# panphlan_profiling.py \
#    --i_dna results/04_panphlan/map \
#    --pangenome databases/panphlan_db/Bacteroides_dorei_pangenome.tsv \
#    --o_matrix results/04_panphlan/Pdorei_gene_matrix.tsv \
#    --min_coverage 1