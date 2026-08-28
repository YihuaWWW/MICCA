lm_meta.R
# ============================================================

# MMUPHin-based batch effect correction and meta-analysis

# ============================================================

#

# Description:

# This script performs:

# 1. Batch-effect correction using MMUPHin

# 2. Meta-analytical differential abundance analysis

# 3. Differential genus identification for:

# - ADLS vs HC

# - ADNLS vs HC

# - ADLS vs ADNLS

# 4. Visualization of significant genera

#

# Input files:

# data/all/group.csv

# data/all/tax_L6.txt

#

# Required R packages:

# magrittr

# dplyr

# ggplot2

# vegan

# ggpubr

# ggsci

# patchwork

# MMUPHin

#

DA.R
# ============================================================

# Differential genus analysis between two groups

# ============================================================

#

# Description:

# This script identifies differentially abundant genera between

# two groups using permutation tests and Wilcoxon rank-sum tests.

#

# Genera are filtered according to:

# 1. Minimum mean relative abundance

# 2. Minimum prevalence across samples

#

# Input files:

# - tax_L6.txt

# - meta_all.txt

#

# Required R packages:

# - dplyr

lefse.R
# ============================================================

# LEfSe analysis using microeco

# ============================================================

#

# Description:

# This script:

# 1. Matches sample metadata with the abundance table

# 2. Generates study-specific metadata files

# 3. Converts feature abundance to relative abundance

# 4. Performs LEfSe analysis using the microeco package

# 5. Identifies differentially abundant genera

#

# Comparisons:

# - ADLS vs HC

# - ADNLS vs HC

# - ADLS vs ADNLS

#

# Input files:

# data/

# ├── tax_L6.txt

# ├── meta_all.txt

# └── lefse/

# ├── otutable.txt

# └── tax_table.txt

#

# Required R packages:

# tidyverse

# microeco

# magrittr

#

# ============================================================


Cross_cohort_Feature_Selection_and_Classification.R
# ============================================================

# MMUPHin Batch Correction, Meta-analysis, and Random Forest

# ============================================================

#

# Description:

#

# This script performs a multi-cohort microbiome analysis

# consisting of three major steps:

#

# 1. Batch-effect correction using MMUPHin

# 2. Cross-cohort meta-analytical differential abundance

# analysis

# 3. Random Forest classification with feature selection

# using RFE and independent-cohort validation

#

#

# Main comparison:

#

# ADLS vs HC

#

#

# Training cohorts:

#

# S2, S17, S19, S47, S71, S75

#

# Independent validation cohort:

#

# S74

#

#

# Input files:

#

# data/

# └── all/

# ├── group.csv

# └── tax_L6.txt

#

#

# Required R packages:

#

# MMUPHin

# dplyr

# ggplot2

# vegan

# ggpubr

# ggsci

# patchwork

# magrittr

# randomForest

# caret

# pROC

#

#

# Output:

#

# results/

# ├── batch_correction/

# ├── meta_analysis/

# └── machine_learning/

#

#

# Reproducibility:

#

# Random seed: 12

#

# ============================================================

# ============================================================
# Random Forest Classification with Independent Cohort Validation
# ============================================================
#
# Description:
# This script performs:
#
# 1. Identification of robust differential genera across cohorts
# 2. Feature selection based on the intersection of:
#      - Permutation test
#      - Wilcoxon rank-sum test
#      - LEfSe
# 3. Study-stratified 5-fold cross-validation
# 4. Recursive Feature Elimination (RFE)
# 5. Random Forest classification
# 6. Independent validation in a held-out cohort
# 7. Evaluation of classification performance using:
#      - Accuracy
#      - Sensitivity
#      - Specificity
#      - Precision
#      - F1-score
#      - Balanced Accuracy
#      - AUC
# 8. Saving fold assignments and validation results
#
# Comparison:
#   ADLS vs HC
#
# Training cohorts:
#   S2, S17, S19, S47, S71, S75
#
# Independent validation cohort:
#   S74
#
# ============================================================

DA_guided_Feature_Selection_Classification.R
# ============================================================
# 0. Configuration
# ============================================================

# ------------------------------------------------------------
# Project directory
# ------------------------------------------------------------
# IMPORTANT:
# Run this script from the root directory of the GitHub project.
#
# Example project structure:
#
# project/
# ├── data/
# │   ├── all/
# │   │   ├── group.csv
# │   │   └── revised_results/
# │   │       └── genus_NPC_ra_abd_adj.csv
# │   │
# │   ├── 2/
# │   │   └── results/diff/
# │   ├── 17/
# │   │   └── results/diff/
# │   ├── 19/
# │   │   └── results/diff/
# │   ├── 47/
# │   │   └── results/diff/
# │   ├── 71/
# │   │   └── results/diff/
# │   ├── 74/
# │   │   └── results/diff/
# │   └── 75/
# │       └── results/diff/
# │
# └── scripts/
#     └── random_forest_independent_validation.R
#
# If the script is stored in "scripts/", the project root can be
# specified as the parent directory of the current working directory.
#
# For maximum portability, users can also simply run the script
# after setting the project root manually.
#
# ------------------------------------------------------------