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

#

# Author: Wang Yihua


# ============================================================

############################

# 1. Load packages

############################

library(dplyr)

############################

# 2. User-defined parameters

############################

# Working directory

# Please change this to the directory containing:

# tax_L6.txt

# meta_all.txt

data_dir <- "data"

# Input files

genus_file <- file.path(data_dir, "tax_L6.txt")
metadata_file <- file.path(data_dir, "meta_all.txt")

# Output directory

output_dir <- file.path(data_dir, "results", "diff")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################

# 3. Analysis parameters

############################

# Minimum mean relative abundance

# 0.0002 = 0.02%

min_abundance <- 0.0002

# Minimum prevalence

# 0.30 = 30% of samples

min_prevalence <- 0.30

# Number of permutations

n_perm <- 10000

# Significance threshold

fdr_cutoff <- 0.05

# Groups to compare

group_X <- "ADLS"
group_Y <- "HC"

############################

# 4. Read input data

############################

genus <- read.table(
  genus_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)

meta <- read.delim(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################

# 5. Match metadata and abundance data

############################

# Select metadata columns

meta2 <- meta[, c("SampleID", "Project", "Group2")]

# Match metadata to abundance matrix

meta3 <- meta2[match(colnames(genus), meta2$SampleID), ]

# Check whether all samples have metadata

if (any(is.na(meta3$SampleID))) {
  warning("Some samples in the abundance matrix do not have matching metadata.")
}

# Keep samples with valid group information

valid_samples <- !is.na(meta3$Group2)

genus <- genus[, valid_samples, drop = FALSE]
meta3 <- meta3[valid_samples, , drop = FALSE]

############################

# 6. Extract the two groups

############################

group_X_cols <- meta3$SampleID[
  meta3$Group2 == group_X
]

group_Y_cols <- meta3$SampleID[
  meta3$Group2 == group_Y
]

if (length(group_X_cols) == 0) {
  stop(paste("No samples found for group:", group_X))
}

if (length(group_Y_cols) == 0) {
  stop(paste("No samples found for group:", group_Y))
}

genus_X <- genus[, group_X_cols, drop = FALSE]
genus_Y <- genus[, group_Y_cols, drop = FALSE]

cat("Group", group_X, ":", ncol(genus_X), "samples\n")
cat("Group", group_Y, ":", ncol(genus_Y), "samples\n")

############################

# 7. Filter low-abundance genera

############################

# Calculate mean relative abundance

mean_abundance_X <- rowMeans(genus_X, na.rm = TRUE)
mean_abundance_Y <- rowMeans(genus_Y, na.rm = TRUE)

# Calculate prevalence

prevalence_X <- rowSums(genus_X > 0, na.rm = TRUE) /
  ncol(genus_X)

prevalence_Y <- rowSums(genus_Y > 0, na.rm = TRUE) /
  ncol(genus_Y)

# A genus is retained if it satisfies the filtering criteria

# in at least one of the two groups.

keep_genus <- (
  mean_abundance_X >= min_abundance |
    mean_abundance_Y >= min_abundance
) &
  (
    prevalence_X >= min_prevalence |
      prevalence_Y >= min_prevalence
  )

genus_clean <- genus[keep_genus, , drop = FALSE]

cat(
  "Number of genera before filtering:",
  nrow(genus), "\n"
)

cat(
  "Number of genera after filtering:",
  nrow(genus_clean), "\n"
)

############################

# 8. Prepare filtered data

############################

genus_clean_X <- genus_clean[, group_X_cols, drop = FALSE]
genus_clean_Y <- genus_clean[, group_Y_cols, drop = FALSE]

############################

# 9. Permutation test

############################

# Function for a two-sided permutation test

# based on the difference in group means.

permutation_test <- function(x, y, B = 10000) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  observed <- mean(x, na.rm = TRUE) -
    mean(y, na.rm = TRUE)
  
  combined <- c(x, y)
  
  n_x <- length(x)
  
  perm_stats <- numeric(B)
  
  for (b in seq_len(B)) {
    

    permuted <- sample(combined, length(combined))
    
    perm_x <- permuted[seq_len(n_x)]
    perm_y <- permuted[(n_x + 1):length(combined)]
    
    perm_stats[b] <-
      mean(perm_x, na.rm = TRUE) -
      mean(perm_y, na.rm = TRUE)
   
    
  }
  
  # Add one to numerator and denominator
  
  # to avoid obtaining a p-value of zero.
  
  p_value <- (
    sum(abs(perm_stats) >= abs(observed)) + 1
  ) / (B + 1)
  
  return(p_value)
}

############################

# 10. Perform permutation tests

############################

perm_p <- numeric(nrow(genus_clean))

for (i in seq_len(nrow(genus_clean))) {
  
  perm_p[i] <- permutation_test(
    genus_clean_X[i, ],
    genus_clean_Y[i, ],
    B = n_perm
  )
}

perm_p_adj <- p.adjust(
  perm_p,
  method = "BH"
)

perm_results <- data.frame(
  genus = rownames(genus_clean),
  mean_X = rowMeans(genus_clean_X),
  mean_Y = rowMeans(genus_clean_Y),
  prevalence_X = rowSums(genus_clean_X > 0) /
    ncol(genus_clean_X),
  prevalence_Y = rowSums(genus_clean_Y > 0) /
    ncol(genus_clean_Y),
  permutation_p = perm_p,
  permutation_p_adj = perm_p_adj,
  stringsAsFactors = FALSE
)

############################

# 11. Perform Wilcoxon rank-sum tests

############################

wilcox_p <- numeric(nrow(genus_clean))

for (i in seq_len(nrow(genus_clean))) {
  
  wilcox_p[i] <- wilcox.test(
    as.numeric(genus_clean_X[i, ]),
    as.numeric(genus_clean_Y[i, ]),
    alternative = "two.sided",
    exact = FALSE
  )$p.value
}

wilcox_p_adj <- p.adjust(
  wilcox_p,
  method = "BH"
)

wilcox_results <- data.frame(
  genus = rownames(genus_clean),
  mean_X = rowMeans(genus_clean_X),
  mean_Y = rowMeans(genus_clean_Y),
  prevalence_X = rowSums(genus_clean_X > 0) /
    ncol(genus_clean_X),
  prevalence_Y = rowSums(genus_clean_Y > 0) /
    ncol(genus_clean_Y),
  wilcox_p = wilcox_p,
  wilcox_p_adj = wilcox_p_adj,
  stringsAsFactors = FALSE
)

############################

# 12. Merge statistical results

############################

results <- merge(
  perm_results,
  wilcox_results,
  by = c(
    "genus",
    "mean_X",
    "mean_Y",
    "prevalence_X",
    "prevalence_Y"
  )
)

############################

# 13. Add effect direction

############################

results$direction <- ifelse(
  results$mean_X > results$mean_Y,
  paste0(group_X, "_enriched"),
  paste0(group_Y, "_enriched")
)

############################

# 14. Save all results

############################

comparison_name <- paste(
  group_X,
  "vs",
  group_Y,
  sep = "_"
)

# Save complete results

write.csv(
  results,
  file.path(
    output_dir,
    paste0("genus_", comparison_name, "_all.csv")
  ),
  row.names = FALSE
)

# Significant genera based on permutation test

perm_sig <- results %>%
  filter(permutation_p_adj <= fdr_cutoff)

write.csv(
  perm_sig,
  file.path(
    output_dir,
    paste0(
      "genus_",
      comparison_name,
      "_permutation_significant.csv"
    )
  ),
  row.names = FALSE
)

# Significant genera based on Wilcoxon test

wilcox_sig <- results %>%
  filter(wilcox_p_adj <= fdr_cutoff)

write.csv(
  wilcox_sig,
  file.path(
    output_dir,
    paste0(
      "genus_",
      comparison_name,
      "_wilcoxon_significant.csv"
    )
  ),
  row.names = FALSE
)

############################

# 15. Save analysis information

############################

analysis_info <- data.frame(
  parameter = c(
    "Group X",
    "Group Y",
    "Number of samples in Group X",
    "Number of samples in Group Y",
    "Minimum mean abundance",
    "Minimum prevalence",
    "Number of permutations",
    "FDR cutoff",
    "Genera before filtering",
    "Genera after filtering",
    "Permutation significant genera",
    "Wilcoxon significant genera"
  ),
  value = c(
    group_X,
    group_Y,
    ncol(genus_X),
    ncol(genus_Y),
    min_abundance,
    min_prevalence,
    n_perm,
    fdr_cutoff,
    nrow(genus),
    nrow(genus_clean),
    nrow(perm_sig),
    nrow(wilcox_sig)
  )
)

write.csv(
  analysis_info,
  file.path(
    output_dir,
    paste0(
      "genus_",
      comparison_name,
      "_analysis_summary.csv"
    )
  ),
  row.names = FALSE
)

############################

# 16. Finished

############################

cat("\n========================================\n")
cat("Analysis completed successfully.\n")
cat("Comparison:", group_X, "vs", group_Y, "\n")
cat("Results saved to:", output_dir, "\n")
cat("========================================\n")
