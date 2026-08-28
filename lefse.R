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

############################

# 1. Install and load packages

############################

cran_packages <- c(
  "tidyverse",
  "magrittr"
)

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
  library(
    pkg,
    character.only = TRUE,
    quietly = TRUE,
    warn.conflicts = FALSE
  )
}

# microeco is installed from CRAN

if (!requireNamespace("microeco", quietly = TRUE)) {
  install.packages("microeco")
}

suppressPackageStartupMessages(
  library(microeco)
)

############################

# 2. User-defined parameters

############################

# ------------------------------------------------------------

# Project directory

#

# Recommended repository structure:

#

# repository/

# ├── data/

# │   ├── tax_L6.txt

# │   ├── meta_all.txt

# │   └── lefse/

# │       ├── otutable.txt

# │       └── tax_table.txt

# │

# ├── R/

# │   └── lefse_microeco.R

# │

# └── results/

# └── diff/

# ------------------------------------------------------------

project_dir <- "."

data_dir <- file.path(
  project_dir,
  "data"
)

lefse_dir <- file.path(
  data_dir,
  "lefse"
)

result_dir <- file.path(
  project_dir,
  "results",
  "diff"
)

############################

# 3. Create output directory

############################

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################

# 4. Input files

############################

genus_file <- file.path(
  data_dir,
  "tax_L6.txt"
)

metadata_file <- file.path(
  data_dir,
  "meta_all.txt"
)

feature_file <- file.path(
  lefse_dir,
  "otutable.txt"
)

tax_file <- file.path(
  lefse_dir,
  "tax_table.txt"
)

############################

# 5. Analysis parameters

############################

# LEfSe significance threshold

lefse_alpha <- 0.05

# Minimum LDA score

lda_cutoff <- 2

# Minimum number of subsamples

lefse_min_subsam <- 8

############################

# 6. Read genus abundance data

############################

genus <- read.table(
  genus_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

############################

# 7. Read metadata

############################

meta.all <- read.delim(
  metadata_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################

# 8. Match metadata to genus abundance data

############################

if (!"SampleID" %in% colnames(meta.all)) {
  
  stop(
    "The metadata file must contain a column named 'SampleID'."
  )
}

matched_samples <- intersect(
  colnames(genus),
  meta.all$SampleID
)

if (length(matched_samples) == 0) {
  
  stop(
    "No matching samples were found between abundance data and metadata."
  )
}

# Keep only matched samples

genus <- genus[
  ,
  matched_samples,
  drop = FALSE
]

meta_study <- meta.all[
  match(
    matched_samples,
    meta.all$SampleID
  ),
  ,
  drop = FALSE
]

############################

# 9. Save matched metadata

############################

write.table(
  meta_study,
  file = file.path(
    data_dir,
    "meta_study.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

write.csv(
  meta_study,
  file = file.path(
    data_dir,
    "meta_study.csv"
  ),
  row.names = FALSE
)

############################

# 10. Read LEfSe input data

############################

feature_table <- read.table(
  feature_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

tax_table <- read.table(
  tax_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t"
)

############################

# 11. Match features and taxonomy

############################

common_features <- intersect(
  rownames(feature_table),
  rownames(tax_table)
)

if (length(common_features) == 0) {
  
  stop(
    "No common features were found between the feature table and taxonomy table."
  )
}

feature_table <- feature_table[
  common_features,
  ,
  drop = FALSE
]

tax_table <- tax_table[
  common_features,
  ,
  drop = FALSE
]

############################

# 12. Convert abundance to relative abundance

############################

sample_sums <- colSums(
  feature_table,
  na.rm = TRUE
)

# Remove samples with zero total abundance

valid_samples <- sample_sums > 0

if (!all(valid_samples)) {
  
  warning(
    paste(
      "Removing",
      sum(!valid_samples),
      "samples with zero total abundance."
    )
  )
  
  feature_table <- feature_table[
    ,
    valid_samples,
    drop = FALSE
  ]
  
  sample_sums <- sample_sums[
    valid_samples
  ]
}

relative_abundance <- sweep(
  feature_table,
  2,
  sample_sums,
  "/"
)

feature_table <- relative_abundance

############################

# 13. Function for LEfSe analysis

############################

run_lefse <- function(
    group1,
    group2,
    output_prefix
) {
  
  cat("\n========================================\n")
  cat(
    "Running LEfSe:",
    group1,
    "vs",
    group2,
    "\n"
  )
  cat("========================================\n")
  
  # ----------------------------------------------------------
  
  # Select samples
  
  # ----------------------------------------------------------
  
  sample_table <- meta_study %>%
    filter(
      Group2 %in% c(group1, group2)
    )
  
  if (nrow(sample_table) == 0) {
    
    ```
    warning(
      paste(
        "No samples found for:",
        group1,
        "vs",
        group2
      )
    )
    
    return(NULL)
    ```
    
  }
  
  # ----------------------------------------------------------
  
  # Set SampleID as row names
  
  # ----------------------------------------------------------
  
  rownames(sample_table) <- sample_table$SampleID
  
  # ----------------------------------------------------------
  
  # Match abundance table to metadata
  
  # ----------------------------------------------------------
  
  common_samples <- intersect(
    rownames(sample_table),
    colnames(feature_table)
  )
  
  if (length(common_samples) == 0) {
    
    ```
    warning(
      paste(
        "No matched samples found for:",
        group1,
        "vs",
        group2
      )
    )
    
    return(NULL)
    ```
    
  }
  
  sample_table <- sample_table[
    common_samples,
    ,
    drop = FALSE
  ]
  
  feature_table_sub <- feature_table[
    ,
    common_samples,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  
  # Set group factor
  
  # ----------------------------------------------------------
  
  sample_table$Group2 <- factor(
    sample_table$Group2,
    levels = c(group1, group2)
  )
  
  # ----------------------------------------------------------
  
  # Construct microtable object
  
  # ----------------------------------------------------------
  
  dataset <- microtable$new(
    sample_table = sample_table,
    otu_table = feature_table_sub,
    tax_table = tax_table
  )
  
  # ----------------------------------------------------------
  
  # Run LEfSe
  
  # ----------------------------------------------------------
  
  lefse <- trans_diff$new(
    dataset = dataset,
    method = "lefse",
    group = "Group2",
    alpha = lefse_alpha,
    lefse_subgroup = NULL,
    taxa_level = "Genus",
    lefse_min_subsam = lefse_min_subsam
  )
  
  # ----------------------------------------------------------
  
  # Extract complete results
  
  # ----------------------------------------------------------
  
  all_results <- lefse$res_diff
  
  write.csv(
    all_results,
    file = file.path(
      result_dir,
      paste0(
        output_prefix,
        "_all.csv"
      )
    ),
    row.names = TRUE
  )
  
  # ----------------------------------------------------------
  
  # Filter significant genera
  
  # ----------------------------------------------------------
  
  significant_results <- all_results %>%
    filter(
      LDA > lda_cutoff
    )
  
  # ----------------------------------------------------------
  
  # Add direction
  
  # ----------------------------------------------------------
  
  if (nrow(significant_results) > 0) {
    
    ```
    # The exact column containing the enriched group can vary
    # between microeco versions. Therefore, the original
    # statistical results are retained without modifying
    # the package-generated group labels.
    ```
    
  }
  
  # ----------------------------------------------------------
  
  # Save significant results
  
  # ----------------------------------------------------------
  
  write.csv(
    significant_results,
    file = file.path(
      result_dir,
      paste0(
        output_prefix,
        "_LDA",
        lda_cutoff,
        ".csv"
      )
    ),
    row.names = TRUE
  )
  
  # ----------------------------------------------------------
  
  # Print summary
  
  # ----------------------------------------------------------
  
  cat(
    "Number of samples:",
    nrow(sample_table),
    "\n"
  )
  
  cat(
    "Number of significant taxa:",
    nrow(significant_results),
    "\n"
  )
  
  return(
    list(
      dataset = dataset,
      lefse = lefse,
      all_results = all_results,
      significant_results = significant_results
    )
  )
}

############################

# 14. Run LEfSe comparisons

############################

# ------------------------------------------------------------

# ADLS vs HC

# ------------------------------------------------------------

result_ADLS_HC <- run_lefse(
  group1 = "ADLS",
  group2 = "HC",
  output_prefix = "lefse_genus_ADLS_HC"
)

# ------------------------------------------------------------

# ADNLS vs HC

# ------------------------------------------------------------

result_ADNLS_HC <- run_lefse(
  group1 = "ADNLS",
  group2 = "HC",
  output_prefix = "lefse_genus_ADNLS_HC"
)

# ------------------------------------------------------------

# ADLS vs ADNLS

# ------------------------------------------------------------

result_ADLS_ADNLS <- run_lefse(
  group1 = "ADLS",
  group2 = "ADNLS",
  output_prefix = "lefse_genus_ADLS_ADNLS"
)

############################

# 15. Save session information

############################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    result_dir,
    "sessionInfo_LEfSe.txt"
  )
)

############################

# 16. Finished

############################

cat("\n========================================\n")
cat("LEfSe analysis completed successfully.\n")
cat("Results saved to:\n")
cat(result_dir, "\n")
cat("========================================\n")


