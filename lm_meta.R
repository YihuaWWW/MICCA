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

# ============================================================

############################

# 1. Install and load packages

############################

cran_packages <- c(
  "magrittr",
  "dplyr",
  "ggplot2",
  "vegan",
  "ggpubr",
  "ggsci",
  "patchwork"
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

# Install MMUPHin from GitHub if necessary

if (!requireNamespace("MMUPHin", quietly = TRUE)) {
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  remotes::install_github(
    "biobakery/MMUPHin"
  )
}

suppressPackageStartupMessages(
  library(MMUPHin)
)

############################

# 2. User-defined parameters

############################

# ------------------------------------------------------------

# Project directory

#

# Recommended GitHub structure:

#

# repository/

# ├── data/

# │   └── all/

# │       ├── group.csv

# │       └── tax_L6.txt

# ├── R/

# │   └── mmuphin_meta_analysis.R

# └── results/

# ------------------------------------------------------------

project_dir <- "."

data_dir <- file.path(
  project_dir,
  "data",
  "all"
)

result_dir <- file.path(
  project_dir,
  "results"
)

figure_dir <- file.path(
  project_dir,
  "figures"
)

############################

# 3. Create output directories

############################

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(figure_dir, "lm_meta"),
  recursive = TRUE,
  showWarnings = FALSE
)

############################

# 4. Input files

############################

metadata_file <- file.path(
  data_dir,
  "group.csv"
)

abundance_file <- file.path(
  data_dir,
  "tax_L6.txt"
)

############################

# 5. Read metadata

############################

meta.all <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

############################

# 6. Check required metadata columns

############################

required_columns <- c(
  "Project",
  "Project_rawID",
  "Group2"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta.all)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "The following required metadata columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

############################

# 7. Read genus-level abundance data

############################

feat.abu <- read.table(
  abundance_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################

# 8. Remove problematic samples

############################

# These samples were excluded in the original analysis.

# Modify this vector if different samples need to be removed.

samples_to_remove <- c(
  "ERR6650487",
  "ERR6650884",
  "ERR6650887"
)

samples_to_remove <- intersect(
  samples_to_remove,
  colnames(feat.abu)
)

if (length(samples_to_remove) > 0) {
  
  feat.abu <- feat.abu %>%
    select(
      -all_of(samples_to_remove)
    )
}

############################

# 9. Aggregate genus-level abundance

############################

# If multiple rows correspond to the same genus,

# their abundances are summed.

if (!"genus" %in% colnames(feat.abu)) {
  
  stop(
    "The abundance table must contain a column named 'genus'."
  )
}

feat.abu <- aggregate(
  . ~ genus,
  data = feat.abu,
  sum
)

rownames(feat.abu) <- feat.abu$genus

feat.abu <- feat.abu[, -1, drop = FALSE]

############################

# 10. Match metadata and abundance data

############################

common_samples <- intersect(
  rownames(meta.all),
  colnames(feat.abu)
)

if (length(common_samples) == 0) {
  
  stop(
    "No overlapping samples were found between metadata and abundance data."
  )
}

meta.all <- meta.all[
  common_samples,
  ,
  drop = FALSE
]

feat.abu <- feat.abu[
  ,
  common_samples,
  drop = FALSE
]

############################

# 11. Remove missing values

############################

feat.abu[is.na(feat.abu)] <- 0

############################

# 12. Convert percentages to proportions

############################

# The original abundance table was expressed as percentage.

# MMUPHin requires compositional abundance data.

feat.abu <- feat.abu / 100

############################

# 13. Define batch factor

############################

# The original analysis used Project as the batch variable.

meta.all$Project <- factor(
  meta.all$Project
)

############################

# 14. Batch-effect correction using MMUPHin

############################

fit_adjust_batch <- adjust_batch(
  feature_abd = feat.abu,
  batch = "Project",
  data = meta.all,
  control = list(
    verbose = FALSE
  )
)

############################

# 15. Extract adjusted abundance matrix

############################

npc_abd_adj <- fit_adjust_batch$feature_abd_adj

# Convert proportions back to percentages

npc_abd_adj <- npc_abd_adj * 100

############################

# 16. Function for meta-analytical analysis

############################

run_mmuphin_analysis <- function(
    group1,
    group2,
    project_ids,
    project_levels,
    output_prefix,
    make_plot = TRUE
) {
  
  cat("\n========================================\n")
  cat("Running:", group1, "vs", group2, "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  
  # Select groups
  
  # ----------------------------------------------------------
  
  meta_sub <- meta.all %>%
    filter(
      Group2 %in% c(group1, group2)
    )
  
  # ----------------------------------------------------------
  
  # Select cohorts
  
  # ----------------------------------------------------------
  
  meta_sub <- meta_sub %>%
    filter(
      Project_rawID %in% project_ids
    )
  
  # ----------------------------------------------------------
  
  # Check sample number
  
  # ----------------------------------------------------------
  
  if (nrow(meta_sub) == 0) {
    
    ```
    warning(
      paste(
        "No samples available for",
        group1,
        "vs",
        group2
      )
    )
    
    return(NULL)
    ```
    
  }
  
  # ----------------------------------------------------------
  
  # Match abundance matrix
  
  # ----------------------------------------------------------
  
  sample_ids <- intersect(
    rownames(meta_sub),
    colnames(npc_abd_adj)
  )
  
  if (length(sample_ids) == 0) {
    
    ```
    warning(
      paste(
        "No matched abundance data for",
        group1,
        "vs",
        group2
      )
    )
    
    return(NULL)
    ```
    
  }
  
  meta_sub <- meta_sub[
    sample_ids,
    ,
    drop = FALSE
  ]
  
  abundance_sub <- npc_abd_adj[
    ,
    sample_ids,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  
  # Set factor levels
  
  # ----------------------------------------------------------
  
  meta_sub$Project <- factor(
    meta_sub$Project,
    levels = project_levels
  )
  
  meta_sub$Group2 <- factor(
    meta_sub$Group2,
    levels = c(group1, group2)
  )
  
  # Remove unused batch levels
  
  meta_sub$Project <- droplevels(
    meta_sub$Project
  )
  
  # ----------------------------------------------------------
  
  # Meta-analytical differential abundance analysis
  
  # ----------------------------------------------------------
  
  fit_lm_meta <- lm_meta(
    feature_abd = abundance_sub,
    batch = "Project",
    exposure = "Group2",
    data = meta_sub,
    control = list(
      verbose = FALSE
    )
  )
  
  meta_fits <- fit_lm_meta$meta_fits
  
  # ----------------------------------------------------------
  
  # Save complete results
  
  # ----------------------------------------------------------
  
  write.csv(
    meta_fits,
    file.path(
      result_dir,
      paste0(
        output_prefix,
        "_all.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  
  # Significant features
  
  # ----------------------------------------------------------
  
  signif_res <- meta_fits %>%
    filter(
      qval.fdr < 0.05
    ) %>%
    arrange(
      coef
    ) %>%
    mutate(
      direction = ifelse(
        coef > 0,
        group2,
        group1
      )
    )
  
  write.csv(
    signif_res,
    file.path(
      result_dir,
      paste0(
        output_prefix,
        "_FDR0.05.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  
  # Print summary
  
  # ----------------------------------------------------------
  
  cat(
    "Number of samples:",
    nrow(meta_sub),
    "\n"
  )
  
  cat(
    "Number of significant genera:",
    nrow(signif_res),
    "\n"
  )
  
  ############################
  
  # 17. Visualization
  
  ############################
  
  if (make_plot && nrow(signif_res) > 0) {
    
    ```
    # --------------------------------------------------------
    # Bar plot
    # --------------------------------------------------------
    
    bar_data <- signif_res %>%
      mutate(
        feature = factor(
          feature,
          levels = feature
        ),
        direction = ifelse(
          coef > 0,
          group2,
          group1
        )
      )
    
    
    pdf(
      file.path(
        figure_dir,
        "lm_meta",
        paste0(
          "barplot_FDR0.05_",
          output_prefix,
          ".pdf"
        ),
      ),
      width = 8,
      height = 6
    )
    
    
    print(
      ggplot(
        bar_data,
        aes(
          y = coef,
          x = feature,
          fill = direction
        )
      ) +
        geom_bar(
          stat = "identity"
        ) +
        coord_flip() +
        labs(
          x = NULL,
          y = "Meta-analysis coefficient",
          fill = "Direction"
        ) +
        theme_minimal()
    )
    
    
    dev.off()
    
    
    # --------------------------------------------------------
    # Lollipop plot
    # --------------------------------------------------------
    
    pdf(
      file.path(
        figure_dir,
        "lm_meta",
        paste0(
          "lollipop_FDR0.05_",
          output_prefix,
          ".pdf"
        )
      ),
      width = 8,
      height = 6
    )
    
    
    print(
      ggplot(
        bar_data,
        aes(
          x = feature,
          y = coef,
          color = direction
        )
      ) +
        
        geom_segment(
          aes(
            x = feature,
            xend = feature,
            y = 0,
            yend = coef
          ),
          color = "grey"
        ) +
        
        geom_hline(
          yintercept = 0,
          linewidth = 0.1,
          linetype = 2
        ) +
        
        geom_point(
          size = 4
        ) +
        
        labs(
          x = NULL,
          y = "Meta-analysis coefficient",
          color = "Direction"
        ) +
        
        coord_flip() +
        
        theme_light() +
        
        theme(
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank()
        )
    )
    
    
    dev.off()
    ```
    
  }
  
  # Return results
  
  return(
    list(
      fit = fit_lm_meta,
      all_results = meta_fits,
      significant_results = signif_res
    )
  )
}

############################

# 18. Run all comparisons

############################

# ------------------------------------------------------------

# ADLS vs HC

# ------------------------------------------------------------

result_ADLS_HC <- run_mmuphin_analysis(
  group1 = "ADLS",
  group2 = "HC",
  
  project_ids = c(
    "S2",
    "S17",
    "S19",
    "S47",
    "S71",
    "S74",
    "S75"
  ),
  
  project_levels = c(
    "S1",
    "S2",
    "S3",
    "S7",
    "S8",
    "S9",
    "S10"
  ),
  
  output_prefix = "lm_meta_ADLS_HC"
)

# ------------------------------------------------------------

# ADNLS vs HC

# ------------------------------------------------------------

result_ADNLS_HC <- run_mmuphin_analysis(
  group1 = "ADNLS",
  group2 = "HC",
  
  project_ids = c(
    "S2",
    "S17",
    "S71",
    "S74"
  ),
  
  project_levels = c(
    "S1",
    "S2",
    "S8",
    "S9"
  ),
  
  output_prefix = "lm_meta_ADNLS_HC"
)

# ------------------------------------------------------------

# ADLS vs ADNLS

# ------------------------------------------------------------

result_ADLS_ADNLS <- run_mmuphin_analysis(
  group1 = "ADLS",
  group2 = "ADNLS",
  
  project_ids = c(
    "S2",
    "S17",
    "S24",
    "S33",
    "S71",
    "S74"
  ),
  
  project_levels = c(
    "S1",
    "S2",
    "S4",
    "S6",
    "S8",
    "S9"
  ),
  
  output_prefix = "lm_meta_ADLS_ADNLS"
)

############################

# 19. Save session information

############################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    result_dir,
    "sessionInfo.txt"
  )
)

############################

# 20. Finished

############################

cat("\n========================================\n")
cat("MMUPHin analysis completed.\n")
cat("Results:", result_dir, "\n")
cat("Figures:", figure_dir, "\n")
cat("========================================\n")


