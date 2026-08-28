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
  "patchwork",
  "randomForest",
  "caret",
  "pROC"
)

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
  suppressPackageStartupMessages(
    library(
      pkg,
      character.only = TRUE
    )
  )
}

# MMUPHin

#

# MMUPHin is installed from GitHub because the analysis

# requires the development version used in the original study.

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

# 2. Project directories

############################

# The script is designed to be run from the root directory

# of the GitHub repository.

#

# Recommended structure:

#

# project/

# ├── R/

# ├── data/

# ├── results/

# └── figures/

#

# In RStudio, open the repository as an RStudio Project and

# run this script from the project root.

project_dir <- "."

############################

# 3. Define input directories

############################

data_dir <- file.path(
  project_dir,
  "data",
  "all"
)

############################

# 4. Define output directories

############################

batch_dir <- file.path(
  project_dir,
  "results",
  "batch_correction"
)

meta_dir <- file.path(
  project_dir,
  "results",
  "meta_analysis"
)

ml_dir <- file.path(
  project_dir,
  "results",
  "machine_learning"
)

############################

# 5. Create output directories

############################

dir.create(
  batch_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  meta_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ml_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################

# 6. Input files

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

# 7. Check input files

############################

input_files <- c(
  metadata_file,
  abundance_file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste(
      "The following input files were not found:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

############################

# 8. Analysis parameters

############################

# ------------------------------------------------------------

# Disease groups

# ------------------------------------------------------------

case_group <- "ADLS"

control_group <- "HC"

# ------------------------------------------------------------

# Training cohorts

# ------------------------------------------------------------

train_projects <- c(
  "S2",
  "S17",
  "S19",
  "S47",
  "S71",
  "S75"
)

# ------------------------------------------------------------

# Independent validation cohort

# ------------------------------------------------------------

validation_project <- "S74"

# ------------------------------------------------------------

# Cross-validation

# ------------------------------------------------------------

K <- 5

# ------------------------------------------------------------

# Random seed

# ------------------------------------------------------------

random_seed <- 12

# ------------------------------------------------------------

# RFE

#

# By default, all candidate feature sizes are evaluated.

# For datasets with a very large number of features, this

# can be computationally intensive.

# ------------------------------------------------------------

rfe_use_all_sizes <- TRUE

# ------------------------------------------------------------

# Random Forest tuning

# ------------------------------------------------------------

rf_tune_length <- 3

############################

# 9. Load metadata

############################

meta.all <- read.csv(
  file = metadata_file,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

############################

# 10. Check required metadata columns

############################

required_metadata_columns <- c(
  "Project",
  "Project_rawID",
  "Group2"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(meta.all)
)

if (length(missing_metadata_columns) > 0) {
  
  stop(
    paste(
      "The metadata table is missing the following columns:",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    )
  )
}

############################

# 11. Order cohorts

############################

project_order <- c(
  train_projects,
  validation_project
)

project_order <- c(
  project_order,
  setdiff(
    unique(
      as.character(
        meta.all$Project
      )
    ),
    project_order
  )
)

meta.all$Project <- factor(
  meta.all$Project,
  levels = project_order
)

meta.all <- meta.all[
  order(meta.all$Project),
  ,
  drop = FALSE
]

############################

# 12. Load genus-level abundance data

############################

feat.abu <- read.table(
  file = abundance_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

############################

# 13. Aggregate duplicate genera

############################

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

feat.abu <- feat.abu[
  ,
  setdiff(
    colnames(feat.abu),
    "genus"
  ),
  drop = FALSE
]

############################

# 14. Match metadata and abundance table

############################

common_samples <- intersect(
  rownames(meta.all),
  colnames(feat.abu)
)

if (length(common_samples) == 0) {
  
  stop(
    "No common samples were found between metadata and abundance data."
  )
}

meta.all <- meta.all[
  common_samples,
  ,
  drop = FALSE
]

feat.abu <- feat.abu[
  ,
  rownames(meta.all),
  drop = FALSE
]

############################

# 15. Replace missing values

############################

feat.abu[is.na(feat.abu)] <- 0

############################

# 16. Remove samples with zero total abundance

############################

zero.samples <- colnames(feat.abu)[
  colSums(feat.abu) == 0
]

if (length(zero.samples) > 0) {
  
  message(
    "Removing ",
    length(zero.samples),
    " samples with zero total abundance."
  )
  
  feat.abu <- feat.abu[
    ,
    !colnames(feat.abu) %in% zero.samples,
    drop = FALSE
  ]
  
  meta.all <- meta.all[
    !rownames(meta.all) %in% zero.samples,
    ,
    drop = FALSE
  ]
}

############################

# 17. Clean factor levels

############################

meta.all$Project <- droplevels(
  meta.all$Project
)

############################

# 18. Reorder abundance matrix

############################

feat.abu <- feat.abu[
  ,
  rownames(meta.all),
  drop = FALSE
]

############################

# 19. Check data consistency

############################

stopifnot(
  identical(
    colnames(feat.abu),
    rownames(meta.all)
  )
)

stopifnot(
  all(
    colSums(feat.abu) > 0
  )
)

############################

# 20. MMUPHin batch correction

############################

message("")
message("==========================================")
message("MMUPHin batch-effect correction")
message("==========================================")

# Convert percentage abundance to proportions

feat.abu_adj <- feat.abu / 100

fit_adjust_batch <- adjust_batch(
  feature_abd = feat.abu_adj,
  batch = "Project",
  data = meta.all,
  control = list(
    verbose = FALSE
  )
)

npc_abd_adj <- fit_adjust_batch$feature_abd_adj

# Convert proportions back to percentage

npc_abd_adj <- npc_abd_adj * 100

############################

# 21. Save batch-corrected abundance

############################

batch_corrected_file <- file.path(
  batch_dir,
  "genus_MMUPHin_batch_corrected.csv"
)

write.csv(
  npc_abd_adj,
  file = batch_corrected_file,
  row.names = TRUE
)

############################

# 22. ADLS vs HC meta-analysis

############################

message("")
message("==========================================")
message("ADLS vs HC meta-analysis")
message("==========================================")

meta.all2 <- meta.all %>%
  filter(
    Group2 %in% c(
      case_group,
      control_group
    )
  ) %>%
  filter(
    Project_rawID %in% train_projects
  )

############################

# 23. Match corrected abundance data

############################

selected_samples <- intersect(
  rownames(meta.all2),
  colnames(npc_abd_adj)
)

if (length(selected_samples) == 0) {
  
  stop(
    "No samples were found for the ADLS vs HC meta-analysis."
  )
}

npc_abd_adj2 <- npc_abd_adj[
  ,
  selected_samples,
  drop = FALSE
]

meta.all2 <- meta.all2[
  selected_samples,
  ,
  drop = FALSE
]

############################

# 24. Set factor levels

############################

meta.all2$Project <- factor(
  meta.all2$Project,
  levels = unique(
    meta.all2$Project
  )
)

meta.all2$Group2 <- factor(
  meta.all2$Group2,
  levels = c(
    case_group,
    control_group
  )
)

############################

# 25. Meta-analytical differential abundance

############################

fit_lm_meta <- lm_meta(
  feature_abd = npc_abd_adj2,
  batch = "Project",
  exposure = "Group2",
  data = meta.all2
)

meta_fits <- fit_lm_meta$meta_fits

############################

# 26. Save meta-analysis results

############################

meta_result_file <- file.path(
  meta_dir,
  "lm_meta_ADLS_HC.csv"
)

write.csv(
  meta_fits,
  file = meta_result_file,
  row.names = TRUE
)

############################

# 27. Save significant meta-analysis results

############################

significant_meta <- meta_fits %>%
  filter(
    qval.fdr < 0.05
  )

write.csv(
  significant_meta,
  file = file.path(
    meta_dir,
    "lm_meta_ADLS_HC_FDR05.csv"
  ),
  row.names = TRUE
)

############################

# 28. Prepare data for machine learning

############################

message("")
message("==========================================")
message("Random Forest classification")
message("==========================================")

############################

# 29. Restrict to ADLS and HC

############################

meta.ml <- meta.all %>%
  filter(
    Group2 %in% c(
      case_group,
      control_group
    )
  )

############################

# 30. Define training and validation cohorts

############################

meta.train <- meta.ml %>%
  filter(
    Project_rawID %in% train_projects
  )

meta.valid <- meta.ml %>%
  filter(
    Project_rawID == validation_project
  )

############################

# 31. Check training and validation samples

############################

cat("\n==============================\n")
cat("Training cohort distribution:\n")
cat("==============================\n")

print(
  table(
    meta.train$Project_rawID,
    meta.train$Group2
  )
)

cat("\n==============================\n")
cat("Independent validation cohort:\n")
cat("==============================\n")

print(
  table(
    meta.valid$Project_rawID,
    meta.valid$Group2
  )
)

############################

# 32. Check sample overlap

############################

overlap_samples <- intersect(
  rownames(meta.train),
  rownames(meta.valid)
)

if (length(overlap_samples) > 0) {
  
  stop(
    "Training and validation datasets contain overlapping samples."
  )
  
} else {
  
  cat(
    "\nTraining and validation datasets are completely independent.\n"
  )
}

############################

# 33. Extract training and validation abundance

############################

genus_train <- npc_abd_adj[
  ,
  rownames(meta.train),
  drop = FALSE
]

genus_valid <- npc_abd_adj[
  ,
  rownames(meta.valid),
  drop = FALSE
]

############################

# 34. Check sample order

############################

stopifnot(
  identical(
    colnames(genus_train),
    rownames(meta.train)
  )
)

stopifnot(
  identical(
    colnames(genus_valid),
    rownames(meta.valid)
  )
)

############################

# 35. Extract significant features

############################

diff_features <- unique(
  significant_meta$feature
)

diff_features <- diff_features[
  diff_features %in% rownames(npc_abd_adj)
]

cat(
  "\nNumber of candidate differential genera: ",
  length(diff_features),
  "\n",
  sep = ""
)

if (length(diff_features) == 0) {
  
  stop(
    "No significant genera were found for machine-learning analysis."
  )
}

############################

# 36. Extract differential genera

############################

abun_diff_train <- genus_train[
  diff_features,
  ,
  drop = FALSE
]

abun_diff_valid <- genus_valid[
  diff_features,
  ,
  drop = FALSE
]

############################

# 37. Ensure identical feature sets

############################

common_features <- intersect(
  rownames(abun_diff_train),
  rownames(abun_diff_valid)
)

abun_diff_train <- abun_diff_train[
  common_features,
  ,
  drop = FALSE
]

abun_diff_valid <- abun_diff_valid[
  common_features,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(abun_diff_train),
    rownames(abun_diff_valid)
  )
)

############################

# 38. Construct metadata

############################

meta_train <- data.frame(
  SampleID = rownames(meta.train),
  Group = meta.train$Group2,
  Project = meta.train$Project_rawID,
  row.names = rownames(meta.train),
  check.names = FALSE
)

meta_valid <- data.frame(
  SampleID = rownames(meta.valid),
  Group = meta.valid$Group2,
  Project = meta.valid$Project_rawID,
  row.names = rownames(meta.valid),
  check.names = FALSE
)

############################

# 39. Transpose abundance matrix

############################

abun_diff_train_transposed <- t(
  abun_diff_train
)

abun_diff_valid_transposed <- t(
  abun_diff_valid
)

############################

# 40. Combine metadata and abundance

############################

train_data <- cbind(
  meta_train,
  abun_diff_train_transposed
)

test_data <- cbind(
  meta_valid,
  abun_diff_valid_transposed
)

############################

# 41. Set class labels

############################

train_data$Group <- factor(
  train_data$Group,
  levels = c(
    case_group,
    control_group
  )
)

test_data$Group <- factor(
  test_data$Group,
  levels = c(
    case_group,
    control_group
  )
)

train_data$Project <- as.character(
  train_data$Project
)

test_data$Project <- as.character(
  test_data$Project
)

############################

# 42. Study-stratified 5-fold CV

############################

set.seed(
  random_seed
)

fold_test <- vector(
  mode = "list",
  length = K
)

names(fold_test) <- paste0(
  "Fold",
  seq_len(K)
)

############################

# 43. Generate folds within each cohort

############################

for (study in train_projects) {
  
  cat(
    "\n------------------------------------\n"
  )
  
  cat(
    "Generating stratified folds for cohort: ",
    study,
    "\n",
    sep = ""
  )
  
  cat(
    "------------------------------------\n"
  )
  
  study_idx <- which(
    train_data$Project == study
  )
  
  if (length(study_idx) == 0) {
    
    ```
    warning(
      paste0(
        "No samples found for cohort ",
        study,
        ". Skipping."
      )
    )
    
    next
    ```
    
  }
  
  study_group <- train_data$Group[
    study_idx
  ]
  
  print(
    table(
      study_group
    )
  )
  
  if (length(study_idx) < K) {
    
    ```
    stop(
      paste0(
        "Cohort ",
        study,
        " contains fewer than ",
        K,
        " samples."
      )
    )
    ```
    
  }
  
  class_counts <- table(
    study_group
  )
  
  if (any(class_counts < K)) {
    
    ```
    warning(
      paste0(
        "At least one class in cohort ",
        study,
        " contains fewer than ",
        K,
        " samples. ",
        "Each fold may not contain both classes."
      )
    )
    ```
    
  }
  
  study_folds <- createFolds(
    y = study_group,
    k = K,
    list = TRUE,
    returnTrain = FALSE
  )
  
  for (k in seq_len(K)) {
    
    ```
    local_test_idx <- study_folds[[k]]
    
    global_test_idx <- study_idx[
      local_test_idx
    ]
    
    fold_test[[k]] <- c(
      fold_test[[k]],
      global_test_idx
    )
    ```
    
  }
}

############################

# 44. Sort fold indices

############################

fold_test <- lapply(
  fold_test,
  sort
)

############################

# 45. Generate training indices

############################

all_train_rows <- seq_len(
  nrow(train_data)
)

fold_train <- lapply(
  fold_test,
  function(test_idx) {
    
    ```
    setdiff(
      all_train_rows,
      test_idx
    )
    ```
    
  }
)

names(fold_train) <- names(
  fold_test
)

############################

# 46. Check fold assignment

############################

all_validation_indices <- unlist(
  fold_test
)

if (
  length(all_validation_indices) !=
  nrow(train_data)
) {
  
  stop(
    "The total number of validation samples does not match the training set."
  )
}

if (
  length(unique(all_validation_indices)) !=
  nrow(train_data)
) {
  
  stop(
    "Some samples appear in multiple validation folds or are missing."
  )
}

cat(
  "\nStudy-stratified 5-fold CV check passed.\n"
)

############################

# 47. Save fold assignment

############################

fold_assignment <- data.frame(
  SampleID = rownames(train_data),
  Project = train_data$Project,
  Group = train_data$Group,
  Fold = NA_character_,
  stringsAsFactors = FALSE
)

for (k in seq_len(K)) {
  
  fold_assignment$Fold[
    fold_test[[k]]
  ] <- paste0(
    "Fold",
    k
  )
}

fold_assignment <- fold_assignment %>%
  arrange(
    Project,
    Fold,
    Group
  )

write.csv(
  fold_assignment,
  file = file.path(
    ml_dir,
    "study_stratified_5fold_assignment.csv"
  ),
  row.names = FALSE
)

############################

# 48. Define candidate features

############################

feature_columns <- setdiff(
  colnames(train_data),
  c(
    "SampleID",
    "Group",
    "Project"
  )
)

cat(
  "\nCandidate features for RFE: ",
  length(feature_columns),
  "\n",
  sep = ""
)

############################

# 49. RFE control

############################

rfe_control <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = K,
  index = fold_train,
  verbose = FALSE,
  returnResamp = "final"
)

############################

# 50. Define RFE feature sizes

############################

if (rfe_use_all_sizes) {
  
  rfe_sizes <- seq_len(
    length(feature_columns)
  )
  
} else {
  
  rfe_sizes <- unique(
    c(
      1,
      5,
      10,
      20
    )
  )
  
  rfe_sizes <- rfe_sizes[
    rfe_sizes <= length(feature_columns)
  ]
}

############################

# 51. Run RFE

############################

set.seed(
  random_seed
)

rfe_result <- rfe(
  x = train_data[
    ,
    feature_columns,
    drop = FALSE
  ],
  y = train_data$Group,
  sizes = rfe_sizes,
  rfeControl = rfe_control
)

############################

# 52. Save RFE results

############################

write.csv(
  rfe_result$results,
  file = file.path(
    ml_dir,
    "RFE_results.csv"
  ),
  row.names = FALSE
)

############################

# 53. Extract optimal features

############################

best_features <- rfe_result$optVariables

cat(
  "\n============================================\n"
)

cat(
  "Number of RFE-selected features: ",
  length(best_features),
  "\n",
  sep = ""
)

cat(
  "============================================\n"
)

print(
  best_features
)

############################

# 54. Save selected features

############################

write.csv(
  data.frame(
    Feature = best_features
  ),
  file = file.path(
    ml_dir,
    "RFE_selected_features.csv"
  ),
  row.names = FALSE
)

############################

# 55. Check selected features

############################

missing_train_features <- setdiff(
  best_features,
  colnames(train_data)
)

missing_test_features <- setdiff(
  best_features,
  colnames(test_data)
)

if (length(missing_train_features) > 0) {
  
  stop(
    "Some RFE-selected features are missing from the training data."
  )
}

if (length(missing_test_features) > 0) {
  
  stop(
    paste(
      "Some RFE-selected features are missing from the validation data:",
      paste(
        missing_test_features,
        collapse = ", "
      )
    )
  )
}

############################

# 56. Final Random Forest model

############################

train_control <- trainControl(
  method = "cv",
  number = K,
  index = fold_train,
  indexOut = fold_test,
  classProbs = TRUE,
  savePredictions = "final",
  summaryFunction = twoClassSummary,
  verboseIter = FALSE
)

############################

# 57. Train Random Forest

############################

set.seed(
  random_seed
)

rf_model <- train(
  x = train_data[
    ,
    best_features,
    drop = FALSE
  ],
  y = train_data$Group,
  method = "rf",
  trControl = train_control,
  tuneLength = rf_tune_length,
  metric = "ROC",
  importance = TRUE
)

############################

# 58. Save Random Forest model

############################

saveRDS(
  rf_model,
  file = file.path(
    ml_dir,
    "random_forest_model.rds"
  )
)

############################

# 59. Save variable importance

############################

variable_importance <- varImp(
  rf_model
)$importance

variable_importance$Feature <- rownames(
  variable_importance
)

variable_importance <- variable_importance %>%
  arrange(
    desc(
      Overall
    )
  )

write.csv(
  variable_importance,
  file = file.path(
    ml_dir,
    "random_forest_variable_importance.csv"
  ),
  row.names = FALSE
)

############################

# 60. Internal CV predictions

############################

cv_predictions <- rf_model$pred

if ("mtry" %in% colnames(cv_predictions)) {
  
  cv_predictions_best <- cv_predictions[
    cv_predictions$mtry ==
      rf_model$bestTune$mtry,
    ,
    drop = FALSE
  ]
  
} else {
  
  cv_predictions_best <- cv_predictions
}

write.csv(
  cv_predictions_best,
  file = file.path(
    ml_dir,
    "random_forest_5fold_CV_predictions.csv"
  ),
  row.names = FALSE
)

############################

# 61. Internal CV ROC

############################

cv_roc <- roc(
  response = cv_predictions_best$obs,
  predictor = cv_predictions_best[[case_group]],
  levels = c(
    control_group,
    case_group
  ),
  direction = "<",
  quiet = TRUE
)

cv_auc <- auc(
  cv_roc
)

cat(
  "\nStudy-stratified 5-fold CV AUC: ",
  as.numeric(cv_auc),
  "\n",
  sep = ""
)

############################

# 62. Determine optimal threshold

############################

best_threshold_info <- coords(
  cv_roc,
  x = "best",
  best.method = "youden",
  ret = c(
    "threshold",
    "sensitivity",
    "specificity"
  ),
  transpose = FALSE
)

best_threshold <- as.numeric(
  best_threshold_info["threshold"]
)

cat(
  "\nSelected threshold from training CV: ",
  best_threshold,
  "\n",
  sep = ""
)

############################

# 63. Independent validation

############################

test_probs <- predict(
  rf_model,
  newdata = test_data[
    ,
    best_features,
    drop = FALSE
  ],
  type = "prob"
)

############################

# 64. Independent validation classification

############################

test_predictions <- ifelse(
  test_probs[
    ,
    case_group
  ] >= best_threshold,
  case_group,
  control_group
)

test_predictions <- factor(
  test_predictions,
  levels = c(
    case_group,
    control_group
  )
)

############################

# 65. Confusion matrix

############################

conf_matrix <- confusionMatrix(
  data = test_predictions,
  reference = test_data$Group,
  positive = case_group
)

print(
  conf_matrix
)

############################

# 66. Save confusion matrix

############################

conf_matrix_df <- as.data.frame(
  conf_matrix$table
)

colnames(conf_matrix_df) <- c(
  "Prediction",
  "Reference",
  "Frequency"
)

write.csv(
  conf_matrix_df,
  file = file.path(
    ml_dir,
    paste0(
      validation_project,
      "_independent_validation_confusion_matrix.csv"
    )
  ),
  row.names = FALSE
)

############################

# 67. Independent validation metrics

############################

accuracy <- unname(
  conf_matrix$overall[
    "Accuracy"
  ]
)

precision <- unname(
  conf_matrix$byClass[
    "Pos Pred Value"
  ]
)

recall <- unname(
  conf_matrix$byClass[
    "Sensitivity"
  ]
)

specificity <- unname(
  conf_matrix$byClass[
    "Specificity"
  ]
)

balanced_accuracy <- unname(
  conf_matrix$byClass[
    "Balanced Accuracy"
  ]
)

f1_score <- ifelse(
  is.na(precision) |
    is.na(recall) |
    precision + recall == 0,
  NA,
  2 * precision * recall /
    (precision + recall)
)

############################

# 68. Independent validation ROC

############################

validation_roc <- roc(
  response = test_data$Group,
  predictor = test_probs[
    ,
    case_group
  ],
  levels = c(
    control_group,
    case_group
  ),
  direction = "<",
  quiet = TRUE
)

############################

# 69. Independent validation AUC

############################

validation_auc <- auc(
  validation_roc
)

validation_auc_ci <- ci.auc(
  validation_roc
)

############################

# 70. Compile performance metrics

############################

performance_df <- data.frame(
  Validation_Cohort = validation_project,
  Threshold = best_threshold,
  Accuracy = accuracy,
  Sensitivity = recall,
  Specificity = specificity,
  Precision = precision,
  F1_Score = f1_score,
  Balanced_Accuracy = balanced_accuracy,
  AUC = as.numeric(validation_auc),
  AUC_CI_Lower = as.numeric(
    validation_auc_ci[1]
  ),
  AUC_CI_Upper = as.numeric(
    validation_auc_ci[3]
  ),
  stringsAsFactors = FALSE
)

############################

# 71. Save validation performance

############################

write.csv(
  performance_df,
  file = file.path(
    ml_dir,
    paste0(
      validation_project,
      "_independent_validation_performance.csv"
    )
  ),
  row.names = FALSE
)

############################

# 72. Save validation predictions

############################

validation_predictions <- data.frame(
  SampleID = rownames(test_data),
  Project = test_data$Project,
  Observed_Group = test_data$Group,
  Predicted_Group = test_predictions,
  ADLS_Probability = test_probs[
    ,
    case_group
  ],
  HC_Probability = test_probs[
    ,
    control_group
  ],
  stringsAsFactors = FALSE
)

write.csv(
  validation_predictions,
  file = file.path(
    ml_dir,
    paste0(
      validation_project,
      "_independent_validation_predictions.csv"
    )
  ),
  row.names = FALSE
)

############################

# 73. Save threshold information

############################

threshold_df <- data.frame(
  Threshold = best_threshold,
  CV_Sensitivity = as.numeric(
    best_threshold_info[
      "sensitivity"
    ]
  ),
  CV_Specificity = as.numeric(
    best_threshold_info[
      "specificity"
    ]
  ),
  CV_AUC = as.numeric(
    cv_auc
  )
)

write.csv(
  threshold_df,
  file = file.path(
    ml_dir,
    "training_CV_threshold.csv"
  ),
  row.names = FALSE
)

############################

# 74. Save session information

############################

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    ml_dir,
    "sessionInfo_MMUPHin_RF.txt"
  )
)

############################

# 75. Final summary

############################

message("")
message("==========================================")
message("Analysis completed successfully.")
message("==========================================")

message(
  "Batch-corrected abundance: ",
  batch_corrected_file
)

message(
  "Meta-analysis results: ",
  meta_result_file
)

message(
  "Number of candidate genera: ",
  length(diff_features)
)

message(
  "Number of RFE-selected genera: ",
  length(best_features)
)

message(
  "Training CV AUC: ",
  round(
    as.numeric(cv_auc),
    4
  )
)

message(
  "Independent validation AUC: ",
  round(
    as.numeric(validation_auc),
    4
  )
)

message(
  "Independent validation cohort: ",
  validation_project
)

message(
  "All machine-learning results are saved in: ",
  ml_dir
)

