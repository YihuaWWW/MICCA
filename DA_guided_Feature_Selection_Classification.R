```r
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

project_dir <- normalizePath(
  "..",
  winslash = "/",
  mustWork = FALSE
)

# If running this script directly from the project root,
# replace the previous line with:
#
# project_dir <- "."
#
# Alternatively, users can specify an absolute path:
#
# project_dir <- "/path/to/your/project"


# ============================================================
# 1. Install and load required packages
# ============================================================

required_packages <- c(
  "randomForest",
  "caret",
  "pROC",
  "dplyr"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
  library(
    pkg,
    character.only = TRUE
  )
}


# ============================================================
# 2. Analysis parameters
# ============================================================

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
# Disease groups
# ------------------------------------------------------------

case_group <- "ADLS"
control_group <- "HC"

group_levels <- c(
  case_group,
  control_group
)


# ------------------------------------------------------------
# Number of CV folds
# ------------------------------------------------------------

K <- 5


# ------------------------------------------------------------
# Random seeds
# ------------------------------------------------------------

seed_cv <- 42
seed_model <- 12


# ============================================================
# 3. Define input and output paths
# ============================================================

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

metadata_file <- file.path(
  project_dir,
  "data",
  "all",
  "group.csv"
)


genus_file <- file.path(
  project_dir,
  "data",
  "all",
  "revised_results",
  "genus_NPC_ra_abd_adj.csv"
)


# ------------------------------------------------------------
# Differential genus result directories
# ------------------------------------------------------------

diff_dir <- function(project) {
  
  file.path(
    project_dir,
    "data",
    project,
    "results",
    "diff"
  )
}


# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------

output_dir <- file.path(
  project_dir,
  "data",
  "all",
  "revised_results"
)


if (!dir.exists(output_dir)) {
  
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ============================================================
# 4. Identify robust differential genera
# ============================================================
#
# For each training cohort, three differential-abundance
# approaches are considered:
#
#   1. Permutation test
#   2. Wilcoxon rank-sum test
#   3. LEfSe
#
# A genus is considered robust within a cohort if it is
# identified by all available methods.
#
# The final feature set is obtained by taking the union of
# cohort-specific robust genera.
#
# ============================================================


# ------------------------------------------------------------
# Function to read differential genus results
# ------------------------------------------------------------

read_diff_results <- function(
    project,
    comparison = "ADLS&HC"
) {
  
  current_dir <- diff_dir(project)
  
  
  permutation_file <- file.path(
    current_dir,
    "perm_genus_LS&HC.csv"
  )
  
  
  wilcoxon_file <- file.path(
    current_dir,
    "wilcox_genus_LS&HC.csv"
  )
  
  
  lefse_file <- file.path(
    current_dir,
    paste0(
      "lefse_genus_",
      comparison,
      ".csv"
    )
  )
  
  
  result <- list()
  
  
  # ----------------------------------------------------------
  # Permutation test
  # ----------------------------------------------------------
  
  if (file.exists(permutation_file)) {
    
    result$permutation <- read.csv(
      permutation_file,
      row.names = 1,
      check.names = FALSE
    )
    
  } else {
    
    warning(
      "Permutation result not found: ",
      permutation_file
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Wilcoxon test
  # ----------------------------------------------------------
  
  if (file.exists(wilcoxon_file)) {
    
    result$wilcoxon <- read.csv(
      wilcoxon_file,
      row.names = 1,
      check.names = FALSE
    )
    
  } else {
    
    warning(
      "Wilcoxon result not found: ",
      wilcoxon_file
    )
    
  }
  
  
  # ----------------------------------------------------------
  # LEfSe
  # ----------------------------------------------------------
  
  if (file.exists(lefse_file)) {
    
    result$lefse <- read.csv(
      lefse_file,
      row.names = 1,
      check.names = FALSE
    )
    
  } else {
    
    warning(
      "LEfSe result not found: ",
      lefse_file
    )
    
  }
  
  
  return(result)
}


# ------------------------------------------------------------
# Find robust genera for each training cohort
# ------------------------------------------------------------

cohort_features <- list()


for (project in train_projects) {
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "Processing cohort:",
    project,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  diff_results <- read_diff_results(
    project = project,
    comparison = "ADLS&HC"
  )
  
  
  genus_vectors <- list()
  
  
  if (!is.null(diff_results$permutation)) {
    
    genus_vectors$permutation <-
      diff_results$permutation$genus
    
  }
  
  
  if (!is.null(diff_results$wilcoxon)) {
    
    genus_vectors$wilcoxon <-
      diff_results$wilcoxon$genus
    
  }
  
  
  if (!is.null(diff_results$lefse)) {
    
    genus_vectors$lefse <-
      diff_results$lefse$genus
    
  }
  
  
  # ----------------------------------------------------------
  # If all three methods are available, take their
  # intersection.
  # ----------------------------------------------------------
  
  if (length(genus_vectors) >= 2) {
    
    robust_features <- Reduce(
      intersect,
      genus_vectors
    )
    
  } else if (length(genus_vectors) == 1) {
    
    robust_features <- genus_vectors[[1]]
    
  } else {
    
    robust_features <- character(0)
    
  }
  
  
  cohort_features[[project]] <- robust_features
  
  
  cat(
    "Number of robust genera:",
    length(robust_features),
    "\n"
  )
}


# ------------------------------------------------------------
# Combine cohort-specific robust genera
# ------------------------------------------------------------

fulljoin <- Reduce(
  union,
  cohort_features
)


cat(
  "\n============================================\n"
)

cat(
  "Total candidate genera:",
  length(fulljoin),
  "\n"
)

cat(
  "============================================\n"
)


# Save candidate genera
write.csv(
  data.frame(
    genus = fulljoin
  ),
  file = file.path(
    output_dir,
    "candidate_differential_genera.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. Load metadata
# ============================================================

meta.all <- read.csv(
  file = metadata_file,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# Define Project as factor
# ------------------------------------------------------------

project_order <- c(
  train_projects,
  validation_project
)


meta.all$Project <- factor(
  meta.all$Project,
  levels = c(
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
)


# Sort metadata according to cohort
meta.all <- meta.all[
  order(meta.all$Project),
  ,
  drop = FALSE
]


# ============================================================
# 6. Load MMUPHin batch-adjusted genus abundance matrix
# ============================================================

genus <- read.csv(
  file = genus_file,
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# Keep samples available in the abundance matrix
# ------------------------------------------------------------

meta.all <- meta.all[
  rownames(meta.all) %in% colnames(genus),
  ,
  drop = FALSE
]


# ============================================================
# 7. Select ADLS and HC samples
# ============================================================

meta.all2 <- meta.all %>%
  filter(
    Group2 %in% group_levels
  )


# ============================================================
# 8. Separate training and independent validation cohorts
# ============================================================

meta.train <- meta.all2 %>%
  filter(
    Project_rawID %in% train_projects
  )


meta.valid <- meta.all2 %>%
  filter(
    Project_rawID == validation_project
  )


# ============================================================
# 9. Check training and validation cohorts
# ============================================================

cat(
  "\n==============================\n"
)

cat(
  "Training cohort distribution:\n"
)

cat(
  "==============================\n"
)

print(
  table(
    meta.train$Project_rawID,
    meta.train$Group2
  )
)


cat(
  "\n==============================\n"
)

cat(
  "Independent validation cohort:\n"
)

cat(
  "==============================\n"
)

print(
  table(
    meta.valid$Project_rawID,
    meta.valid$Group2
  )
)


# ============================================================
# 10. Check sample independence
# ============================================================

overlap_samples <- intersect(
  rownames(meta.train),
  rownames(meta.valid)
)


if (length(overlap_samples) > 0) {
  
  stop(
    "Training and validation cohorts contain overlapping samples."
  )
  
} else {
  
  cat(
    "\nTraining and independent validation sets contain no overlapping samples.\n"
  )
  
}


# ============================================================
# 11. Extract training and validation abundance matrices
# ============================================================

genus_train <- genus[
  ,
  rownames(meta.train),
  drop = FALSE
]


genus_valid <- genus[
  ,
  rownames(meta.valid),
  drop = FALSE
]


# Check sample order
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


# ============================================================
# 12. Filter candidate genera according to abundance matrix
# ============================================================

diff_features <- fulljoin


diff_features <- diff_features[
  diff_features %in% rownames(genus)
]


cat(
  "\nNumber of candidate genera available in the abundance matrix:",
  length(diff_features),
  "\n"
)


if (length(diff_features) == 0) {
  
  stop(
    "No candidate genera were found in the abundance matrix."
  )
  
}


# ============================================================
# 13. Extract common features from training and validation sets
# ============================================================

abun_diff_train <- genus_train[
  rownames(genus_train) %in% diff_features,
  ,
  drop = FALSE
]


abun_diff_valid <- genus_valid[
  rownames(genus_valid) %in% diff_features,
  ,
  drop = FALSE
]


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


cat(
  "Number of common genera used for modeling:",
  length(common_features),
  "\n"
)


# ============================================================
# 14. Construct training and validation metadata
# ============================================================

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


# ============================================================
# 15. Construct modeling datasets
# ============================================================

# Rows = samples
# Columns = genera

abun_diff_train_transposed <- t(
  abun_diff_train
)


abun_diff_valid_transposed <- t(
  abun_diff_valid
)


train_data <- cbind(
  meta_train,
  abun_diff_train_transposed
)


test_data <- cbind(
  meta_valid,
  abun_diff_valid_transposed
)


# ------------------------------------------------------------
# Define outcome variable
# ------------------------------------------------------------

train_data$Group <- factor(
  train_data$Group,
  levels = group_levels
)


test_data$Group <- factor(
  test_data$Group,
  levels = group_levels
)


# ------------------------------------------------------------
# Keep Project as character
# ------------------------------------------------------------

train_data$Project <- as.character(
  train_data$Project
)


test_data$Project <- as.character(
  test_data$Project
)


# Check row order
stopifnot(
  identical(
    rownames(train_data),
    rownames(meta.train)
  )
)


# ============================================================
# 16. Study-stratified 5-fold cross-validation
# ============================================================

set.seed(seed_cv)


fold_test <- vector(
  mode = "list",
  length = K
)


names(fold_test) <- paste0(
  "Fold",
  seq_len(K)
)


# ------------------------------------------------------------
# Perform stratified CV separately within each study
# ------------------------------------------------------------

for (study in train_projects) {
  
  cat(
    "\n------------------------------------\n"
  )
  
  cat(
    "Generating 5-fold CV for cohort:",
    study,
    "\n"
  )
  
  cat(
    "------------------------------------\n"
  )
  
  
  study_idx <- which(
    train_data$Project == study
  )
  
  
  if (length(study_idx) == 0) {
    
    warning(
      "No samples found for cohort ",
      study,
      ". Skipping."
    )
    
    next
  }
  
  
  study_group <- train_data$Group[
    study_idx
  ]
  
  
  print(
    table(study_group)
  )
  
  
  if (length(study_idx) < K) {
    
    stop(
      paste0(
        "Cohort ",
        study,
        " contains fewer than ",
        K,
        " samples."
      )
    )
  }
  
  
  class_counts <- table(
    study_group
  )
  
  
  if (any(class_counts < K)) {
    
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
  }
  
  
  study_folds <- createFolds(
    y = study_group,
    k = K,
    list = TRUE,
    returnTrain = FALSE
  )
  
  
  for (k in seq_len(K)) {
    
    local_test_idx <- study_folds[[k]]
    
    
    global_test_idx <- study_idx[
      local_test_idx
    ]
    
    
    fold_test[[k]] <- c(
      fold_test[[k]],
      global_test_idx
    )
  }
}


# Sort fold indices
fold_test <- lapply(
  fold_test,
  sort
)


# ------------------------------------------------------------
# Generate training indices
# ------------------------------------------------------------

all_train_rows <- seq_len(
  nrow(train_data)
)


fold_train <- lapply(
  fold_test,
  function(test_idx) {
    
    setdiff(
      all_train_rows,
      test_idx
    )
  }
)


names(fold_train) <- names(
  fold_test
)


# ============================================================
# 17. Check CV fold assignment
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Study-stratified 5-fold CV distribution\n"
)

cat(
  "============================================\n"
)


for (k in seq_len(K)) {
  
  cat(
    "\n========== ",
    names(fold_test)[k],
    ": Validation ==========\n",
    sep = ""
  )
  
  
  fold_info <- data.frame(
    
    SampleID = rownames(train_data)[
      fold_test[[k]]
    ],
    
    Project = train_data$Project[
      fold_test[[k]]
    ],
    
    Group = train_data$Group[
      fold_test[[k]]
    ]
    
  )
  
  
  print(
    table(
      fold_info$Project,
      fold_info$Group
    )
  )
  
  
  cat(
    "\nTotal samples:",
    nrow(fold_info),
    "\n"
  )
}


# ------------------------------------------------------------
# Check whether every training sample appears exactly once
# ------------------------------------------------------------

all_validation_indices <- unlist(
  fold_test
)


if (
  length(all_validation_indices) !=
  nrow(train_data)
) {
  
  stop(
    "The total number of validation assignments does not match the number of training samples."
  )
}


if (
  length(unique(all_validation_indices)) !=
  nrow(train_data)
) {
  
  stop(
    "Some samples appear in multiple validation folds or are missing from validation."
  )
}


cat(
  "\nCheck passed: every training sample appears in exactly one validation fold.\n"
)


# ============================================================
# 18. Save fold assignment
# ============================================================

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
    output_dir,
    "study_stratified_5fold_assignment.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. Recursive Feature Elimination (RFE)
# ============================================================

rfe_control <- rfeControl(
  
  functions = rfFuncs,
  
  method = "cv",
  
  number = K,
  
  index = fold_train,
  
  verbose = FALSE,
  
  returnResamp = "final"
)


# ------------------------------------------------------------
# Define candidate feature columns
# ------------------------------------------------------------

feature_columns <- setdiff(
  
  colnames(train_data),
  
  c(
    "SampleID",
    "Group",
    "Project"
  )
)


cat(
  "\nNumber of candidate genera entering RFE:",
  length(feature_columns),
  "\n"
)


# ------------------------------------------------------------
# Candidate feature sizes
# ------------------------------------------------------------
#
# For reproducibility, the original analysis evaluates
# all possible feature sizes.
#
# If the number of genera is large, users may consider
# evaluating a smaller sequence, for example:
#
# rfe_sizes <- seq(
#   5,
#   length(feature_columns),
#   by = 5
# )
#
# ------------------------------------------------------------

rfe_sizes <- seq_len(
  length(feature_columns)
)


set.seed(seed_cv)


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


# ============================================================
# 20. Save RFE results
# ============================================================

write.csv(
  rfe_result$results,
  file = file.path(
    output_dir,
    "random_forest_RFE_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 21. Select optimal features
# ============================================================

best_features <- rfe_result$optVariables


cat(
  "\n============================================\n"
)

cat(
  "Number of features selected by RFE:",
  length(best_features),
  "\n"
)

cat(
  "============================================\n"
)


print(
  best_features
)


# Save selected features
write.csv(
  data.frame(
    genus = best_features
  ),
  file = file.path(
    output_dir,
    "random_forest_RFE_selected_genera.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. Check selected features
# ============================================================

missing_train_features <- setdiff(
  best_features,
  colnames(train_data)
)


missing_test_features <- setdiff(
  best_features,
  colnames(test_data)
)


if (
  length(missing_train_features) > 0
) {
  
  stop(
    "Some RFE-selected features are missing from the training dataset."
  )
}


if (
  length(missing_test_features) > 0
) {
  
  stop(
    paste0(
      "Some RFE-selected features are missing from the independent validation dataset: ",
      paste(
        missing_test_features,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 23. Train final Random Forest model
# ============================================================

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


set.seed(seed_model)


rf_model <- train(
  
  x = train_data[
    ,
    best_features,
    drop = FALSE
  ],
  
  y = train_data$Group,
  
  method = "rf",
  
  trControl = train_control,
  
  tuneLength = 3,
  
  metric = "ROC",
  
  importance = TRUE
)


print(
  rf_model
)


cat(
  "\nBest Random Forest parameter:\n"
)

print(
  rf_model$bestTune
)


# ============================================================
# 24. Save Random Forest model
# ============================================================

saveRDS(
  rf_model,
  file = file.path(
    output_dir,
    "random_forest_model.rds"
  )
)


# ============================================================
# 25. Internal 5-fold CV performance
# ============================================================

cv_predictions <- rf_model$pred


# ------------------------------------------------------------
# Keep predictions corresponding to the optimal mtry
# ------------------------------------------------------------

if (
  "mtry" %in% colnames(cv_predictions)
) {
  
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
    output_dir,
    "random_forest_internal_CV_predictions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 26. Internal CV ROC and AUC
# ============================================================

cv_roc <- roc(
  
  response = cv_predictions_best$obs,
  
  predictor = cv_predictions_best$ADLS,
  
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
  "\nStudy-stratified 5-fold CV AUC:",
  as.numeric(cv_auc),
  "\n"
)


# ============================================================
# 27. Determine classification threshold using Youden index
# ============================================================

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


cat(
  "\n============================================\n"
)

cat(
  "Threshold selected from training CV\n"
)

cat(
  "============================================\n"
)

print(
  best_threshold_info
)


best_threshold <- as.numeric(
  best_threshold_info["threshold"]
)


cat(
  "\nSelected threshold:",
  best_threshold,
  "\n"
)


# ============================================================
# 28. Independent validation
# ============================================================

# ------------------------------------------------------------
# Predict class probabilities
# ------------------------------------------------------------

test_probs <- predict(
  
  rf_model,
  
  newdata = test_data[
    ,
    best_features,
    drop = FALSE
  ],
  
  type = "prob"
)


# ------------------------------------------------------------
# Classify validation samples using the threshold
# determined exclusively from training CV
# ------------------------------------------------------------

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
  
  levels = group_levels
)


# ============================================================
# 29. Confusion matrix
# ============================================================

conf_matrix <- confusionMatrix(
  
  data = test_predictions,
  
  reference = test_data$Group,
  
  positive = case_group
)


print(
  conf_matrix
)


# Convert confusion matrix to data frame
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
    output_dir,
    paste0(
      validation_project,
      "_independent_validation_confusion_matrix.csv"
    )
  ),
  
  row.names = FALSE
)


# ============================================================
# 30. Independent validation performance
# ============================================================

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


# ============================================================
# 31. Independent validation ROC and AUC
# ============================================================

roc_curve <- roc(
  
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


auc_value <- auc(
  roc_curve
)


auc_ci <- ci.auc(
  roc_curve
)


# ============================================================
# 32. Save independent validation performance
# ============================================================

performance_df <- data.frame(
  
  Validation_Cohort = validation_project,
  
  Threshold = best_threshold,
  
  Accuracy = accuracy,
  
  Sensitivity = recall,
  
  Specificity = specificity,
  
  Precision = precision,
  
  F1_Score = f1_score,
  
  Balanced_Accuracy = balanced_accuracy,
  
  AUC = as.numeric(
    auc_value
  ),
  
  AUC_CI_Lower = as.numeric(
    auc_ci[1]
  ),
  
  AUC_CI_Median = as.numeric(
    auc_ci[2]
  ),
  
  AUC_CI_Upper = as.numeric(
    auc_ci[3]
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  performance_df,
  
  file = file.path(
    
    output_dir,
    
    paste0(
      validation_project,
      "_independent_validation_performance.csv"
    )
  ),
  
  row.names = FALSE
)


# ============================================================
# 33. Save individual validation predictions
# ============================================================

validation_predictions <- data.frame(
  
  SampleID = rownames(test_data),
  
  Project = test_data$Project,
  
  Observed_Group = test_data$Group,
  
  Predicted_Group = test_predictions,
  
  Probability_ADLS = test_probs[
    ,
    case_group
  ],
  
  Threshold = best_threshold,
  
  stringsAsFactors = FALSE
)


write.csv(
  
  validation_predictions,
  
  file = file.path(
    
    output_dir,
    
    paste0(
      validation_project,
      "_independent_validation_predictions.csv"
    )
  ),
  
  row.names = FALSE
)


# ============================================================
# 34. Save summary of analysis parameters
# ============================================================

analysis_parameters <- data.frame(
  
  Parameter = c(
    
    "Case group",
    
    "Control group",
    
    "Training cohorts",
    
    "Independent validation cohort",
    
    "Number of CV folds",
    
    "Candidate genera",
    
    "RFE-selected genera",
    
    "Random Forest mtry",
    
    "CV AUC",
    
    "Validation AUC",
    
    "Validation threshold"
    
  ),
  
  Value = c(
    
    case_group,
    
    control_group,
    
    paste(
      train_projects,
      collapse = ", "
    ),
    
    validation_project,
    
    K,
    
    length(common_features),
    
    length(best_features),
    
    rf_model$bestTune$mtry,
    
    as.numeric(cv_auc),
    
    as.numeric(auc_value),
    
    best_threshold
    
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  analysis_parameters,
  
  file = file.path(
    
    output_dir,
    
    "random_forest_analysis_summary.csv"
    
  ),
  
  row.names = FALSE
)


# ============================================================
# 35. Final message
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "Analysis completed successfully.\n"
)

cat(
  "============================================\n"
)

cat(
  "Training cohorts:",
  paste(
    train_projects,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Independent validation cohort:",
  validation_project,
  "\n"
)

cat(
  "Candidate genera:",
  length(common_features),
  "\n"
)

cat(
  "RFE-selected genera:",
  length(best_features),
  "\n"
)

cat(
  "Internal CV AUC:",
  round(
    as.numeric(cv_auc),
    4
  ),
  "\n"
)

cat(
  "Independent validation AUC:",
  round(
    as.numeric(auc_value),
    4
  ),
  "\n"
)

cat(
  "Results saved to:",
  output_dir,
  "\n"
)

cat(
  "============================================\n"
)
```
