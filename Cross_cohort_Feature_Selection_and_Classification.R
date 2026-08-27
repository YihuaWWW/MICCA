##软件包的安装
# 基于CRAN安装R包，检测没有则安装
# p_list = c("magrittr","dplyr","ggplot2","vegan","ggpubr","ggsci","patchwork")
# for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
#   library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
# 
# 
# # 基于github安装
# library(devtools)
# if(!requireNamespace("MMUPHin", quietly = TRUE))
#   install_github("biobakery/mmuphin@master")


# 加载R包 Load the package
suppressWarnings(suppressMessages(library(MMUPHin)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(vegan)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggsci)))
suppressWarnings(suppressMessages(library(patchwork)))

############data preprocess##########
# Load data
# 导入metadata数据
setwd("G:/博士学习/AAA/cy/AD/data")
meta.all <- read.csv(file = 'all/group.csv',stringsAsFactors = FALSE, header = TRUE, row.names = 1,
                     check.name = FALSE)

# 将 Project 列转换为 factor，并指定排序的级别
meta.all$Project <- factor(meta.all$Project, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))

# 根据 Project 列的顺序对 meta.all 数据框进行排序
meta.all <- meta.all[order(meta.all$Project), ]

# 导入细菌物种相对丰度数据
# Import relative abundance data of bacteria
feat.abu <- read.table(file = "all/tax_L6.txt", sep = "\t", header = T, check.names = FALSE)
library(dplyr)
# sum of Species
# 计算每个Species微生物相对丰度之和，避免有重复Species统计
feat.abu<-aggregate(.~ genus,data=feat.abu,sum)
rownames(feat.abu) = feat.abu$genus
feat.abu = feat.abu[, -1]
# 保留 meta.all 中行名在 feat.abu 列名中能找到的行
meta.all <- meta.all[rownames(meta.all) %in% colnames(feat.abu), ]
#write.csv(meta.all,"all/meta_all_new.csv")
# 根据新的 meta.all 的行名顺序重新排列 feat.abu 的列
feat.abu <- feat.abu[, rownames(meta.all)]
feat.abu[is.na(feat.abu)] <- 0
#feat.abu <- feat.abu/100
## 删除所有微生物丰度均为0的样本
zero.samples <- colnames(feat.abu)[
  colSums(feat.abu) == 0
]

#cat("全零样本数量：", length(zero.samples), "\n")
#print(zero.samples)

# 同时从丰度矩阵和元数据中删除
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

# 删除没有样本的因子水平
meta.all$Project <- droplevels(meta.all$Project)

# 再次按照元数据样本顺序排列丰度矩阵
feat.abu <- feat.abu[, rownames(meta.all), drop = FALSE]

# 检查丰度矩阵和元数据是否完全对应
stopifnot(
  identical(colnames(feat.abu), rownames(meta.all)),
  all(colSums(feat.abu) > 0)
)

#cat("删除后剩余样本数量：", ncol(feat.abu), "\n")
#table(meta.all$Project)
######校正批次效应#####

feat.abu_adj <- feat.abu/100
# Zero-inflated empirical Bayes adjustment of batch effect in compositional feature abundance data
# 成分特征丰度数据中批次效应的零膨胀经验贝叶斯调整
fit_adjust_batch <- adjust_batch(feature_abd = feat.abu_adj,
                                 batch = "Project",
                                 #covariates = c("Gender", "Age","Part"),
                                 data = meta.all,
                                 control = list(verbose = FALSE))
npc_abd_adj <- fit_adjust_batch$feature_abd_adj
npc_abd_adj <- npc_abd_adj*100

write.csv(
  npc_abd_adj,
  file = "G:/博士学习/AAA/cy/AD/data/all/revised_results/genus_NPC_ra_abd_adj.csv",
  row.names = TRUE
)



########ADLS&HC特征选择####
meta.all2 <- meta.all %>% 
  filter(Group2 %in% c("ADLS", "HC"))
# meta.all2 <- meta.all2 %>% 
#   filter(Project_rawID %in% c("S2", "S17","S19","S47","S71","S74","S75"))
#meta.all2 <- meta.all2 %>% 
#  filter(Project_rawID %in% c("S2", "S17","S19","S47","S71","S74"))
meta.all2 <- meta.all2 %>% 
  filter(Project_rawID %in% c("S2", "S17","S19","S47","S71","S75"))
npc_abd_adj2 <- npc_abd_adj[, colnames(npc_abd_adj) %in% rownames(meta.all2)]
new_order <- colnames(npc_abd_adj2)
meta.all2 <- meta.all2[new_order, ]

#meta.all2$Project <- factor(meta.all2$Project, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))
#meta.all2$Project <- factor(meta.all2$Project, levels = c("S1", "S2", "S3", "S7", "S8", "S9", "S10"))
#meta.all2$Project <- factor(meta.all2$Project, levels = c("S1", "S2", "S3", "S7", "S8", "S9"))
meta.all2$Project <- factor(meta.all2$Project, levels = c("S1", "S2", "S3", "S7", "S8","S10"))

meta.all2$Group2 <- factor(meta.all2$Group2, levels = c("ADLS", "HC"))


fit_lm_meta <- lm_meta(feature_abd = npc_abd_adj2,
                       batch = "Project",
                       exposure = "Group2",
                       #covariates = c("gender", "age", "BMI"),
                       data = meta.all2
)

meta_fits<-fit_lm_meta$meta_fits
write.csv(meta_fits, './all/revised_results/lm_meta_ADLS&HC_S9validation.csv')


##### 分类预测 #########

setwd("G:/博士学习/AAA/cy/AD/data")

############################################################
### 加载包
############################################################

library(randomForest)
library(caret)
library(pROC)
library(dplyr)


############################################################
### 1. 读取 metadata
############################################################

meta.all <- read.csv(
  file = "all/group.csv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


############################################################
### 2. 明确训练队列和独立验证队列
############################################################

train_projects <- c(
  "S2",
  "S17",
  "S19",
  "S47",
  "S71",
  "S75"
)

validation_project <- "S74"


# 如果后面需要按照 Project 排序
project_order <- c(
  train_projects,
  validation_project
)


meta.all$Project <- factor(
  meta.all$Project,
  levels = c(
    project_order,
    setdiff(
      unique(as.character(meta.all$Project)),
      project_order
    )
  )
)


# 根据 Project 排序
meta.all <- meta.all[
  order(meta.all$Project),
  ,
  drop = FALSE
]


############################################################
### 3. 读取 MMUPHin 批次校正后的丰度矩阵
############################################################

genus <- read.csv(
  file = "all/revised_results/genus_NPC_ra_abd_adj.csv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# 只保留丰度矩阵中存在的样本
meta.all <- meta.all[
  rownames(meta.all) %in% colnames(genus),
  ,
  drop = FALSE
]


############################################################
### 4. 只保留 ADLS 和 HC
############################################################

meta.all2 <- meta.all %>%
  filter(Group2 %in% c("ADLS", "HC"))


############################################################
### 5. 训练集与独立验证集完全分开
############################################################

# 训练集：
# S2、S17、S19、S47、S71、S75

meta.train <- meta.all2 %>%
  filter(Project_rawID %in% train_projects)


# 独立验证集：
# S74

meta.valid <- meta.all2 %>%
  filter(Project_rawID == validation_project)


############################################################
### 检查训练集和验证集
############################################################

cat("\n==============================\n")
cat("训练集样本分布：\n")
cat("==============================\n")

print(
  table(
    meta.train$Project_rawID,
    meta.train$Group2
  )
)


cat("\n==============================\n")
cat("独立验证集样本分布：\n")
cat("==============================\n")

print(
  table(
    meta.valid$Project_rawID,
    meta.valid$Group2
  )
)


############################################################
### 检查训练集和验证集是否存在重复样本
############################################################

overlap_samples <- intersect(
  rownames(meta.train),
  rownames(meta.valid)
)

if (length(overlap_samples) > 0) {
  
  stop(
    "训练集和独立验证集存在重复样本，请检查！"
  )
  
} else {
  
  cat(
    "\n训练集与独立验证集完全独立，无重复样本。\n"
  )
}


############################################################
### 6. 提取训练集和验证集丰度矩阵
############################################################

# 必须严格按照 metadata 的样本顺序排列

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


############################################################
### 再次检查样本顺序
############################################################

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


############################################################
### 7. 读取训练集筛选得到的差异菌
############################################################

diff_genus <- read.csv(
  file = "all/revised_results/lm_meta_0.05_ADLS&HC_S9validation.csv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


diff_features <- unique(
  diff_genus$feature
)


# 只保留丰度矩阵中真正存在的菌
diff_features <- diff_features[
  diff_features %in% rownames(genus)
]


cat(
  "\n用于建模的差异菌数量：",
  length(diff_features),
  "\n"
)


if (length(diff_features) == 0) {
  
  stop(
    "没有差异菌能够在 genus 丰度矩阵中找到，",
    "请检查 diff_genus$feature。"
  )
}


############################################################
### 8. 提取训练集和验证集相同的差异菌
############################################################

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


############################################################
### 保证训练集和验证集特征完全相同、顺序完全一致
############################################################

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


############################################################
### 9. 构建 metadata
############################################################

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


############################################################
### 10. 构建训练数据和独立验证数据
############################################################

# 转置：
# 行 = 样本
# 列 = 菌

abun_diff_train_transposed <- t(
  abun_diff_train
)

abun_diff_valid_transposed <- t(
  abun_diff_valid
)


############################################################
### 合并 metadata 与丰度
############################################################

train_data <- cbind(
  meta_train,
  abun_diff_train_transposed
)


test_data <- cbind(
  meta_valid,
  abun_diff_valid_transposed
)


############################################################
### 设定 Group 因子
### ADLS 为阳性类别
############################################################

train_data$Group <- factor(
  train_data$Group,
  levels = c("ADLS", "HC")
)


test_data$Group <- factor(
  test_data$Group,
  levels = c("ADLS", "HC")
)


############################################################
### Project 保留为字符，不能进入模型
############################################################

train_data$Project <- as.character(
  train_data$Project
)

test_data$Project <- as.character(
  test_data$Project
)


############################################################
### 检查 train_data 与 meta.train 行顺序
############################################################

stopifnot(
  identical(
    rownames(train_data),
    rownames(meta.train)
  )
)


############################################################
############################################################
### 11. 每个队列内部先做5折，再把相同编号的小fold合并
############################################################
############################################################

set.seed(12)

K <- 5


############################################################
### 保存最终5个 validation folds
###
### fold_test[[1]] =
###     S2的Fold1 +
###     S17的Fold1 +
###     S19的Fold1 +
###     S47的Fold1 +
###     S71的Fold1 +
###     S75的Fold1
###
### 依此类推
############################################################

fold_test <- vector(
  mode = "list",
  length = K
)

names(fold_test) <- paste0(
  "Fold",
  1:K
)

############################################################
### 对每一个训练队列分别进行5折划分
############################################################

for (study in train_projects) {
  
  cat(
    "\n------------------------------------\n"
  )
  
  cat(
    "正在对队列 ",
    study,
    " 内部进行5折划分\n",
    sep = ""
  )
  
  cat(
    "------------------------------------\n"
  )
  
  
  ##########################################################
  ### 找到该队列在 train_data 中的行号
  ##########################################################
  
  study_idx <- which(
    train_data$Project == study
  )
  
  
  if (length(study_idx) == 0) {
    
    warning(
      paste0(
        "训练集中没有找到队列 ",
        study,
        "，已跳过。"
      )
    )
    
    next
  }
  
  
  ##########################################################
  ### 该队列疾病分组
  ##########################################################
  
  study_group <- train_data$Group[
    study_idx
  ]
  
  
  ##########################################################
  ### 打印该队列 ADLS / HC 数量
  ##########################################################
  
  print(
    table(study_group)
  )
  
  
  ##########################################################
  ### 检查是否至少有5个样本
  ##########################################################
  
  if (length(study_idx) < K) {
    
    stop(
      paste0(
        "队列 ",
        study,
        " 总样本数少于5，无法进行5折交叉验证。"
      )
    )
  }
  
  
  ##########################################################
  ### 如果某一组少于5个样本：
  ### createFolds仍可能运行，
  ### 但无法保证5个fold中都有该类别
  ##########################################################
  
  class_counts <- table(
    study_group
  )
  
  
  if (any(class_counts < K)) {
    
    warning(
      paste0(
        "队列 ",
        study,
        " 中至少一个疾病类别样本数少于5。",
        "因此无法保证该队列的每一个fold中",
        "都同时存在ADLS和HC。"
      )
    )
  }
  
  
  ##########################################################
  ### 在该队列内部：
  ### 根据 ADLS / HC 进行分层5折
  ##########################################################
  
  study_folds <- createFolds(
    y = study_group,
    k = K,
    list = TRUE,
    returnTrain = FALSE
  )
  
  
  ##########################################################
  ### 把队列内部的局部行号
  ### 转换成整个 train_data 的全局行号
  ###
  ### 然后把相同编号的小fold合并
  ##########################################################
  
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


############################################################
### 对每个fold排序
############################################################

fold_test <- lapply(
  fold_test,
  sort
)


############################################################
### 生成对应 training indices
###
### caret 的 index = 每一折用于训练的样本
### indexOut = 每一折用于验证的样本
############################################################

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


############################################################
############################################################
### 12. 检查五折划分结果
############################################################
############################################################

cat(
  "\n\n============================================\n"
)

cat(
  "Study-stratified 5-fold CV 分布检查\n"
)

cat(
  "============================================\n"
)


for (k in seq_len(K)) {
  
  cat(
    "\n========== ",
    names(fold_test)[k],
    "：Validation ==========\n",
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
    "\n该fold总样本数：",
    nrow(fold_info),
    "\n"
  )
}


############################################################
### 检查：
### 所有样本是否恰好作为 validation 出现一次
############################################################

all_validation_indices <- unlist(
  fold_test
)


cat(
  "\n============================================\n"
)

cat(
  "每个样本进入validation fold的次数：\n"
)

print(
  table(
    table(all_validation_indices)
  )
)


if (
  length(all_validation_indices) !=
  nrow(train_data)
) {
  
  stop(
    "五折validation样本总数与训练集总样本数不一致，",
    "请检查fold划分。"
  )
}


if (
  length(unique(all_validation_indices)) !=
  nrow(train_data)
) {
  
  stop(
    "存在样本重复进入多个validation fold，",
    "或存在样本未进入任何validation fold。"
  )
}


cat(
  "\n检查通过：每个训练样本恰好进入一个validation fold。\n"
)


############################################################
### 进一步打印每个队列在各Fold中的样本数
############################################################

fold_project_distribution <- data.frame()


for (k in seq_len(K)) {
  
  tmp <- as.data.frame(
    table(
      Project = train_data$Project[
        fold_test[[k]]
      ]
    )
  )
  
  
  tmp$Fold <- paste0(
    "Fold",
    k
  )
  
  
  fold_project_distribution <- rbind(
    fold_project_distribution,
    tmp
  )
}


cat(
  "\n============================================\n"
)

cat(
  "各队列在五个validation folds中的样本数：\n"
)

cat(
  "============================================\n"
)


print(
  fold_project_distribution
)


############################################################
### 各Fold疾病类别比例
############################################################

cat(
  "\n============================================\n"
)

cat(
  "各validation fold的ADLS/HC分布：\n"
)

cat(
  "============================================\n"
)


for (k in seq_len(K)) {
  
  cat(
    "\nFold ",
    k,
    ":\n",
    sep = ""
  )
  
  print(
    table(
      train_data$Group[
        fold_test[[k]]
      ]
    )
  )
}


############################################################
############################################################
### 13. RFE
###
### 使用上面完全相同的 study-stratified folds
############################################################
############################################################

rfe_control <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = K,
  
  # 关键：
  # 自定义每一折的 training indices
  index = fold_train,
  
  verbose = FALSE,
  returnResamp = "final"
)


############################################################
### 模型输入变量
###
### 注意：
### SampleID、Group、Project
### 均不能作为微生物特征进入模型
############################################################

feature_columns <- setdiff(
  colnames(train_data),
  c(
    "SampleID",
    "Group",
    "Project"
  )
)


cat(
  "\n进入RFE的候选菌数量：",
  length(feature_columns),
  "\n"
)


############################################################
### 如果特征很多，不建议 sizes = 1:n
### 这里先保持与你原来的分析逻辑一致
############################################################

rfe_sizes <- seq_len(
  length(feature_columns)
)


set.seed(12)


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


############################################################
### 查看RFE结果
############################################################

print(
  rfe_result
)


############################################################
### 14. 使用RFE筛选的最佳特征
############################################################

best_features <- rfe_result$optVariables

# 如果希望固定top10，可改为：
# best_features <- rfe_result$optVariables[1:10]

# 如果希望固定top20，可改为：
# best_features <- rfe_result$optVariables[1:20]

# 如果希望使用全部差异菌：
# best_features <- colnames(train_data)[-c(1,2)]

cat(
  "\n============================================\n"
)

cat(
  "RFE最终选择的菌数量：",
  length(best_features),
  "\n"
)

cat(
  "============================================\n"
)


print(
  best_features
)


############################################################
### 检查最佳特征在训练集和验证集中均存在
############################################################

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
    "部分RFE特征不在训练集中。"
  )
}


if (
  length(missing_test_features) > 0
) {
  
  stop(
    paste0(
      "部分RFE特征不在独立验证集中：",
      paste(
        missing_test_features,
        collapse = ", "
      )
    )
  )
}


############################################################
############################################################
### 15. 最终随机森林：
###
### 同样使用 study-stratified 5-fold CV
############################################################
############################################################

train_control <- trainControl(
  
  method = "cv",
  
  number = K,
  
  # 每折训练样本
  index = fold_train,
  
  # 每折验证样本
  indexOut = fold_test,
  
  classProbs = TRUE,
  
  savePredictions = "final",
  
  summaryFunction = twoClassSummary,
  
  verboseIter = FALSE
)


############################################################
### 使用全部训练数据建立最终随机森林模型
###
### CV用于：
### 1. 超参数选择
### 2. 内部性能评估
###
### 最终选定参数后，
### caret会用全部训练数据重新拟合最终模型
############################################################

set.seed(12)


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


############################################################
### 查看随机森林结果
############################################################

print(
  rf_model
)


############################################################
### 查看最佳 mtry
############################################################

cat(
  "\n最佳随机森林参数：\n"
)

print(
  rf_model$bestTune
)


############################################################
### 16. 查看内部五折CV预测结果
############################################################

cv_predictions <- rf_model$pred


cat(
  "\n内部五折CV预测结果前几行：\n"
)

print(
  head(cv_predictions)
)


############################################################
### 如果存在多个mtry，
### 只保留最终最佳mtry对应预测
############################################################

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


# ############################################################
# ### 内部5-fold CV ROC
# ###
# ### ADLS是阳性类别，
# ### 必须明确使用ADLS概率
# ############################################################
# 
# cv_roc <- roc(
#   
#   response = cv_predictions_best$obs,
#   
#   predictor = cv_predictions_best$ADLS,
#   
#   levels = c(
#     "HC",
#     "ADLS"
#   ),
#   
#   direction = "<",
#   
#   quiet = TRUE
# )
# 
# 
# cv_auc <- auc(
#   cv_roc
# )
# 
# 
# cat(
#   "\nStudy-stratified 5-fold CV AUC: ",
#   as.numeric(cv_auc),
#   "\n",
#   sep = ""
# )


############################################################
### 根据训练集内部 study-stratified 5-fold CV
### 确定最佳分类阈值
###
### 注意：
### threshold只能由训练集内部CV确定
### 不能使用S74真实标签优化threshold
############################################################

cv_roc <- roc(
  response = cv_predictions_best$obs,
  predictor = cv_predictions_best$ADLS,
  levels = c("HC", "ADLS"),
  direction = "<",
  quiet = TRUE
)


############################################################
### 内部5-fold CV AUC
############################################################

cv_auc <- auc(cv_roc)

cat(
  "\nStudy-stratified 5-fold CV AUC: ",
  as.numeric(cv_auc),
  "\n",
  sep = ""
)


############################################################
### 使用Youden index寻找最佳threshold
###
### 最大化：
### Sensitivity + Specificity - 1
############################################################

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

print(best_threshold_info)


############################################################
### 提取最佳threshold
############################################################

best_threshold <- as.numeric(
  best_threshold_info["threshold"]
)


cat(
  "\nSelected threshold: ",
  best_threshold,
  "\n",
  sep = ""
)





############################################################
############################################################
### 17. S74 独立验证
############################################################
############################################################

test_predictions <- predict(
  
  rf_model,
  
  newdata = test_data[
    ,
    best_features,
    drop = FALSE
  ]
)


############################################################
### S74独立验证预测概率
############################################################

test_probs <- predict(
  rf_model,
  newdata = test_data[
    ,
    best_features,
    drop = FALSE
  ],
  type = "prob"
)


############################################################
### 使用训练CV确定的threshold
### 对S74进行分类
############################################################

test_predictions <- ifelse(
  test_probs[, "ADLS"] >= best_threshold,
  "ADLS",
  "HC"
)


test_predictions <- factor(
  test_predictions,
  levels = c("ADLS", "HC")
)


############################################################
### 混淆矩阵
### ADLS作为阳性类别
############################################################

conf_matrix <- confusionMatrix(
  
  data = test_predictions,
  
  reference = test_data$Group,
  
  positive = "ADLS"
)


print(
  conf_matrix
)

# confusionMatrix$table 转为 data.frame
conf_matrix_df <- as.data.frame(
  conf_matrix$table
)

# 重命名列，便于阅读
colnames(conf_matrix_df) <- c(
  "Prediction",
  "Reference",
  "Frequency"
)

# 保存
write.csv(
  conf_matrix_df,
  file = "all/revised_results/S74_independent_validation_confusion_matrix.csv",
  row.names = FALSE
)

# write.csv(
#   conf_matrix_df,
#   file = "all/revised_results/S74_independent_validation_confusion_matrix_top10.csv",
#   row.names = FALSE
# )

# write.csv(
#   conf_matrix_df,
#   file = "all/revised_results/S74_independent_validation_confusion_matrix_top20.csv",
#   row.names = FALSE
# )

# write.csv(
#   conf_matrix_df,
#   file = "all/revised_results/S74_independent_validation_confusion_matrix_all.csv",
#   row.names = FALSE
# )

############################################################
### 18. 独立验证集性能指标
############################################################

accuracy <- unname(
  conf_matrix$overall["Accuracy"]
)

precision <- unname(
  conf_matrix$byClass["Pos Pred Value"]
)

recall <- unname(
  conf_matrix$byClass["Sensitivity"]
)

specificity <- unname(
  conf_matrix$byClass["Specificity"]
)

balanced_accuracy <- unname(
  conf_matrix$byClass["Balanced Accuracy"]
)


f1_score <- ifelse(
  is.na(precision) |
    is.na(recall) |
    precision + recall == 0,
  NA,
  2 * precision * recall /
    (precision + recall)
)

#获取独立验证集预测概率

test_probs <- predict(
  
  rf_model,
  
  newdata = test_data[
    ,
    best_features,
    drop = FALSE
  ],
  
  type = "prob"
)



### 独立验证 ROC
roc_curve <- roc(test_data$Group, test_probs[, 2])  # 假设 "Group2" 是阳性类别


### 独立验证AUC
auc_value <- auc(
  roc_curve
)


###可选：输出95% CI
auc_ci <- ci.auc(
  roc_curve
)


performance_df <- data.frame(
  Validation_Cohort = validation_project,
  
  Threshold = best_threshold,
  
  Accuracy = accuracy,
  
  Sensitivity = recall,
  
  Specificity = specificity,
  
  Precision = precision,
  
  F1_Score = f1_score,
  
  Balanced_Accuracy = balanced_accuracy,
  
  AUC = as.numeric(auc_value),
  
  stringsAsFactors = FALSE
)


write.csv(
  performance_df,
  file = paste0(
    "all/revised_results/",
    validation_project,
    "_independent_validation_performance.csv"
  ),
  row.names = FALSE
)

# write.csv(
#   performance_df,
#   file = paste0(
#     "all/revised_results/",
#     validation_project,
#     "_independent_validation_performance_top10.csv"
#   ),
#   row.names = FALSE
# )

# write.csv(
#   performance_df,
#   file = paste0(
#     "all/revised_results/",
#     validation_project,
#     "_independent_validation_performance_all.csv"
#   ),
#   row.names = FALSE
# )

write.csv(
  performance_df,
  file = paste0(
    "all/revised_results/",
    validation_project,
    "_independent_validation_performance_top20.csv"
  ),
  row.names = FALSE
)
############################################################

############################################################
### 可选：
### 保存每个训练样本属于哪个validation fold
### 方便论文补充材料和审稿回复
############################################################

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


############################################################
### 按 Project、Fold 排序查看
############################################################

fold_assignment <- fold_assignment %>%
  arrange(
    Project,
    Fold,
    Group
  )


cat(
  "\n最终fold assignment：\n"
)

print(
  head(
    fold_assignment,
    30
  )
)


############################################################
### 保存fold assignment
############################################################

write.csv(
  
  fold_assignment,
  
  file =
    "all/revised_results/study_stratified_5fold_assignment.csv",
  
  row.names = FALSE
)


cat(
  "\n分析完成。\n"
)


