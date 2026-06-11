#!/usr/bin/env Rscript

# /public/software/VersionHub/R/4.4.3/Rscript
# 单特征循环分析：丰度 ~ 阶段 + (1 | subject_id)
# 修改：将特征名转换为代称（如physiology_feature_1），并输出原始特征名与代称的映射表到主输出目录


library(dplyr)
library(Maaslin2)
library(tidyr)

# ————————————————————————————————————————————————————————
# 1. 读入数据并重命名特征
# ————————————————————————————————————————————————————————

files <- c(
  physiology = "data/processed/physiology.csv",
  eye = "data/processed/eye.csv", 
  cognitive = "data/processed/cognitive.csv",
  ultrasound = "data/processed/ultrasound.csv"
)

# 函数：创建特征名映射并重命名特征列
create_feature_mapping_and_rename <- function(data, dataset_name) {
  # 获取特征列（排除subject_id和time）
  feature_cols <- names(data)[!names(data) %in% c("subject_id", "time")]
  
  # 创建代称（例如：physiology_feature_1）
  new_feature_names <- paste0(dataset_name, "_feature_", seq_along(feature_cols))
  
  # 创建映射表
  mapping <- data.frame(
    original_name = feature_cols,
    renamed_name = new_feature_names
  )
  
  write.csv(mapping, 
            file.path("data/processed", 
                      paste0(dataset_name, "_feature_mapping.csv")), 
            row.names = FALSE)

  # 重命名数据中的特征列
  names(data)[names(data) %in% feature_cols] <- new_feature_names
  
  return(list(data = data, mapping = mapping))
}

datasets <- lapply(names(files), function(dataset_name) {
  path <- files[dataset_name]
  data <- read.csv(path, check.names = FALSE)
  data[, -(1:2)] <- lapply(data[, -(1:2)], as.numeric)

  data <- data[
    (data$subject_id %in% c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13")) &
    (data$time %in% c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "Ta", "Tb")), 
  ]
  
  # 创建行名
  rownames(data) <- paste(data$subject_id, data$time, sep = "_")
  
  # 重新编码时间
  data <- data %>%
    mutate(time = case_when(
      time == "T1" ~ "before",
      time %in% c("T9", "T10") ~ "after", 
      TRUE ~ "climb"
    ))
  
  # 重命名特征并生成映射表
  result <- create_feature_mapping_and_rename(data, dataset_name)
  return(result)
})

# 提取数据集和映射表
dataset_mappings <- lapply(datasets, function(x) x$mapping)
datasets <- lapply(datasets, function(x) x$data)
names(datasets) <- names(files)
list2env(datasets, envir = .GlobalEnv)

# ————————————————————————————————————————————————————————
# 2. 设置输出路径和辅助函数
# ————————————————————————————————————————————————————————

output_dir <- "data/processed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 保存并显示特征名映射表到主输出目录
for (dataset_name in names(dataset_mappings)) {
  cat(sprintf("\n=== 特征名映射表 for %s ===\n", dataset_name))
  print(dataset_mappings[[dataset_name]])
  write.csv(dataset_mappings[[dataset_name]], 
            file.path(output_dir, paste0(dataset_name, "_feature_mapping.csv")), 
            row.names = FALSE)
  cat(sprintf("映射表已保存至: %s\n", file.path(output_dir, paste0(dataset_name, "_feature_mapping.csv"))))
}

# 单个特征分析函数
analyze_single_feature <- function(data, feature_name, output_prefix, mapping) {
  
  # 提取单个特征的数据
  feature_data <- data[, c("subject_id", "time", feature_name), drop = FALSE]
  
  # 删除NA值
  feature_data <- feature_data[complete.cases(feature_data), ]
  
  # 检查是否有足够的数据
  if (nrow(feature_data) < 3) {
    cat(sprintf("跳过特征 %s: 数据不足\n", feature_name))
    return(NULL)
  }
  
  # 检查可用阶段
  available_stages <- unique(feature_data$time)
  if (length(available_stages) < 2) {
    cat(sprintf("跳过特征 %s: 阶段不足（只有%d个阶段）\n", feature_name, length(available_stages)))
    return(NULL)
  }
  
  cat(sprintf("    特征 %s 包含阶段: %s\n", feature_name, paste(available_stages, collapse = ", ")))
  
  # 准备丰度表和元数据表
  abundance_table <- feature_data[, feature_name, drop = FALSE]
  rownames(abundance_table) <- rownames(feature_data)
  abundance_table <- t(abundance_table)  # Maaslin2需要特征为行，样本为列
  
  metadata <- feature_data[, c("subject_id", "time"), drop = FALSE]
  rownames(metadata) <- rownames(feature_data)
  
  results_list <- list()
  
  # 根据可用阶段决定分析策略
  has_before <- "before" %in% available_stages
  has_climb <- "climb" %in% available_stages  
  has_after <- "after" %in% available_stages
  
  # 分析1: 如果有before，以before作为参考
  if (has_before) {
    ref_levels <- c("before")
    if (has_climb) ref_levels <- c(ref_levels, "climb")
    if (has_after) ref_levels <- c(ref_levels, "after")
    
    metadata$time_ref_before <- factor(metadata$time, levels = ref_levels)
    
    tryCatch({
      result_before <- Maaslin2(
        input_data = abundance_table,
        input_metadata = metadata,
        output = file.path(output_prefix, paste0(feature_name, "_ref_before")),
        fixed_effects = "time_ref_before",
        random_effects = "subject_id",
        normalization = "NONE",
        transform = "NONE",
        analysis_method = "LM",
        max_significance = 0.25,
        min_abundance = 0,
        min_prevalence = 0
      )
      
      if (!is.null(result_before$results)) {
        result_before$results$comparison_type <- "ref_before"
        result_before$results$original_feature <- feature_name
        results_list[["ref_before"]] <- result_before$results
        cat(sprintf("    特征 %s (ref_before) 分析成功\n", feature_name))
      }
    }, error = function(e) {
      cat(sprintf("    特征 %s (ref_before) 分析失败: %s\n", feature_name, e$message))
    })
  }
  
  # 分析2: 如果有after，以after作为参考
  if (has_after) {
    ref_levels <- c("after")
    if (has_climb) ref_levels <- c(ref_levels, "climb")
    if (has_before) ref_levels <- c(ref_levels, "before")
    
    metadata$time_ref_after <- factor(metadata$time, levels = ref_levels)
    
    tryCatch({
      result_after <- Maaslin2(
        input_data = abundance_table,
        input_metadata = metadata, 
        output = file.path(output_prefix, paste0(feature_name, "_ref_after")),
        fixed_effects = "time_ref_after",
        random_effects = "subject_id",
        normalization = "NONE",
        transform = "NONE", 
        analysis_method = "LM",
        max_significance = 0.25,
        min_abundance = 0,
        min_prevalence = 0
      )
      
      if (!is.null(result_after$results)) {
        result_after$results$comparison_type <- "ref_after"  
        result_after$results$original_feature <- feature_name
        results_list[["ref_after"]] <- result_after$results
        cat(sprintf("    特征 %s (ref_after) 分析成功\n", feature_name))
      }
    }, error = function(e) {
      cat(sprintf("    特征 %s (ref_after) 分析失败: %s\n", feature_name, e$message))
    })
  }
  
  # 分析3: 如果只有before和after（没有climb），进行直接比较
  if (has_before && has_after && !has_climb) {
    metadata$time_before_after <- factor(metadata$time, levels = c("before", "after"))
    
    tryCatch({
      result_before_after <- Maaslin2(
        input_data = abundance_table,
        input_metadata = metadata,
        output = file.path(output_prefix, paste0(feature_name, "_before_vs_after")),
        fixed_effects = "time_before_after",
        random_effects = "subject_id",
        normalization = "NONE",
        transform = "NONE",
        analysis_method = "LM",
        max_significance = 0.25,
        min_abundance = 0,
        min_prevalence = 0
      )
      
      if (!is.null(result_before_after$results)) {
        result_before_after$results$comparison_type <- "before_vs_after"
        result_before_after$results$original_feature <- feature_name
        results_list[["before_vs_after"]] <- result_before_after$results
        cat(sprintf("    特征 %s (before_vs_after) 分析成功\n", feature_name))
      }
    }, error = function(e) {
      cat(sprintf("    特征 %s (before_vs_after) 分析失败: %s\n", feature_name, e$message))
    })
  }
  
  return(results_list)
}

# 计算统计摘要函数
calculate_summary <- function(results_df, dataset_name) {
  if (is.null(results_df) || nrow(results_df) == 0) {
    return(data.frame(
      dataset = dataset_name,
      total_features = 0,
      sig_q0.05 = 0,
      sig_q0.1 = 0,
      sig_q0.25 = 0
    ))
  }
  
  summary_stats <- data.frame(
    dataset = dataset_name,
    total_features = length(unique(results_df$original_feature)),
    sig_q0.05 = sum(results_df$qval <= 0.05, na.rm = TRUE),
    sig_q0.1 = sum(results_df$qval <= 0.1, na.rm = TRUE), 
    sig_q0.25 = sum(results_df$qval <= 0.25, na.rm = TRUE)
  )
  
  return(summary_stats)
}

# ————————————————————————————————————————————————————————
# 3. 主要分析循环
# ————————————————————————————————————————————————————————

dataset_names <- names(datasets)
all_summaries <- list()

for (dataset_name in dataset_names) {
  
  cat(sprintf("\n=== 开始分析数据集: %s ===\n", dataset_name))
  
  # 获取当前数据集和映射表
  current_data <- get(dataset_name)
  current_mapping <- dataset_mappings[[dataset_name]]
  
  # 获取特征列（排除subject_id和time）
  feature_cols <- names(current_data)[!names(current_data) %in% c("subject_id", "time")]
  
  if (length(feature_cols) == 0) {
    cat(sprintf("数据集 %s 没有特征列，跳过\n", dataset_name))
    next
  }
  
  # 创建输出目录
  dataset_output_dir <- file.path(output_dir, dataset_name)
  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 存储所有特征的结果
  all_results_before <- list()
  all_results_after <- list()
  all_results_before_after <- list()
  
  # 对每个特征进行分析
  for (i in seq_along(feature_cols)) {
    feature_name <- feature_cols[i]
    cat(sprintf("  分析特征 %d/%d: %s\n", i, length(feature_cols), feature_name))
    
    feature_results <- analyze_single_feature(
      current_data, 
      feature_name, 
      dataset_output_dir,
      current_mapping
    )
    
    if (!is.null(feature_results)) {
      if ("ref_before" %in% names(feature_results)) {
        all_results_before[[feature_name]] <- feature_results[["ref_before"]]
      }
      if ("ref_after" %in% names(feature_results)) {
        all_results_after[[feature_name]] <- feature_results[["ref_after"]]
      }
      if ("before_vs_after" %in% names(feature_results)) {
        all_results_before_after[[feature_name]] <- feature_results[["before_vs_after"]]
      }
    }
  }
  
  # 合并结果
  if (length(all_results_before) > 0) {
    combined_results_before <- do.call(rbind, all_results_before)
    write.csv(combined_results_before, 
              file.path(dataset_output_dir, paste0(dataset_name, "_combined_results_ref_before.csv")), 
              row.names = FALSE)
  }
  
  if (length(all_results_after) > 0) {
    combined_results_after <- do.call(rbind, all_results_after)
    write.csv(combined_results_after,
              file.path(dataset_output_dir, paste0(dataset_name, "_combined_results_ref_after.csv")),
              row.names = FALSE)
  }
  
  if (length(all_results_before_after) > 0) {
    combined_results_before_after <- do.call(rbind, all_results_before_after)
    write.csv(combined_results_before_after,
              file.path(dataset_output_dir, paste0(dataset_name, "_combined_results_before_vs_after.csv")),
              row.names = FALSE)
  }
  
  # 计算统计摘要
  if (length(all_results_before) > 0) {
    summary_before <- calculate_summary(combined_results_before, paste0(dataset_name, "_ref_before"))
    all_summaries[[paste0(dataset_name, "_ref_before")]] <- summary_before
  }
  
  if (length(all_results_after) > 0) {
    summary_after <- calculate_summary(combined_results_after, paste0(dataset_name, "_ref_after"))  
    all_summaries[[paste0(dataset_name, "_ref_after")]] <- summary_after
  }
  
  if (length(all_results_before_after) > 0) {
    summary_before_after <- calculate_summary(combined_results_before_after, paste0(dataset_name, "_before_vs_after"))  
    all_summaries[[paste0(dataset_name, "_before_vs_after")]] <- summary_before_after
  }
  
  cat(sprintf("=== 完成数据集: %s ===\n", dataset_name))
}

# ————————————————————————————————————————————————————————
# 4. 生成统计摘要报告
# ————————————————————————————————————————————————————————

if (length(all_summaries) > 0) {
  final_summary <- do.call(rbind, all_summaries)
  write.csv(final_summary, 
            file.path(output_dir, "analysis_summary.csv"), 
            row.names = FALSE)
  
  cat("\n=== 分析统计摘要 ===\n")
  print(final_summary)
}

cat("\n=== 所有分析完成 ===\n")
cat(sprintf("结果保存在: %s\n", output_dir))