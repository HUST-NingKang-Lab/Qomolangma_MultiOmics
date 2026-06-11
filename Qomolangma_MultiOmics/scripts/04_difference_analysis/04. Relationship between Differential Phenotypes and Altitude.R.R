
library(dplyr)
library(tidyr)

# —————————————————————————————————————————————
# 01. 读取并处理数据
# —————————————————————————————————————————————
people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 读取差异分析结果并过滤
physiology_DE <- read.csv("data/processed/physiology_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)
eye_DE <- read.csv("data/processed/eye_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)
cognitive_DE <- read.csv("data/processed/cognitive_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)
ultrasound_DE <- read.csv("data/processed/ultrasound_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)

# 读取表型数据并过滤
physiology_abs <- read.csv("data/processed/physiology.csv", check.names = FALSE) %>%
  select(subject_id, time, all_of(unique(physiology_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12)
eye_abs <- read.csv("data/processed/eye.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(eye_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12) 
cognitive_abs <- read.csv("data/processed/cognitive.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(cognitive_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12)
ultrasound_abs <- read.csv("data/processed/ultrasound.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(ultrasound_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12)

# 读取时间和海拔对应关系
pair_time_altitude <- read.csv("data/metadata/Matedata_Information.csv", check.names = FALSE) %>% 
    select(time12, altitude) %>% 
    unique()
colnames(pair_time_altitude)[1] <- "time"

# —————————————————————————————————————————————
# 02. 定义相关性分析函数
# —————————————————————————————————————————————
calculate_correlation <- function(data, data_type) {
  # 获取表型列名（排除subject_id和time列）
  phenotype_cols <- setdiff(colnames(data), c("subject_id", "time"))
  
  # 存储结果的列表
  results_list <- list()
  
  # 对每个人进行分析
  for (person in people) {
    person_data <- data %>% filter(subject_id == person)
    
    # 如果该人没有数据，跳过
    if (nrow(person_data) == 0) next
    
    # 合并海拔数据
    person_data_with_altitude <- person_data %>%
      left_join(pair_time_altitude, by = "time") %>%
      filter(!is.na(altitude))  # 移除没有海拔信息的行
    
    # 如果合并后数据不足，跳过
    if (nrow(person_data_with_altitude) < 3) next
    
    # 对每个表型计算相关性
    for (phenotype in phenotype_cols) {
      phenotype_data <- person_data_with_altitude %>%
        select(time, altitude, !!sym(phenotype)) %>%
        filter(!is.na(!!sym(phenotype)))  # 移除表型缺失值
      
      # 如果数据点不足，跳过
      if (nrow(phenotype_data) < 3) next
      
      # 计算Pearson相关性（趋势相关性）
      cor_result <- cor.test(phenotype_data$altitude, phenotype_data[[phenotype]], 
                           method = "pearson")
      
      # 计算Spearman相关性（作为补充，更稳健）
      spearman_result <- cor.test(phenotype_data$altitude, phenotype_data[[phenotype]], 
                                method = "spearman")
      
      # 存储结果
      results_list[[length(results_list) + 1]] <- data.frame(
        data_type = data_type,
        subject_id = person,
        phenotype = phenotype,
        n_points = nrow(phenotype_data),
        pearson_correlation = cor_result$estimate,
        pearson_pvalue = cor_result$p.value,
        spearman_correlation = spearman_result$estimate,
        spearman_pvalue = spearman_result$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 如果有结果，合并返回
  if (length(results_list) > 0) {
    return(do.call(rbind, results_list))
  } else {
    return(NULL)
  }
}

# —————————————————————————————————————————————
# 03. 对每种数据类型进行相关性分析
# —————————————————————————————————————————————
cat("正在分析生理表型...\n")
physiology_results <- calculate_correlation(physiology_abs, "physiology")

cat("正在分析眼部表型...\n")
eye_results <- calculate_correlation(eye_abs, "eye")

cat("正在分析认知表型...\n")
cognitive_results <- calculate_correlation(cognitive_abs, "cognitive")

cat("正在分析超声表型...\n")
ultrasound_results <- calculate_correlation(ultrasound_abs, "ultrasound")

# —————————————————————————————————————————————
# 04. 合并所有结果
# —————————————————————————————————————————————
all_results <- rbind(
  physiology_results,
  eye_results,
  cognitive_results,
  ultrasound_results
)

# 添加行名
rownames(all_results) <- NULL

# 对结果按p值排序
all_results <- all_results %>%
  arrange(pearson_pvalue)

# 添加FDR校正
all_results$pearson_fdr <- p.adjust(all_results$pearson_pvalue, method = "fdr")
all_results$spearman_fdr <- p.adjust(all_results$spearman_pvalue, method = "fdr")

# —————————————————————————————————————————————
# 05. 保存结果
# —————————————————————————————————————————————
# 创建输出目录
output_dir <- "data/metadata/Matedata_Information.csv"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存完整结果
write.csv(all_results, 
          file.path(output_dir, "phenotype_altitude_correlation_results.csv"), 
          row.names = FALSE)

# 保存显著结果（Pearson相关性p<0.05）
significant_results <- all_results %>%
  filter(pearson_pvalue < 0.05)

write.csv(significant_results, 
          file.path(output_dir, "significant_phenotype_altitude_correlations.csv"), 
          row.names = FALSE)

# —————————————————————————————————————————————
# 06. 计算每个表型在所有个体中的平均相关性
# —————————————————————————————————————————————
cat("正在计算表型平均相关性...\n")

# 按表型和数据类型分组，计算平均相关性
phenotype_mean_correlation <- all_results %>%
  group_by(data_type, phenotype) %>%
  summarise(
    n_subjects = n(),  # 统计的个体数量
    mean_pearson_correlation = mean(pearson_correlation, na.rm = TRUE),
    mean_spearman_correlation = mean(spearman_correlation, na.rm = TRUE),
    sd_pearson_correlation = sd(pearson_correlation, na.rm = TRUE),
    sd_spearman_correlation = sd(spearman_correlation, na.rm = TRUE),
    mean_n_points = mean(n_points, na.rm = TRUE),  # 平均数据点数
    .groups = 'drop'
  )

# 对平均相关性进行显著性检验（单样本t检验，检验是否显著不为0）
phenotype_mean_correlation$pearson_ttest_pvalue <- NA
phenotype_mean_correlation$spearman_ttest_pvalue <- NA

for (i in 1:nrow(phenotype_mean_correlation)) {
  current_data_type <- phenotype_mean_correlation$data_type[i]
  current_phenotype <- phenotype_mean_correlation$phenotype[i]
  
  # 获取该表型的所有个体相关性数据
  phenotype_correlations <- all_results %>%
    filter(data_type == current_data_type, phenotype == current_phenotype)
  
  # 对Pearson相关性进行单样本t检验
  if (nrow(phenotype_correlations) > 1) {
    pearson_ttest <- t.test(phenotype_correlations$pearson_correlation, mu = 0)
    phenotype_mean_correlation$pearson_ttest_pvalue[i] <- pearson_ttest$p.value
    
    spearman_ttest <- t.test(phenotype_correlations$spearman_correlation, mu = 0)
    phenotype_mean_correlation$spearman_ttest_pvalue[i] <- spearman_ttest$p.value
  }
}

# FDR校正
phenotype_mean_correlation$pearson_ttest_fdr <- p.adjust(phenotype_mean_correlation$pearson_ttest_pvalue, method = "fdr")
phenotype_mean_correlation$spearman_ttest_fdr <- p.adjust(phenotype_mean_correlation$spearman_ttest_pvalue, method = "fdr")

# 按平均相关性的绝对值排序
phenotype_mean_correlation <- phenotype_mean_correlation %>%
  arrange(desc(abs(mean_pearson_correlation)))

# 筛选显著结果（基于t检验p<0.05）
significant_phenotype_mean <- phenotype_mean_correlation %>%
  filter(pearson_ttest_pvalue < 0.05 | spearman_ttest_pvalue < 0.05)

# —————————————————————————————————————————————
# 07. 保存结果
# —————————————————————————————————————————————

# 保存表型平均相关性结果
write.csv(phenotype_mean_correlation, 
          file.path(output_dir, "phenotype_mean_correlation_results.csv"), 
          row.names = FALSE)

# 保存显著的表型平均相关性结果
write.csv(significant_phenotype_mean, 
          file.path(output_dir, "significant_phenotype_mean_correlations.csv"), 
          row.names = FALSE)

# —————————————————————————————————————————————
# 08. 结果汇总
# —————————————————————————————————————————————
cat("\n=== 分析完成 ===\n")
cat("个体水平分析：\n")
cat("  总共分析了", nrow(all_results), "个表型-个体组合\n")
cat("  显著相关（p<0.05）的有", nrow(significant_results), "个\n")
cat("  FDR校正后显著（FDR<0.05）的有", sum(all_results$pearson_fdr < 0.05, na.rm = TRUE), "个\n")

cat("\n表型水平分析：\n")
cat("  总共分析了", nrow(phenotype_mean_correlation), "个表型\n")
cat("  平均相关性显著（t检验p<0.05）的有", nrow(significant_phenotype_mean), "个\n")
cat("  FDR校正后显著（FDR<0.05）的有", 
    sum(phenotype_mean_correlation$pearson_ttest_fdr < 0.05, na.rm = TRUE), "个\n")

# 按数据类型统计（个体水平）
summary_by_type <- all_results %>%
  group_by(data_type) %>%
  summarise(
    total_tests = n(),
    significant_p05 = sum(pearson_pvalue < 0.05, na.rm = TRUE),
    significant_fdr05 = sum(pearson_fdr < 0.05, na.rm = TRUE),
    mean_correlation = mean(abs(pearson_correlation), na.rm = TRUE),
    .groups = 'drop'
  )

# 按数据类型统计（表型水平）
summary_phenotype_by_type <- phenotype_mean_correlation %>%
  group_by(data_type) %>%
  summarise(
    total_phenotypes = n(),
    significant_p05 = sum(pearson_ttest_pvalue < 0.05, na.rm = TRUE),
    significant_fdr05 = sum(pearson_ttest_fdr < 0.05, na.rm = TRUE),
    mean_abs_correlation = mean(abs(mean_pearson_correlation), na.rm = TRUE),
    mean_subjects_per_phenotype = mean(n_subjects, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n按数据类型统计（个体水平）：\n")
print(summary_by_type)

cat("\n按数据类型统计（表型水平）：\n")
print(summary_phenotype_by_type)

# 保存汇总统计
write.csv(summary_by_type, 
          file.path(output_dir, "correlation_summary_by_datatype_individual.csv"), 
          row.names = FALSE)

write.csv(summary_phenotype_by_type, 
          file.path(output_dir, "correlation_summary_by_datatype_phenotype.csv"), 
          row.names = FALSE)

cat("\n结果已保存到:", output_dir, "\n")
cat("\n主要输出文件：\n")
cat("  1. phenotype_altitude_correlation_results.csv - 个体水平完整结果\n")
cat("  2. significant_phenotype_altitude_correlations.csv - 个体水平显著结果\n")
cat("  3. phenotype_mean_correlation_results.csv - 表型水平平均相关性结果\n")
cat("  4. significant_phenotype_mean_correlations.csv - 表型水平显著结果\n")
cat("  5. correlation_summary_by_datatype_*.csv - 按数据类型汇总统计\n")