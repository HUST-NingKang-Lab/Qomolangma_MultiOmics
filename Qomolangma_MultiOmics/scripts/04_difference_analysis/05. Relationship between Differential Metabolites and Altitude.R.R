
library(dplyr)
library(tidyr)

# —————————————————————————————————————————————
# 01. 读取并处理数据
# —————————————————————————————————————————————
people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 读取差异分析结果并过滤
Lipid_DE <- read.csv("data/processed/Lipid_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)
Metabolin_DE <- read.csv("data/processed/Metabolin_DE_results.csv", check.names = FALSE) %>% 
  filter(qval < 0.05)

# 读取代谢数据并过滤
Lipid_abs <- read.csv("ddata/processed/Lipid_log10.csv", check.names = FALSE) %>%
  select(subject_id, time, all_of(unique(Lipid_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12)
Metabolin_abs <- read.csv("data/processed/Metabolin_log10.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(Metabolin_DE$original_name))) %>% 
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
  # 获取代谢列名（排除subject_id和time列）
  metabolite_cols <- setdiff(colnames(data), c("subject_id", "time"))
  
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
    
    # 对每个代谢计算相关性
    for (metabolite in metabolite_cols) {
      metabolite_data <- person_data_with_altitude %>%
        select(time, altitude, !!sym(metabolite)) %>%
        filter(!is.na(!!sym(metabolite)))  # 移除代谢缺失值
      
      # 如果数据点不足，跳过
      if (nrow(metabolite_data) < 3) next
      
      # 计算Pearson相关性（趋势相关性）
      cor_result <- cor.test(metabolite_data$altitude, metabolite_data[[metabolite]], 
                           method = "pearson")
      
      # 计算Spearman相关性（作为补充，更稳健）
      spearman_result <- cor.test(metabolite_data$altitude, metabolite_data[[metabolite]], 
                                method = "spearman")
      
      # 存储结果
      results_list[[length(results_list) + 1]] <- data.frame(
        data_type = data_type,
        subject_id = person,
        metabolite = metabolite,
        n_points = nrow(metabolite_data),
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
cat("正在分析脂质代谢...\n")
Lipid_results <- calculate_correlation(Lipid_abs, "Lipid")

cat("正在分析血液代谢...\n")
Metabolin_results <- calculate_correlation(Metabolin_abs, "eye")

# —————————————————————————————————————————————
# 04. 合并所有结果
# —————————————————————————————————————————————
all_results <- rbind(
  Lipid_results,
  Metabolin_results
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
output_dir <- "data/processed"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存完整结果
write.csv(all_results, 
          file.path(output_dir, "metabolite_altitude_correlation_results.csv"), 
          row.names = FALSE)

# 保存显著结果（Pearson相关性p<0.05）
significant_results <- all_results %>%
  filter(pearson_pvalue < 0.05)

write.csv(significant_results, 
          file.path(output_dir, "significant_metabolite_altitude_correlations.csv"), 
          row.names = FALSE)

# —————————————————————————————————————————————
# 06. 计算每个代谢在所有个体中的平均相关性
# —————————————————————————————————————————————
cat("正在计算代谢平均相关性...\n")

# 按代谢和数据类型分组，计算平均相关性
metabolite_mean_correlation <- all_results %>%
  group_by(data_type, metabolite) %>%
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
metabolite_mean_correlation$pearson_ttest_pvalue <- NA
metabolite_mean_correlation$spearman_ttest_pvalue <- NA

for (i in 1:nrow(metabolite_mean_correlation)) {
  current_data_type <- metabolite_mean_correlation$data_type[i]
  current_metabolite <- metabolite_mean_correlation$metabolite[i]
  
  # 获取该代谢的所有个体相关性数据
  metabolite_correlations <- all_results %>%
    filter(data_type == current_data_type, metabolite == current_metabolite)
  
  # 对Pearson相关性进行单样本t检验
  if (nrow(metabolite_correlations) > 1) {
    pearson_ttest <- t.test(metabolite_correlations$pearson_correlation, mu = 0)
    metabolite_mean_correlation$pearson_ttest_pvalue[i] <- pearson_ttest$p.value
    
    spearman_ttest <- t.test(metabolite_correlations$spearman_correlation, mu = 0)
    metabolite_mean_correlation$spearman_ttest_pvalue[i] <- spearman_ttest$p.value
  }
}

# FDR校正
metabolite_mean_correlation$pearson_ttest_fdr <- p.adjust(metabolite_mean_correlation$pearson_ttest_pvalue, method = "fdr")
metabolite_mean_correlation$spearman_ttest_fdr <- p.adjust(metabolite_mean_correlation$spearman_ttest_pvalue, method = "fdr")

# 按平均相关性的绝对值排序
metabolite_mean_correlation <- metabolite_mean_correlation %>%
  arrange(desc(abs(mean_pearson_correlation)))

# 筛选显著结果（基于t检验p<0.05）
significant_metabolite_mean <- metabolite_mean_correlation %>%
  filter(pearson_ttest_pvalue < 0.05 | spearman_ttest_pvalue < 0.05)

# —————————————————————————————————————————————
# 07. 保存结果
# —————————————————————————————————————————————

# 保存代谢平均相关性结果
write.csv(metabolite_mean_correlation, 
          file.path(output_dir, "metabolite_mean_correlation_results.csv"), 
          row.names = FALSE)

# 保存显著的代谢平均相关性结果
write.csv(significant_metabolite_mean, 
          file.path(output_dir, "significant_metabolite_mean_correlations.csv"), 
          row.names = FALSE)

# —————————————————————————————————————————————
# 08. 结果汇总
# —————————————————————————————————————————————
cat("\n=== 分析完成 ===\n")
cat("个体水平分析：\n")
cat("  总共分析了", nrow(all_results), "个代谢-个体组合\n")
cat("  显著相关（p<0.05）的有", nrow(significant_results), "个\n")
cat("  FDR校正后显著（FDR<0.05）的有", sum(all_results$pearson_fdr < 0.05, na.rm = TRUE), "个\n")

cat("\n代谢水平分析：\n")
cat("  总共分析了", nrow(metabolite_mean_correlation), "个代谢\n")
cat("  平均相关性显著（t检验p<0.05）的有", nrow(significant_metabolite_mean), "个\n")
cat("  FDR校正后显著（FDR<0.05）的有", 
    sum(metabolite_mean_correlation$pearson_ttest_fdr < 0.05, na.rm = TRUE), "个\n")

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

# 按数据类型统计（代谢水平）
summary_metabolite_by_type <- metabolite_mean_correlation %>%
  group_by(data_type) %>%
  summarise(
    total_metabolites = n(),
    significant_p05 = sum(pearson_ttest_pvalue < 0.05, na.rm = TRUE),
    significant_fdr05 = sum(pearson_ttest_fdr < 0.05, na.rm = TRUE),
    mean_abs_correlation = mean(abs(mean_pearson_correlation), na.rm = TRUE),
    mean_subjects_per_metabolite = mean(n_subjects, na.rm = TRUE),
    .groups = 'drop'
  )

cat("\n按数据类型统计（个体水平）：\n")
print(summary_by_type)

cat("\n按数据类型统计（代谢水平）：\n")
print(summary_metabolite_by_type)

# 保存汇总统计
write.csv(summary_by_type, 
          file.path(output_dir, "correlation_summary_by_datatype_individual.csv"), 
          row.names = FALSE)

write.csv(summary_metabolite_by_type, 
          file.path(output_dir, "correlation_summary_by_datatype_metabolite.csv"), 
          row.names = FALSE)

cat("\n结果已保存到:", output_dir, "\n")
cat("\n主要输出文件：\n")
cat("  1. metabolite_altitude_correlation_results.csv - 个体水平完整结果\n")
cat("  2. significant_metabolite_altitude_correlations.csv - 个体水平显著结果\n")
cat("  3. metabolite_mean_correlation_results.csv - 代谢水平平均相关性结果\n")
cat("  4. significant_metabolite_mean_correlations.csv - 代谢水平显著结果\n")
cat("  5. correlation_summary_by_datatype_*.csv - 按数据类型汇总统计\n")