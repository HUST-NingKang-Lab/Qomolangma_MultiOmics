library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)

selected_metabolite <- "4-Guanidinobutanoic acid"

oral_abs <- read.csv("data/processed/oral_species.csv", check.names = FALSE, row.names = 1)
gut_abs <- read.csv("data/processed/gut_species.csv", check.names = FALSE, row.names = 1)
skin_abs <- read.csv("data/processed/skin_species.csv", check.names = FALSE, row.names = 1)

oral_abs <- oral_abs %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time),
         subject_id = gsub("^S0([0-9])$", "S\\1", subject_id)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

gut_abs <- gut_abs %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time),
         subject_id = gsub("^S0([0-9])$", "S\\1", subject_id)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

skin_abs <- skin_abs %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time),
         subject_id = gsub("^S0([0-9])$", "S\\1", subject_id)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

met_data_path <- "/public/home/xiaokechen/01.HuaDa_Qomolangma/raw_data/代谢/Metabolin_plasma.csv"
met_data <- read.csv(met_data_path, check.names = FALSE)
met_data <- met_data %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time),
         subject_id = gsub("^S0([0-9])$", "S\\1", subject_id)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

# 函数：对菌群数据按subject_id和time求均值
aggregate_microbiome <- function(df) {
  df_aggregated <- df %>%
    group_by(subject_id, time) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")
  
  return(df_aggregated)
}

# 对三个菌群数据框进行聚合
oral_abs_agg <- aggregate_microbiome(oral_abs)
gut_abs_agg <- aggregate_microbiome(gut_abs)
skin_abs_agg <- aggregate_microbiome(skin_abs)

# 提取目标代谢物数据
if (!selected_metabolite %in% colnames(met_data)) {
  stop(paste("代谢物", selected_metabolite, "不存在于met_data中"))
}

met_selected <- met_data %>%
  select(subject_id, time, all_of(selected_metabolite))

# 定义年龄组：S1和S10为老年人，其他为年轻人
elderly_ids <- c("S1", "S10")
young_ids <- paste0("S", c(2:9, 11:13))

# 函数：计算Spearman相关性（支持按年龄组过滤）
calculate_spearman <- function(microbiome_df, metabolite_df, source_name, age_group = "all") {
  # 根据年龄组过滤数据
  if (age_group == "elderly") {
    microbiome_df <- microbiome_df %>% filter(subject_id %in% elderly_ids)
    metabolite_df <- metabolite_df %>% filter(subject_id %in% elderly_ids)
  } else if (age_group == "young") {
    microbiome_df <- microbiome_df %>% filter(subject_id %in% young_ids)
    metabolite_df <- metabolite_df %>% filter(subject_id %in% young_ids)
  }
  
  # 合并菌群数据和代谢数据
  merged_data <- inner_join(microbiome_df, metabolite_df, by = c("subject_id", "time"))
  
  # 获取菌群列名（排除subject_id和time）
  bacteria_cols <- setdiff(colnames(microbiome_df), c("subject_id", "time"))
  
  # 存储结果
  results <- data.frame()
  
  # 对每个菌进行Spearman相关性分析
  for (bacteria in bacteria_cols) {
    # 提取有效数据（去除NA值）
    valid_data <- merged_data %>%
      select(all_of(c(bacteria, selected_metabolite))) %>%
      na.omit()
    
    # 如果有效数据点少于3个，跳过
    if (nrow(valid_data) < 3) {
      next
    }
    
    # 计算Spearman相关性
    cor_test <- cor.test(valid_data[[bacteria]], 
                         valid_data[[selected_metabolite]], 
                         method = "spearman", 
                         exact = FALSE)
    
    # 保存结果
    results <- rbind(results, data.frame(
      Metabolite = selected_metabolite,
      Source = source_name,
      Age_group = age_group,
      Bacteria = bacteria,
      Correlation = cor_test$estimate,
      P_value = cor_test$p.value,
      Significant = cor_test$p.value < 0.05,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# 对三个菌群数据分别进行相关性分析（老年人、年轻人、全体）
# 老年人
oral_elderly <- calculate_spearman(oral_abs_agg, met_selected, "oral", "elderly")
gut_elderly <- calculate_spearman(gut_abs_agg, met_selected, "gut", "elderly")
skin_elderly <- calculate_spearman(skin_abs_agg, met_selected, "skin", "elderly")

# 年轻人
oral_young <- calculate_spearman(oral_abs_agg, met_selected, "oral", "young")
gut_young <- calculate_spearman(gut_abs_agg, met_selected, "gut", "young")
skin_young <- calculate_spearman(skin_abs_agg, met_selected, "skin", "young")

# 全体
oral_all <- calculate_spearman(oral_abs_agg, met_selected, "oral", "all")
gut_all <- calculate_spearman(gut_abs_agg, met_selected, "gut", "all")
skin_all <- calculate_spearman(skin_abs_agg, met_selected, "skin", "all")

# 合并所有结果
final_results <- rbind(
  oral_elderly, gut_elderly, skin_elderly,
  oral_young, gut_young, skin_young,
  oral_all, gut_all, skin_all
)

# 按年龄组、source、p值排序
final_results <- final_results %>%
  arrange(Age_group, Source, P_value)

# 查看结果
print(head(final_results, 30))

# 保存结果
write.csv(final_results, 
          "data/processed/spearman_correlation_results.csv", 
          row.names = FALSE)

# 统计显著相关的菌
significant_summary <- final_results %>%
  group_by(Age_group, Source) %>%
  summarise(
    Total = n(),
    Significant = sum(Significant),
    Percentage = round(Significant/Total * 100, 2),
    .groups = "drop"
  ) %>%
  arrange(Age_group, Source)

print("=== 显著相关菌数统计 ===")
print(significant_summary)

# 保存统计表
write.csv(significant_summary, 
          "data/processed/significant_summary.csv", 
          row.names = FALSE)


## 读入差异菌
## 读入bin和物种信息的对应关系
MAG_meta <- read.csv("data/processed/02.MAG_stat_with_metadata_R226.csv", check.names = FALSE) %>%
  select(MAG, classification) %>%
  mutate(
    # 提取 phylum 信息 (第二个层级)
    phylum = sub(".*;p__([^;]+);.*", "\\1", classification),
    # 提取 genus 信息 (第六个层级)
    genus = sub(".*;g__([^;]+);.*", "\\1", classification),
    # 保留 genus 和 species 用于显示
    classification_short = sub(".*;([^;]*;[^;]*)$", "\\1", classification)
  ) %>%
  mutate(MAG = gsub(".fa", "", MAG))

# 读取差异分析结果并过滤
gut_DE <- read.csv("data/processed/gut_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(GTDB_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

oral_DE <- read.csv("data/processed/oral_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(GTDB_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

skin_DE <- read.csv("data/processed/skin_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(GTDB_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

# 提取各差异分析结果中的 original_name 作为集合，便于快速查找
gut_DE_names   <- gut_DE$original_name
oral_DE_names  <- oral_DE$original_name
skin_DE_names  <- skin_DE$original_name

# 初始化三列为 FALSE
final_results$gut_DE  <- FALSE
final_results$oral_DE <- FALSE
final_results$skin_DE <- FALSE

# 根据 Source 分别赋值
final_results <- final_results %>%
  mutate(
    gut_DE  = ifelse(Source == "gut"  & Bacteria %in% gut_DE_names,  TRUE, gut_DE),
    oral_DE = ifelse(Source == "oral" & Bacteria %in% oral_DE_names, TRUE, oral_DE),
    skin_DE = ifelse(Source == "skin" & Bacteria %in% skin_DE_names, TRUE, skin_DE)
  )

final_results <- merge(MAG_meta, final_results, 
                       by.x = "MAG", 
                       by.y = "Bacteria", 
                       all.y = TRUE)

# 可选：保存更新后的结果
write.csv(final_results, 
          "data/processed/spearman_correlation_results_with_DE.csv", 
          row.names = FALSE)