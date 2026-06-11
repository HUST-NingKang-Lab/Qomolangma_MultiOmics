setwd("/public/home/xiaokechen/01.HuaDa_Qomolangma/Catalog")
library(dplyr)
library(stringr)
library(tidyr)
# ———————————————————————————————————
## step1: Phylum注释——01_01.Phylum.csv
# ———————————————————————————————————
data <- read.csv("data/processed/02.MAG_stat_with_metadata_R226.csv", check.names = FALSE)
sc_data <- data %>% 
    select(refMAG, classification) %>%
    mutate(Phylum = str_split_i(classification, ";", 2)) %>%
    mutate(Phylum = str_remove(Phylum, "^p__")) %>% 
    unique() %>% 
    select(-classification)

phylum_counts <- sc_data %>% count(Phylum, sort = TRUE)
top_5_phyla <- phylum_counts$Phylum[1:5]

# 来源：常用科研配色，如 ColorBrewer Set1, 或发表文章常用色
scientific_colors <- c("#84C3B7", "#7DA6C6", "#B7B2D0", "#E68B81", "#EAAA60")
color_map <- setNames(scientific_colors, top_5_phyla)

sc_data <- sc_data %>%
    mutate(Color = ifelse(Phylum %in% top_5_phyla, 
                          color_map[Phylum], 
                          "#999999"), # 灰色用于其他
           Color = as.character(Color)) # 确保不是因子

write.csv(sc_data, "data/02.MAG_analysis/03.绘制树图/01_01.Phylum.csv", row.names = FALSE, quote = FALSE)

# ———————————————————————————————————
## step2: MAG在各部位检出数量——01_2.coverage.csv
# ———————————————————————————————————
data <- read.csv("data/02.MAG_analysis/01.整理MAG信息/06.统计每个部位下每个物种被检测的次数.csv", check.names = FALSE)

library(dplyr)
library(tidyr)

# 读取数据
data <- read.csv("data/02.MAG_analysis/01.整理MAG信息/06.统计每个部位下每个物种被检测的次数.csv", check.names = FALSE)
data <- data %>%
  pivot_wider(
    names_from = Source,        # 将 Source 的唯一值作为新列名
    values_from = MAG_Count,    # 使用 MAG_Count 的值填充
    values_fill = 0             # 缺失值填充为 0
  ) %>%
  select(refMAG, gut, skin, oral)

# 对第2到第4列进行 Min-Max 归一化
data <- data %>% mutate(across(2:4, ~ (.- min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))

write.csv(data, "data/02.MAG_analysis/03.绘制树图/01_02.coverage.csv", row.names = FALSE, quote = FALSE)

# ———————————————————————————————————
## step3: MAG丰度最低或最高的时间点——01_3.abs_MinT_MaxT.csv
# ———————————————————————————————————
abs <- read.csv("00.MAG_Catalog/03.MAG_abs/02.MAG_abs_composition_normalized.csv", check.names = FALSE)
meta <- read.csv("raw_data/metadata_20250819.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time12")]

# Reshape abs to long format
abs_long <- abs %>%
  pivot_longer(cols = -Genome, names_to = "SampleID", values_to = "Abundance")

data <- abs_long %>% left_join(meta, by = "SampleID")

avg_data <- data %>%
  group_by(Genome, Source, time12) %>%
  summarise(Avg_Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

result <- avg_data %>%
  group_by(Genome, Source) %>%
  summarise(
    MinT = if (all(Avg_Abundance == 0, na.rm = TRUE)) NA else time12[which.min(Avg_Abundance)],
    MaxT = if (all(Avg_Abundance == 0, na.rm = TRUE)) NA else time12[which.max(Avg_Abundance)],
    .groups = "drop"
  ) %>%
  # Pivot to get gut, skin, oral columns
  pivot_wider(
    names_from = Source,
    values_from = c(MinT, MaxT),
    names_glue = "{Source}_{.value}"
  ) %>%
  # Ensure all expected columns exist
  complete(Genome, fill = list(
    gut_MinT = NA, gut_MaxT = NA,
    skin_MinT = NA, skin_MaxT = NA,
    oral_MinT = NA, oral_MaxT = NA
  ))

write.csv(result, "data/02.MAG_analysis/03.绘制树图/01_03.abs_MinT_MaxT.csv", row.names = FALSE, quote = FALSE)

# ———————————————————————————————————
## step4: MAG被检出最多或最少的时间点——01_04.coverage_MinT_MaxT.csv
# ———————————————————————————————————
data <- read.csv("data/02.MAG_analysis/01.整理MAG信息/03.每个样本恢复的refMAG覆盖情况[发现每个样本统计的MAG都是唯一物种].csv", check.names = FALSE) %>% 
  select(SampleID, refMAG, Source, time12)
meta <- read.csv("raw_data/metadata_20250819.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time12")]

# 计算meta表中每个Source和time12的样本数量
sample_counts <- meta %>%
  group_by(Source, time12) %>%
  summarise(Sample_Count = n(), .groups = "drop")

View(sample_counts)

# 计算每个refMAG在每个Source和time12的流行率（行数/样本数量）
prevalence_data <- data %>%
  group_by(refMAG, Source, time12) %>%
  summarise(Row_Count = n(), .groups = "drop") %>%
  left_join(sample_counts, by = c("Source", "time12")) %>%
  mutate(Prevalence = Row_Count / Sample_Count)

# 为每个refMAG和Source，找到流行率最低和最高的time12
# 如果某个Source没有行或所有流行率为0，返回NA
# 如果MinT和MaxT相同，设置为"Only"
result <- prevalence_data %>%
  group_by(refMAG, Source) %>%
  summarise( 
    MinT = if (all(Prevalence == 0, na.rm = TRUE) || n() == 0) NA else time12[which.min(Prevalence)],
    MaxT = if (all(Prevalence == 0, na.rm = TRUE) || n() == 0) NA else time12[which.max(Prevalence)],
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Source,
    values_from = c(MinT, MaxT),
    names_glue = "{Source}_{.value}"
  ) %>%
  # 确保所有预期列存在
  complete(refMAG, fill = list(
    gut_MinT = NA, gut_MaxT = NA,
    skin_MinT = NA, skin_MaxT = NA,
    oral_MinT = NA, oral_MaxT = NA
  )) %>%
  rename(Genome = refMAG) %>%
  # 检查每个Source的MinT和MaxT是否相同，若相同则设置为"Only"
  mutate(
    gut_MinT = if_else(gut_MinT == gut_MaxT & !is.na(gut_MinT), "Only", as.character(gut_MinT)),
    gut_MaxT = if_else(gut_MinT == "Only", "Only", as.character(gut_MaxT)),
    skin_MinT = if_else(skin_MinT == skin_MaxT & !is.na(skin_MinT), "Only", as.character(skin_MinT)),
    skin_MaxT = if_else(skin_MinT == "Only", "Only", as.character(skin_MaxT)),
    oral_MinT = if_else(oral_MinT == oral_MaxT & !is.na(oral_MinT), "Only", as.character(oral_MinT)),
    oral_MaxT = if_else(oral_MinT == "Only", "Only", as.character(oral_MaxT))
  )

write.csv(result, "data/02.MAG_analysis/03.绘制树图/01_04.coverage_MinT_MaxT.csv", row.names = FALSE, quote = FALSE)

# ———————————————————————————————————
## step5: MMAG被检出只存在于upland和lowland的情况——01_05.Only_Upland_Lowland.csv
# ———————————————————————————————————
data <- read.csv("data/02.MAG_analysis/01.整理MAG信息/03.每个样本恢复的refMAG覆盖情况[发现每个样本统计的MAG都是唯一物种].csv", check.names = FALSE) %>% 
  select(SampleID, refMAG, Source, time12)
meta <- read.csv("raw_data/metadata_20250819.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time12")]

# 定义lowland和upland的时间点
lowland_times <- c("T1", "T9", "T10")
upland_times <- c("T2", "T3", "T4", "T5", "T6", "T7", "T8", "Ta", "Tb")

# 获取每个refMAG在每个Source中存在的时间点
habitat_data <- data %>%
  distinct(refMAG, Source, time12) %>%
  group_by(refMAG, Source) %>%
  summarise(
    # 检查是否在lowland时间点存在
    in_lowland = any(time12 %in% lowland_times),
    # 检查是否在upland时间点存在
    in_upland = any(time12 %in% upland_times),
    .groups = "drop"
  ) %>%
  mutate(
    # 判定栖息地特异性
    Habitat = case_when(
      in_lowland & in_upland ~ "all",  # 同时存在于两种环境
      in_lowland & !in_upland ~ "lowland",  # 只在lowland
      !in_lowland & in_upland ~ "upland",  # 只在upland
      TRUE ~ "NA"  # 其他情况（理论上不应该出现）
    )
  ) %>%
  select(refMAG, Source, Habitat) %>%
  pivot_wider(
    names_from = Source,
    values_from = Habitat,
    names_glue = "{Source}_Habitat"
  ) %>%
  # 确保所有预期列存在
  complete(refMAG, fill = list(
    gut_Habitat = NA,
    skin_Habitat = NA,
    oral_Habitat = NA
  )) %>%
  rename(Genome = refMAG)

write.csv(habitat_data, "data/02.MAG_analysis/03.绘制树图/01_05.Only_Upland_Lowland.csv", row.names = FALSE, quote = FALSE)

# 查看每个Source中不同栖息地特异性的MAG数量统计
summary_stats <- habitat_data %>%
  pivot_longer(
    cols = ends_with("_Habitat"),
    names_to = "Source",
    values_to = "Habitat"
  ) %>%
  mutate(Source = gsub("_Habitat", "", Source)) %>%
  group_by(Source, Habitat) %>%
  summarise(Count = n(), .groups = "drop")

print("栖息地特异性统计：")
print(summary_stats)

# [1] "栖息地特异性统计："
# # A tibble: 12 × 3
#    Source Habitat Count
#    <chr>  <chr>   <int>
#  1 gut    all       243
#  2 gut    lowland    53
#  3 gut    upland    308
#  4 gut    NA       1469
#  5 oral   all       339
#  6 oral   lowland   127
#  7 oral   upland    693
#  8 oral   NA        914
#  9 skin   all        42
# 10 skin   lowland    74
# 11 skin   upland    360
# 12 skin   NA       1597