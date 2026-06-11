# 加载必要的包
library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# 设置工作目录


# —————————————————————————————————————————————
# 01. 读取并处理数据
# —————————————————————————————————————————————
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

people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 读取差异分析结果并过滤
gut_DE <- read.csv("data/processed/gut_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(original_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

oral_DE <- read.csv("data/processed/oral_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(original_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

skin_DE <- read.csv("data/processed/skin_DE_results.csv", check.names = FALSE) %>%
  filter(qval < 0.05) %>%
  left_join(MAG_meta, by = c("original_name" = "MAG")) %>%
  mutate(original_name = ifelse(is.na(classification_short), original_name, classification_short)) %>%
  filter(!is.na(original_name))

# 读取菌群数据并过滤
gut_abs <- read.csv("data/processed/gut_species.csv", check.names = FALSE)
colnames(gut_abs) <- ifelse(colnames(gut_abs) %in% MAG_meta$MAG,
                            MAG_meta$classification_short[match(colnames(gut_abs), MAG_meta$MAG)],
                            colnames(gut_abs))
gut_abs <- gut_abs %>%
  select(subject_id, time, any_of(unique(gut_DE$original_name))) %>%
  filter(subject_id %in% people & time %in% time12) %>% 
  mutate(time = case_when(
    time == "T1" ~ "before", 
    time %in% c("T9", "T10") ~ "after", 
    TRUE ~ "climb"
  ))

oral_abs <- read.csv("data/2processed/oral_species.csv", check.names = FALSE)
colnames(oral_abs) <- ifelse(colnames(oral_abs) %in% MAG_meta$MAG,
                             MAG_meta$classification_short[match(colnames(oral_abs), MAG_meta$MAG)],
                             colnames(oral_abs))
oral_abs <- oral_abs %>%
  select(subject_id, time, any_of(unique(oral_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12) %>% 
  mutate(time = case_when(
    time == "T1" ~ "before", 
    time %in% c("T9", "T10") ~ "after", 
    TRUE ~ "climb"
  ))

skin_abs <- read.csv("data/processed/skin_species.csv", check.names = FALSE)
colnames(skin_abs) <- ifelse(colnames(skin_abs) %in% MAG_meta$MAG,
                             MAG_meta$classification_short[match(colnames(skin_abs), MAG_meta$MAG)],
                             colnames(skin_abs))
skin_abs <- skin_abs %>%
  select(subject_id, time, any_of(unique(skin_DE$original_name))) %>% 
  filter(subject_id %in% people & time %in% time12) %>% 
  mutate(time = case_when(
    time == "T1" ~ "before", 
    time %in% c("T9", "T10") ~ "after", 
    TRUE ~ "climb"
  ))

# —————————————————————————————————————————————
# 02. 定义函数：生成按 phylum 分组的横向热图（无数字、更窄、PDF不截断）
# —————————————————————————————————————————————
process_and_plot_heatmap <- function(data, title, output_file, de_data, mag_meta) {
  # 1. 长格式 + 计算每个时间段的均值
  long_data <- data %>%
    pivot_longer(cols = -c(subject_id, time), names_to = "Species", values_to = "value")
  
  mean_data <- long_data %>%
    group_by(time, Species) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = time, values_from = mean_value, values_fill = NA)
  
  required_times <- c("before", "climb", "after")
  for (t in required_times) {
    if (!t %in% colnames(mean_data)) mean_data[[t]] <- NA_real_
  }
  mean_data <- mean_data %>% select(Species, all_of(required_times))
  
  # 原始矩阵
  raw_matrix <- as.matrix(mean_data[, required_times])
  rownames(raw_matrix) <- mean_data$Species
  
  # 标准化矩阵
  scaled_matrix <- t(apply(raw_matrix, 1, function(x) {
    if (all(is.na(x))) {
      rep(NA, length(x))
    } else {
      min_val <- min(x, na.rm = TRUE)
      max_val <- max(x, na.rm = TRUE)
      if (max_val == min_val) {
        rep(0.5, length(x))
      } else {
        (x - min_val) / (max_val - min_val)
      }
    }
  }))
  
  # 显著性标记
  sig_info <- de_data %>%
    select(original_name, qval) %>%
    group_by(original_name) %>%
    summarise(qval = min(qval, na.rm = TRUE), .groups = "drop") %>%
    mutate(significance = case_when(
      qval < 0.0001 ~ "****",
      qval < 0.001  ~ "***",
      qval < 0.01   ~ "**",
      qval < 0.05   ~ "*",
      TRUE          ~ ""
    )) %>%
    column_to_rownames("original_name")
  
  sig_info <- sig_info[rownames(scaled_matrix), , drop = FALSE]
  new_rownames <- paste0(rownames(scaled_matrix),
                         ifelse(sig_info$significance == "", "", " "),
                         sig_info$significance)
  rownames(scaled_matrix) <- new_rownames
  
  # ==================== 添加 phylum 注释 ====================
  annotation_df <- data.frame(
    Species = rownames(scaled_matrix),
    stringsAsFactors = FALSE
  ) %>%
    mutate(Species_clean = gsub(" \\*+$", "", Species)) %>%
    left_join(
      mag_meta %>% 
        select(classification_short, phylum) %>% 
        distinct(),
      by = c("Species_clean" = "classification_short")
    ) %>%
    select(Phylum = phylum) %>%
    as.data.frame()
  
  rownames(annotation_df) <- rownames(scaled_matrix)
  
  unique_phyla <- unique(annotation_df$Phylum)
  unique_phyla <- unique_phyla[!is.na(unique_phyla)]
  
  if(length(unique_phyla) <= 8) {
    phylum_colors <- brewer.pal(max(3, length(unique_phyla)), "Set2")
  } else if(length(unique_phyla) <= 12) {
    phylum_colors <- brewer.pal(length(unique_phyla), "Set3")
  } else {
    phylum_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_phyla))
  }
  names(phylum_colors) <- unique_phyla
  
  annotation_colors <- list(Phylum = phylum_colors)
  
  # ==================== 按 phylum 分组并在每个 phylum 内聚类 ====================
  annotation_df_sorted <- annotation_df %>%
    mutate(row_id = rownames(.))
  
  reordered_rows <- c()
  gap_positions <- c()
  
  for(phl in unique_phyla) {
    if(is.na(phl)) next
    
    phylum_rows <- annotation_df_sorted %>% filter(Phylum == phl)
    phylum_row_ids <- phylum_rows$row_id
    
    if(length(phylum_row_ids) == 1) {
      reordered_rows <- c(reordered_rows, phylum_row_ids)
    } else {
      phylum_matrix <- scaled_matrix[phylum_row_ids, , drop = FALSE]
      if(nrow(phylum_matrix) > 1) {
        row_dist <- dist(phylum_matrix)
        row_clust <- hclust(row_dist)
        clustered_order <- phylum_row_ids[row_clust$order]
        reordered_rows <- c(reordered_rows, clustered_order)
      } else {
        reordered_rows <- c(reordered_rows, phylum_row_ids)
      }
    }
    
    if(length(reordered_rows) > 0) {
      gap_positions <- c(gap_positions, length(reordered_rows))
    }
  }
  
  if(length(gap_positions) > 0) {
    gap_positions <- gap_positions[-length(gap_positions)]
  }
  
  scaled_matrix <- scaled_matrix[reordered_rows, ]
  annotation_df <- annotation_df[reordered_rows, , drop = FALSE]
  
  # ==================== 转置为横向热图 ====================
  scaled_matrix_t <- t(scaled_matrix)
  annotation_col <- annotation_df
  
  color_palette <- colorRampPalette(c("#4281A4", "#E4DFDA", "#C1666B"))(1000)
  
  # ==================== 绘图（关键修改：border_color = NA）====================
  p <- pheatmap(scaled_matrix_t,
                color = color_palette,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                border_color = '#F1FAEE',       # ⭐⭐ 关键：去除单元格边框！
                cellwidth = 4.5,           # 稍微放宽一点，避免过窄
                cellheight = 8,
                fontsize = 8,
                fontsize_row = 8,
                fontsize_col = 6,
                main = paste0(title, "\n(mean abundance, grouped by Phylum)"),
                display_numbers = FALSE,
                angle_col = 90,
                na_col = "#d9d9d9",
                legend = TRUE,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                annotation_legend = TRUE,
                annotation_names_col = TRUE,
                fontsize_title = 12,
                gaps_col = gap_positions)
  
  # 调整PDF尺寸，确保不截断
  fig_width <- max(6, ncol(scaled_matrix_t) * 0.13 + 4.5)
  pdf(output_file, width = fig_width, height = 5)
  print(p)
  dev.off()
  
  message("Heatmap saved: ", output_file)
  message("   Total species: ", ncol(scaled_matrix_t))
}

# —————————————————————————————————————————————
# 03. 输出结果
# —————————————————————————————————————————————
output_dir <- getwd()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_and_plot_heatmap(gut_abs, "Gut Microbe",
                         file.path(output_dir, "gut_heatmap_by_phylum_horizontal.pdf"),
                         gut_DE, MAG_meta)

process_and_plot_heatmap(oral_abs, "Oral Microbe",
                         file.path(output_dir, "oral_heatmap_by_phylum_horizontal.pdf"),
                         oral_DE, MAG_meta)

process_and_plot_heatmap(skin_abs, "Skin Microbe",
                         file.path(output_dir, "skin_heatmap_by_phylum_horizontal.pdf"),
                         skin_DE, MAG_meta)

# 检查警告
warnings()