# ==============================
#  优化版：表型阶段热图（横向、无分组、无边框、自动宽度、美观统一）
#  适用于：生理、眼部、认知、超声等所有表型数据
# ==============================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)


# ——————————————————————————————————
# 01. 读取数据
# ——————————————————————————————————
people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 差异表型（qval < 0.05）
physiology_DE <- read.csv("data/processed/physiology_DE_results.csv", check.names = FALSE) %>% filter(qval < 0.05)
eye_DE        <- read.csv("data/processed/eye_DE_results.csv",        check.names = FALSE) %>% filter(qval < 0.05)
cognitive_DE  <- read.csv("data/processed/cognitive_DE_results.csv",  check.names = FALSE) %>% filter(qval < 0.05)
ultrasound_DE <- read.csv("data/processed/ultrasound_DE_results.csv", check.names = FALSE) %>% filter(qval < 0.05)

# 原始表型数据（只保留差异项）
physiology_abs <- read.csv("data/20250901/01.构建二维数据集/03.整理表型/生理表型.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(physiology_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(
    time == "T1"             ~ "before",
    time %in% c("T9", "T10") ~ "after",
    TRUE                     ~ "climb"
  )) %>% select(-time)

eye_abs <- read.csv("data/20250901/01.构建二维数据集/03.整理表型/眼部数据.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(eye_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(time == "T1" ~ "before", time %in% c("T9","T10") ~ "after", TRUE ~ "climb")) %>% select(-time)

cognitive_abs <- read.csv("data/20250901/01.构建二维数据集/03.整理表型/认知数据.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(cognitive_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(time == "T1" ~ "before", time %in% c("T9","T10") ~ "after", TRUE ~ "climb")) %>% select(-time)

ultrasound_abs <- read.csv("data/20250901/01.构建二维数据集/03.整理表型/超声数据.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(ultrasound_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(time == "T1" ~ "before", time %in% c("T9","T10") ~ "after", TRUE ~ "climb")) %>% select(-time)

# ——————————————————————————————————
# 02. 核心函数：统一风格横向表型热图
# ——————————————————————————————————
process_phenotype_heatmap <- function(data, title_prefix, output_file, de_data) {
  
  # 1. 长格式 → 计算阶段均值
  long_data <- data %>%
    pivot_longer(cols = -c(subject_id, stage), names_to = "Phenotype", values_to = "value")
  
  mean_data <- long_data %>%
    group_by(stage, Phenotype) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = mean_value, values_fill = NA)
  
  required_stages <- c("before", "climb", "after")
  for (s in required_stages) {
    if (!s %in% colnames(mean_data)) mean_data[[s]] <- NA_real_
  }
  mean_data <- mean_data %>% select(Phenotype, all_of(required_stages))
  
  # 2. 矩阵
  mat <- as.matrix(mean_data[, required_stages])
  rownames(mat) <- mean_data$Phenotype
  
  # 3. 行标准化（0-1）
  scaled_mat <- t(apply(mat, 1, function(x) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    minv <- min(x, na.rm = TRUE)
    maxv <- max(x, na.rm = TRUE)
    if (maxv == minv) return(rep(0.5, length(x)))
    (x - minv) / (maxv - minv)
  }))
  
  # 4. 显著性星号
  sig_df <- de_data %>%
    select(original_name, qval) %>%
    group_by(original_name) %>%
    summarise(qval = min(qval, na.rm = TRUE), .groups = "drop") %>%
    mutate(sig = case_when(
      qval < 0.0001 ~ "****",
      qval < 0.001  ~ "***",
      qval < 0.01   ~ "**",
      qval < 0.05   ~ "*",
      TRUE          ~ ""
    )) %>%
    distinct(original_name, .keep_all = TRUE)
  
  sig_vec <- sig_df$sig[match(rownames(scaled_mat), sig_df$original_name)]
  sig_vec[is.na(sig_vec)] <- ""
  new_rownames <- ifelse(sig_vec == "", rownames(scaled_mat), paste0(rownames(scaled_mat), " ", sig_vec))
  rownames(scaled_mat) <- new_rownames
  
  # 5. 整体行聚类（让相似的表型靠在一起）
  row_dist <- dist(scaled_mat, method = "euclidean")
  row_clust <- hclust(row_dist, method = "ward.D2")
  row_order <- row_clust$order
  scaled_mat <- scaled_mat[row_order, ]
  
  # 6. 转置为横向热图
  scaled_mat_t <- t(scaled_mat)
  
  # 7. 配色
  col_palette <- colorRampPalette(c("#4281A4", "#E4DFDA", "#C1666B"))(1000)
  
  # 8. 自动计算PDF宽度
  n_traits <- ncol(scaled_mat_t)
  pdf_width <- max(6, n_traits * 0.13 + 4.5)   # 与微生物、代谢物完全一致的策略
  
  # 9. 绘图（极简现代风格）
  p <- pheatmap(scaled_mat_t,
                color            = col_palette,
                cluster_rows     = FALSE,       # 只有3行，不聚类
                cluster_cols     = FALSE,       # 已手动排序
                show_rownames    = TRUE,
                show_colnames    = TRUE,
                border_color     = '#F1FAEE',          # 无边框，干净！
                cellwidth        = 5.5,
                cellheight       = 20,
                fontsize         = 9,
                fontsize_row     = 12,
                fontsize_col     = 7.5,
                main             = paste0(title_prefix, "\n(mean value across stages)"),
                display_numbers  = FALSE,       # 不显示数字
                na_col           = "white",
                legend           = TRUE,
                angle_col        = 90,
                gaps_col         = NULL)
  
  # 10. 保存
  pdf(output_file, width = pdf_width, height = 5.5)
  print(p)
  dev.off()
  
  message("√ 已保存: ", output_file)
  message("   表型数量: ", n_traits, "   PDF宽度: ", round(pdf_width, 2), " inches")
}

# ——————————————————————————————————
# 03. 执行绘图（全部统一风格）
# ——————————————————————————————————
output_dir <- "/public/home/xiaokechen/01.HuaDa_Qomolangma/data/20250901/05.差异分析/03.表型差异分析/03.表型差异分析热图/01.stage_heatmap"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_phenotype_heatmap(
  data         = physiology_abs,
  title_prefix = "Physiological Phenotypes",
  output_file  = file.path(output_dir, "Physiology_stage_heatmap_horizontal.pdf"),
  de_data      = physiology_DE
)

process_phenotype_heatmap(
  data         = eye_abs,
  title_prefix = "Eye-related Phenotypes",
  output_file  = file.path(output_dir, "Eye_stage_heatmap_horizontal.pdf"),
  de_data      = eye_DE
)

process_phenotype_heatmap(
  data         = cognitive_abs,
  title_prefix = "Cognitive Phenotypes",
  output_file  = file.path(output_dir, "Cognitive_stage_heatmap_horizontal.pdf"),
  de_data      = cognitive_DE
)

process_phenotype_heatmap(
  data         = ultrasound_abs,
  title_prefix = "Ultrasound Phenotypes",
  output_file  = file.path(output_dir, "Ultrasound_stage_heatmap_horizontal.pdf"),
  de_data      = ultrasound_DE
)

message("所有表型阶段热图绘制完成！风格与微生物、代谢物热图完全统一")