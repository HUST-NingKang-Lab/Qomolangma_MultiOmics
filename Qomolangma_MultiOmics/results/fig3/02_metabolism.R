# ==============================
#  优化版：代谢物阶段热图（横向、无分组、无边框、自动宽度）
#  作者：基于你的原始代码 + 参考微生物热图最佳实践
# ==============================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)



# ——————————————————————————————————
# 01. 读取数据
# ——————————————————————————————————
people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 差异代谢物
Lipid_DE    <- read.csv("data/processed/Lipid_DE_results.csv",    check.names = FALSE) %>% filter(qval < 0.05)
Metabolin_DE<- read.csv("data/processed/Metabolin_DE_results.csv",check.names = FALSE) %>% filter(qval < 0.05)

# 丰度表（已log10）
Lipid_abs <- read.csv("data/processed/Lipid_log10.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(Lipid_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(
    time == "T1"                  ~ "before",
    time %in% c("T9", "T10")      ~ "after",
    TRUE                          ~ "climb"
  )) %>% select(-time)

Metabolin_abs <- read.csv("data/processed/Metabolin_log10.csv", check.names = FALSE) %>%
  select(subject_id, time, any_of(Metabolin_DE$original_name)) %>%
  filter(subject_id %in% people, time %in% time12) %>%
  mutate(stage = case_when(
    time == "T1"                  ~ "before",
    time %in% c("T9", "T10")      ~ "after",
    TRUE                          ~ "climb"
  )) %>% select(-time)

# ——————————————————————————————————
# 02. 核心函数：横向、标准化、无边框、带星号、自动宽度
# ——————————————————————————————————
process_metabolite_heatmap <- function(data, title_prefix, output_file, de_data) {
  
  # 1. 长格式 → 计算每个阶段均值
  long_data <- data %>%
    pivot_longer(cols = -c(subject_id, stage), names_to = "Metabolite", values_to = "value")
  
  mean_data <- long_data %>%
    group_by(stage, Metabolite) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = stage, values_from = mean_value, values_fill = NA)
  
  # 确保三个阶段都在
  required_stages <- c("before", "climb", "after")
  for (s in required_stages) {
    if (!s %in% colnames(mean_data)) mean_data[[s]] <- NA_real_
  }
  mean_data <- mean_data %>% select(Metabolite, all_of(required_stages))
  
  # 2. 构造矩阵
  mat <- as.matrix(mean_data[, required_stages])
  rownames(mat) <- mean_data$Metabolite
  
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
  
  # 匹配到矩阵行
  sig_vec <- sig_df$sig[match(rownames(scaled_mat), sig_df$original_name)]
  sig_vec[is.na(sig_vec)] <- ""
  new_rownames <- ifelse(sig_vec == "", 
                         rownames(scaled_mat),
                         paste0(rownames(scaled_mat), " ", sig_vec))
  rownames(scaled_mat) <- new_rownames
  
  # 5. 整体聚类（可选：这里保留聚类，更美观；如不需要可设 cluster_cols = FALSE）
  #    由于只有3列，不聚类列；行聚类让相似代谢物靠在一起
  row_dist <- dist(scaled_mat, method = "euclidean")
  row_clust <- hclust(row_dist, method = "ward.D2")
  row_order <- row_clust$order
  scaled_mat <- scaled_mat[row_order, ]
  
  # 6. 转置 → 横向热图（代谢物在横轴）
  scaled_mat_t <- t(scaled_mat)
  
  # 7. 颜色
  col_palette <- colorRampPalette(c("#4281A4", "#E4DFDA", "#C1666B"))(1000)
  
  # 8. 计算合适的PDF宽度（每列约0.13英寸 + 固定部分）
  n_metab <- ncol(scaled_mat_t)
  pdf_width <- max(6, n_metab * 0.13 + 4.5)   # 与参考代码完全一致的策略
  
  # 9. 绘图
  p <- pheatmap(scaled_mat_t,
                color            = col_palette,
                cluster_rows     = FALSE,           # 只有3行，不聚类
                cluster_cols     = FALSE,           # 我们已经手动排序过了
                show_rownames    = TRUE,
                show_colnames    = TRUE,
                border_color     = '#F1FAEE',              # 关键：去除格子边框，更清爽！
                cellwidth        = 5,               # 横向时稍微宽一点更好看
                cellheight       = 18,
                fontsize         = 9,
                fontsize_row     = 12,
                fontsize_col     = 7,
                main             = paste0(title_prefix, "\n(mean log10 abundance across stages)"),
                display_numbers  = FALSE,           # 不显示数字
                na_col           = "white",
                legend           = TRUE,
                angle_col        = 90,
                gaps_col         = NULL)
  
  # 10. 保存
  pdf(output_file, width = pdf_width, height = 5.5)
  print(p)
  dev.off()
  
  message("√ 已保存: ", output_file)
  message("   代谢物数量: ", n_metab, "   PDF宽度: ", round(pdf_width, 2), " inches")
}

# ——————————————————————————————————
# 03. 执行绘图
# ——————————————————————————————————
output_dir <- "01.stage_heatmap"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

process_metabolite_heatmap(
  data        = Lipid_abs,
  title_prefix= "Lipid Metabolites",
  output_file = file.path(output_dir, "Lipid_stage_heatmap_horizontal.pdf"),
  de_data     = Lipid_DE
)

process_metabolite_heatmap(
  data        = Metabolin_abs,
  title_prefix= "Polar Metabolites",   # 根据你的实际命名调整
  output_file = file.path(output_dir, "Metabolin_stage_heatmap_horizontal.pdf"),
  de_data     = Metabolin_DE
)

message("所有代谢物阶段热图绘制完成！")