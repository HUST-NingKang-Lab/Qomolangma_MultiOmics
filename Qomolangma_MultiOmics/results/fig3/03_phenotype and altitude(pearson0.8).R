
# 加载必要的库
library(dplyr)
library(reshape2)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)

# —————————————————————————————————————————————
# 01. 读取并处理数据
# —————————————————————————————————————————————
people <- paste0("S", 1:13)
time12 <- c("Ta", "Tb", paste0("T", 1:10))

# 读取和海拔高度相关的表型
sig_phenotype_altitude <- read.csv("data/processed/significant_phenotype_mean_correlations.csv", check.names = FALSE) %>% 
    filter(abs(mean_pearson_correlation) > 0.8 & pearson_ttest_fdr < 0.05) %>% 
    select(data_type, phenotype, n_subjects, mean_pearson_correlation)

# 读取表型数据并过滤
physiology_abs <- read.csv("data/processed/physiology.csv", check.names = FALSE) %>%
  select(subject_id, time, all_of(unique(sig_phenotype_altitude[sig_phenotype_altitude$data_type == "physiology",]$phenotype))) %>% 
  filter(subject_id %in% people & time %in% time12)
eye_abs <- read.csv("data/processed/eye.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(sig_phenotype_altitude[sig_phenotype_altitude$data_type == "eye",]$phenotype))) %>% 
  filter(subject_id %in% people & time %in% time12) 
cognitive_abs <- read.csv("data/processed/cognitive.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(sig_phenotype_altitude[sig_phenotype_altitude$data_type == "cognitive",]$phenotype))) %>% 
  filter(subject_id %in% people & time %in% time12)
ultrasound_abs <- read.csv("data/processed/ultrasound.csv", check.names = FALSE) %>% 
  select(subject_id, time, all_of(unique(sig_phenotype_altitude[sig_phenotype_altitude$data_type == "ultrasound",]$phenotype))) %>% 
  filter(subject_id %in% people & time %in% time12)

# 读取时间和海拔对应关系
pair_time_altitude <- read.csv("data/metadata/Matedata_Information.csv", check.names = FALSE) %>% 
    select(time12, altitude) %>% 
    unique()
colnames(pair_time_altitude)[1] <- "time"

# —————————————————————————————————————————————
# 02. 合并所有表型数据
# —————————————————————————————————————————————
# 合并所有数据集
all_data <- bind_rows(
  physiology_abs %>% mutate(data_type = "physiology"),
  eye_abs %>% mutate(data_type = "eye"),
  cognitive_abs %>% mutate(data_type = "cognitive"),
  ultrasound_abs %>% mutate(data_type = "ultrasound")
) %>%
  distinct()

# 检查数据
cat("合并后的数据维度:", dim(all_data), "\n")
cat("可用的表型:", setdiff(names(all_data), c("subject_id", "time", "data_type")), "\n")

# —————————————————————————————————————————————
# 03. 计算每个时间点的平均表型值
# —————————————————————————————————————————————
# 定义时间顺序
time_order <- c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6", "T7", "T8", "T9", "T10")

# 计算每个时间点的平均表型值
calculate_mean_values <- function(data) {
  # 确保时间按正确顺序排列
  data$time <- factor(data$time, levels = time_order)
  
  # 获取表型列名（除了subject_id、time和data_type）
  phenotype_cols <- names(data)[!names(data) %in% c("subject_id", "time", "data_type")]
  
  # 计算每个表型在每个时间点的均值
  mean_values <- data.frame()
  
  for (time_point in time_order) {
    time_data <- data %>% filter(time == time_point)
    
    if (nrow(time_data) > 0) {
      for (phenotype in phenotype_cols) {
        # 跳过没有数据的表型
        if (!phenotype %in% names(time_data)) {
          next
        }
        
        phenotype_values <- time_data[[phenotype]]
        
        # 确保数据是数值型，过滤掉NA值和非数值
        phenotype_values <- as.numeric(phenotype_values)
        valid_values <- phenotype_values[!is.na(phenotype_values)]
        
        if (length(valid_values) > 0) {
          mean_value <- mean(valid_values, na.rm = TRUE)
          
          mean_values <- rbind(mean_values, data.frame(
            phenotype = phenotype,
            time = time_point,
            mean_value = mean_value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  return(mean_values)
}

# 计算平均值
mean_values <- calculate_mean_values(all_data)

# 检查结果
cat("计算出的平均值数量:", nrow(mean_values), "\n")
cat("涉及的表型:", unique(mean_values$phenotype), "\n")
cat("涉及的时间点:", unique(mean_values$time), "\n")

# 如果没有数据，停止执行
if (nrow(mean_values) == 0) {
  stop("没有计算出任何平均值，请检查数据筛选条件")
}

# —————————————————————————————————————————————
# 04. 数据标准化和热图矩阵准备
# —————————————————————————————————————————————
# 创建热图矩阵
heatmap_matrix <- mean_values %>%
  dcast(phenotype ~ time, value.var = "mean_value")

# 检查矩阵
cat("热图矩阵维度:", dim(heatmap_matrix), "\n")
cat("矩阵中的表型:", heatmap_matrix$phenotype, "\n")

# 将表型名设为行名
rownames(heatmap_matrix) <- heatmap_matrix$phenotype
heatmap_matrix$phenotype <- NULL

# 确保列按时间顺序排列
available_times <- intersect(time_order, colnames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[, available_times, drop = FALSE]

cat("最终矩阵维度:", dim(heatmap_matrix), "\n")
cat("可用时间点:", colnames(heatmap_matrix), "\n")

# 检查是否有空矩阵
if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
  stop("热图矩阵为空，请检查数据")
}

# 标准化到[-1, 1]
normalize_to_range <- function(x) {
  if (all(is.na(x))) return(x)
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (min_val == max_val) return(rep(0, length(x)))
  2 * (x - min_val) / (max_val - min_val) - 1
}

# 对每行（表型）进行标准化
heatmap_matrix_normalized <- t(apply(heatmap_matrix, 1, normalize_to_range))

# 检查标准化后的矩阵
cat("标准化后矩阵维度:", dim(heatmap_matrix_normalized), "\n")
cat("是否有无穷值:", any(is.infinite(heatmap_matrix_normalized)), "\n")
cat("是否有NA值:", any(is.na(heatmap_matrix_normalized)), "\n")

# 处理可能的NA或无穷值
heatmap_matrix_normalized[is.na(heatmap_matrix_normalized)] <- 0
heatmap_matrix_normalized[is.infinite(heatmap_matrix_normalized)] <- 0

# 移除所有值都相同的行（这些行无法进行聚类）
row_vars <- apply(heatmap_matrix_normalized, 1, var, na.rm = TRUE)
valid_rows <- !is.na(row_vars) & row_vars > 0
heatmap_matrix_normalized <- heatmap_matrix_normalized[valid_rows, , drop = FALSE]

cat("最终用于绘图的矩阵维度:", dim(heatmap_matrix_normalized), "\n")

# —————————————————————————————————————————————
# 05. 准备注释数据
# —————————————————————————————————————————————
# 准备海拔注释数据 - 直接使用时间点对应的海拔
available_times <- colnames(heatmap_matrix_normalized)
altitude_annotation <- pair_time_altitude %>%
  filter(time %in% available_times) %>%
  arrange(match(time, available_times))

# 海拔数据归一化
altitude_annotation$altitude_norm <- (altitude_annotation$altitude - min(altitude_annotation$altitude, na.rm = TRUE)) / 
  (max(altitude_annotation$altitude, na.rm = TRUE) - min(altitude_annotation$altitude, na.rm = TRUE))

# 准备表型数量和相关性注释
phenotype_count <- sig_phenotype_altitude %>%
  select(phenotype, n_subjects, mean_pearson_correlation, data_type) %>%
  filter(phenotype %in% rownames(heatmap_matrix_normalized))

# 检查注释数据
cat("sig_phenotype_altitude中的表型:", unique(sig_phenotype_altitude$phenotype), "\n")
cat("热图矩阵中的表型:", rownames(heatmap_matrix_normalized), "\n")
cat("匹配的表型数量:", nrow(phenotype_count), "\n")

# 确保顺序一致
phenotype_count <- phenotype_count[match(rownames(heatmap_matrix_normalized), phenotype_count$phenotype), ]

# 检查是否有NA值
cat("phenotype_count中是否有NA:", any(is.na(phenotype_count)), "\n")
cat("n_subjects:", phenotype_count$n_subjects, "\n")
cat("mean_pearson_correlation:", phenotype_count$mean_pearson_correlation, "\n")
cat("data_type:", phenotype_count$data_type, "\n")

# 处理可能的NA值并确保数据完整性
phenotype_count$n_subjects[is.na(phenotype_count$n_subjects)] <- 0
phenotype_count$mean_pearson_correlation[is.na(phenotype_count$mean_pearson_correlation)] <- 0
phenotype_count$data_type[is.na(phenotype_count$data_type)] <- "unknown"

# 移除第一行如果它对应的是无效数据
if (nrow(phenotype_count) > nrow(heatmap_matrix_normalized)) {
  phenotype_count <- phenotype_count[1:nrow(heatmap_matrix_normalized), ]
}

# —————————————————————————————————————————————
# 06. 绘制热图
# —————————————————————————————————————————————
# 定义配色
color_palette <- colorRampPalette(c("#4281A4", "#E4DFDA", "#C1666B"))(1000)
altitude_palette <- colorRampPalette(c("#2E8B57", "#FFD700", "#CD853F"))(100)

# 创建输出目录
output_dir <- getwd()
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 使用ComplexHeatmap绘制
library(ComplexHeatmap)

# 创建海拔顶部注释（作为列注释）
ha_top <- HeatmapAnnotation(
  "Altitude" = altitude_annotation$altitude_norm,
  col = list("Altitude" = colorRamp2(c(0, 0.5, 1), c("#2E8B57", "#FFD700", "#CD853F"))),
  annotation_name_gp = gpar(fontsize = 14),
  annotation_legend_param = list("Altitude" = list(title = "Normalized\nAltitude", title_gp = gpar(fontsize = 12), labels_gp = gpar(fontsize = 11))),
  height = unit(0.8, "cm")
)

# 为data_type创建颜色映射
data_type_colors <- c(
  "physiology" = "#E74C3C",    # 红色
  "eye" = "#3498DB",           # 蓝色
  "cognitive" = "#2ECC71",     # 绿色
  "ultrasound" = "#9B59B6"     # 紫色
)

# 创建表型数量和相关性柱状图注释
ha_left <- rowAnnotation(
  "Correlation" = create_barplot_with_labels(phenotype_count$mean_pearson_correlation, "#FF7F50"),
  "Cor Values" = anno_text(
    sprintf("%.3f", phenotype_count$mean_pearson_correlation),
    gp = gpar(fontsize = 9, col = "black"),
    width = unit(1.5, "cm")
  ),
  "N Subjects" = create_barplot_with_labels(phenotype_count$n_subjects, "#69B3A2"),
  "N Values" = anno_text(
    as.character(phenotype_count$n_subjects),
    gp = gpar(fontsize = 9, col = "black"),
    width = unit(1, "cm")
  ),
  "Data Type" = phenotype_count$data_type,
  col = list("Data Type" = data_type_colors),
  annotation_name_gp = gpar(fontsize = 12),
  annotation_legend_param = list(
    "Data Type" = list(
      title = "Data Type",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  ),
  gap = unit(2, "mm")
)

# 检查注释数据长度
cat("注释数据长度 - n_subjects:", length(phenotype_count$n_subjects), "\n")
cat("注释数据长度 - correlation:", length(phenotype_count$mean_pearson_correlation), "\n")
cat("注释数据长度 - data_type:", length(phenotype_count$data_type), "\n")
cat("热图行数:", nrow(heatmap_matrix_normalized), "\n")

# 创建主热图
ht_main <- Heatmap(
  heatmap_matrix_normalized,
  name = "Mean Value",
  col = colorRamp2(c(-1, 0, 1), c("#4281A4", "#E4DFDA", "#C1666B")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,  # 不显示行聚类树
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 14),
  row_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Normalized\nMean Value",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 11)
  ),
  top_annotation = ha_top,
  left_annotation = ha_left,
  rect_gp = gpar(col = "grey85", lwd = 0.5)  # 添加灰色边框
)

# 组合热图
ht_combined <- ht_main

# 保存热图
png(file.path(output_dir, "phenotype_mean_values_heatmap.png"), 
    width = 16, height = 12, units = "in", res = 300)
draw(ht_combined, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(file.path(output_dir, "phenotype_mean_values_heatmap.pdf"), 
    width = 16, height = 12)
draw(ht_combined, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

# —————————————————————————————————————————————
# 07. 保存相关数据
# —————————————————————————————————————————————
# 保存平均值矩阵
write.csv(heatmap_matrix, file.path(output_dir, "phenotype_mean_values.csv"))
write.csv(heatmap_matrix_normalized, file.path(output_dir, "phenotype_mean_values_normalized.csv"))

# 保存注释数据
write.csv(altitude_annotation, file.path(output_dir, "altitude_annotation.csv"), row.names = FALSE)
write.csv(phenotype_count, file.path(output_dir, "phenotype_subject_counts.csv"), row.names = FALSE)

# 输出摘要信息
cat("热图绘制完成！\n")
cat("输出文件保存在:", output_dir, "\n")
cat("包含表型数量:", nrow(heatmap_matrix_normalized), "\n")
cat("时间点数量:", ncol(heatmap_matrix_normalized), "\n")
cat("热图文件: phenotype_mean_values_heatmap.png 和 phenotype_mean_values_heatmap.pdf\n")