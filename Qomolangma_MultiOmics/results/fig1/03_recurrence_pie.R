library(vegan)
library(ggplot2)
library(patchwork)

# 设置工作目录


# 加载数据
meta <- read.csv("data/metadata/Matedata_Information.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time12", "people")]
colnames(meta)[1] <- "SampleID"
meta$time12 <- factor(meta$time12, levels = c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6", "T7", "T8", "T9", "T10"))
meta$people <- factor(meta$people)

# 读取微生物群落丰度数据（标准化后的MAG绝对丰度）
abs <- read.csv("data/processed/02.MAG_abs_composition_normalized.csv", check.names = FALSE, row.names = 1)

# 提取皮肤样本的丰度数据
abs_skin <- abs[, colnames(abs) %in% meta$SampleID[meta$Source == "skin"]]
# 过滤低丰度和低发生率的MAG：总丰度 > 0.001 且在至少20%的样本中出现
abs_skin <- abs_skin[rowSums(abs_skin) > 0.001 & rowSums(abs_skin > 0) >= 0.2 * ncol(abs_skin), ]

# 提取肠道样本的丰度数据
abs_gut <- abs[, colnames(abs) %in% meta$SampleID[meta$Source == "gut"]]
# 过滤低丰度和低发生率的MAG
abs_gut <- abs_gut[rowSums(abs_gut) > 0.001 & rowSums(abs_gut > 0) >= 0.2 * ncol(abs_gut), ]

# 提取口腔样本的丰度数据
abs_oral <- abs[, colnames(abs) %in% meta$SampleID[meta$Source == "oral"]]
# 过滤低丰度和低发生率的MAG
abs_oral <- abs_oral[rowSums(abs_oral) > 0.001 & rowSums(abs_oral > 0) >= 0.2 * ncol(abs_oral), ]

# 定义CAP分析及绘图函数
cap_plot <- function(abundance, meta_subset, source_name) {
  # 确保样本顺序与丰度矩阵列名一致
  meta_subset <- meta_subset[match(colnames(abundance), meta_subset$SampleID), ]
  
  # 计算Bray-Curtis距离矩阵
  dist_matrix <- vegdist(t(abundance), method = "bray")
  
  # 执行CAP分析：控制个体差异（people），时间点（time12）为解释变量
  cap_model <- capscale(dist_matrix ~ time12 + Condition(people), data = meta_subset)
  
  # 显著性检验（置换检验）
  perm_test <- anova(cap_model, permutations = 999)
  p_value <- perm_test[["Pr(>F)"]][1]
  
  # 提取CAP轴解释的方差比例
  cap_summary <- summary(cap_model)
  cap1_var <- round(100 * cap_summary$cont$importance[2, 1], 1)  # CAP1解释比例
  cap2_var <- round(100 * cap_summary$cont$importance[2, 2], 1)  # CAP2解释比例
  
  # 提取样本在CAP空间中的坐标
  scores <- as.data.frame(scores(cap_model)$sites)
  scores$time12 <- meta_subset$time12
  scores$people <- meta_subset$people
  
  # 定义显著性标记
  sig_label <- ifelse(p_value < 0.001, "***",
                      ifelse(p_value < 0.01, "**",
                             ifelse(p_value < 0.05, "*", "ns")))
  # 组合P值和显著性符号
  sig_text <- paste("P =", sprintf("%.3f", p_value), sig_label)
  
  # 使用ggplot2绘制CAP图
  p <- ggplot(scores, aes(x = CAP1, y = CAP2, color = time12, shape = people)) +
    geom_point(size = 3, alpha = 1) +
    labs(title = source_name,
         x = paste("CAP1 (", cap1_var, "%)"),
         y = paste("CAP2 (", cap2_var, "%)")) +
    scale_color_viridis_d(name = "Time12") +  # 使用viridis调色板，图例标题为“时间点”
    scale_shape_manual(name = "", values = 1:length(unique(meta_subset$people))) +
    # 在右上角添加显著性文本
    annotate("text", x = Inf, y = Inf, label = sig_text, hjust = 1.1, vjust = 1.1, size = 8) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 15, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 14),
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 14),
      legend.box = "horizontal",
      legend.margin = margin(t = 0, b = 10)
    )
  
  # 返回结果：绘图对象和P值
  return(list(plot = p, p_value = p_value, source = source_name))
}

# 为每种样本类型生成CAP分析图
meta_skin <- meta[meta$Source == "skin", ]
meta_gut <- meta[meta$Source == "gut", ]
meta_oral <- meta[meta$Source == "oral", ]

# 分别对皮肤、肠道、口腔数据进行分析
result_skin <- cap_plot(abs_skin, meta_skin, "skin")
result_gut <- cap_plot(abs_gut, meta_gut, "gut")
result_oral <- cap_plot(abs_oral, meta_oral, "oral")

# 创建显著性结果的数据框
significance_data <- data.frame(
  样本类型 = c(result_skin$source, result_gut$source, result_oral$source),
  P值 = c(result_skin$p_value, result_gut$p_value, result_oral$p_value),
  显著性 = c(
    ifelse(result_skin$p_value < 0.001, "***",
           ifelse(result_skin$p_value < 0.01, "**",
                  ifelse(result_skin$p_value < 0.05, "*", "ns"))),
    ifelse(result_gut$p_value < 0.001, "***",
           ifelse(result_gut$p_value < 0.01, "**",
                  ifelse(result_gut$p_value < 0.05, "*", "ns"))),
    ifelse(result_oral$p_value < 0.001, "***",
           ifelse(result_oral$p_value < 0.01, "**",
                  ifelse(result_oral$p_value < 0.05, "*", "ns")))
  )
)

# 提取各个子图
p_skin <- result_skin$plot
p_gut <- result_gut$plot
p_oral <- result_oral$plot

# 使用patchwork将三个图组合在一起
layout_design <- "
ABC
"

# 组合图形，统一图例并设置布局
combined_plot <- p_oral + p_gut + p_skin +
  plot_layout(design = layout_design, heights = c(4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

# 保存组合图（PDF格式）
ggsave("12个时间点的CAP.pdf", combined_plot, width = 18, height = 8, units = "in")

# 输出显著性检验结果
print("显著性检验结果:")
print(significance_data)