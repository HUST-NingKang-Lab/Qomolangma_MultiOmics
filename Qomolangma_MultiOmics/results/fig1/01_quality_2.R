# 设置工作目录


# 加载必要的R包
library(dplyr)
library(ggplot2)
library(cowplot)
library(gghalves)
library(ggdist)
library(patchwork)

# 读取数据
plotdata <- read.csv("data/processed/09.MAG_stat_all_20250821.csv") %>% 
    select(MAG, Source, Completeness, Contamination, GC, avg_contig_len, N50) %>%
    mutate(
        avg_contig_len = as.numeric(gsub(",", "", avg_contig_len)),
        N50 = as.numeric(gsub(",", "", N50))
    ) %>% 
    mutate(
        avg_contig_len = log10(avg_contig_len),
        N50 = log10(N50)
    )

# 定义颜色方案
colors <- c("oral" = "#46B8DAFF", "gut" = "#D43F3AFF", "skin" = "#EEA236FF")

# 创建云雨图绘制函数
create_raincloud <- function(data, y_var, title) {
  ggplot(data, aes_string("Source", y_var, fill = "Source", color = "Source")) +
    # 绘制半边小提琴图（密度分布）
    geom_half_violin(color = "white", side = "right", alpha = 0.7, position = position_nudge(x = -2.5)) +
    # 在小提琴图下边缘添加很窄的箱线图，属于小提琴图的一部分
    geom_boxplot(color = "black", fill = "white", width = 0.07, staplewidth = 0.2, 
                 outlier.shape = NA, position = position_nudge(x = -2.5, y = 0),
                 show.legend = FALSE) +
    # 绘制原有的箱线图（保持原位置）
    geom_boxplot(color = "black", width = 0.1, staplewidth = 0.4, 
                 outlier.shape = NA, position = position_nudge(x = -2.8)) +
    # 绘制散点图
    geom_half_point(range_scale = 0.3, alpha = 0.6,
                    position = position_nudge(x = -3)) +
    # 应用颜色方案
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    # 翻转坐标系
    coord_flip(clip = "off") + 
    # 设置标签和主题
    labs(x = NULL, y = y_var, title = title) + 
    theme_minimal() +
    theme(
      # 调整plot.margin，减少上下边距，让图形更靠近x轴
      plot.margin = unit(c(0.2, 0.5, 0.1, 0.5), units = "cm"),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = unname(colors), size = 11), 
      axis.line.x = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      # 减少轴标题与轴的距离
      axis.title.x = element_text(margin = margin(t = 5)),
      # 减少面板与轴的间距
      panel.spacing = unit(0.1, "cm")
    )
}

# 为每个Source和每个指标绘制云雨图
sources <- c("oral", "gut", "skin")
metrics <- c("Completeness", "Contamination", "GC", "avg_contig_len", "N50")
plots_list <- list()

# 生成所有子图 (5行3列：5个指标 x 3个Source)
for (i in 1:length(metrics)) {
  for (j in 1:length(sources)) {
    metric <- metrics[i]
    source <- sources[j]
    
    # 过滤当前Source的数据
    source_data <- plotdata %>% filter(Source == source)
    
    # 创建云雨图
    p <- create_raincloud(source_data, metric, paste(source, "-", metric))
    
    # 将图添加到列表中，使用行列索引命名
    plots_list[[paste0("r", i, "c", j)]] <- p
  }
}

# 使用patchwork组合所有图形（5行3列布局），减少子图间距
combined_plot <- (plots_list[["r1c1"]] | plots_list[["r1c2"]] | plots_list[["r1c3"]]) /
                 (plots_list[["r2c1"]] | plots_list[["r2c2"]] | plots_list[["r2c3"]]) /
                 (plots_list[["r3c1"]] | plots_list[["r3c2"]] | plots_list[["r3c3"]]) /
                 (plots_list[["r4c1"]] | plots_list[["r4c2"]] | plots_list[["r4c3"]]) /
                 (plots_list[["r5c1"]] | plots_list[["r5c2"]] | plots_list[["r5c3"]]) &
                 theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# 保存为PDF文件
ggsave("云雨图_MAG质量.pdf",
       plot = combined_plot,
       width = 13, height = 15, 
       units = "in", dpi = 300)

# 显示最终图形
print(combined_plot)