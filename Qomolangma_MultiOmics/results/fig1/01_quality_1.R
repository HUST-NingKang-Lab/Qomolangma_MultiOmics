# 加载必要的R包
library(dplyr)
library(plotrix)
library(RColorBrewer)

# 设置工作目录

# 读取数据
plotdata <- read.csv("data/processed/09.MAG_stat_all_20250821.csv", check.names=FALSE) %>% 
  select(MAG, Source, "Quality(MIMAG)")

# 定义颜色映射
color_map <- list(
  oral = "#46B8DAFF",
  gut = "#D43F3AFF",
  skin = "#EEA236FF"
)

# 创建浅一度和浅两度的颜色函数
lighten_color <- function(color, factor) {
  col <- col2rgb(color)
  col <- rgb(t(pmin(255, col + (255 - col) * factor)), maxColorValue = 255)
  return(col)
}

# 为禁止在 xAI 提示中定义颜色
colors <- list(
  oral = list(
    High = lighten_color(color_map$oral, 0.6),   # 浅一度
    Medium = lighten_color(color_map$oral, 0.3)  # 浅两度
  ),
  gut = list(
    High = lighten_color(color_map$gut, 0.6),
    Medium = lighten_color(color_map$gut, 0.3)
  ),
  skin = list(
    High = lighten_color(color_map$skin, 0.6),
    Medium = lighten_color(color_map$skin, 0.3)
  )
)

# 创建输出目录（如果不存在）
output_dir <- "04.MAG目录概述"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 设置PDF输出
pdf(file.path(output_dir, "01.饼状图_绘制三个部位高质量和中等质量MAG数量.pdf"), width = 15, height = 5)

# 设置绘图布局为1行3列
par(mfrow = c(1, 3), mar = c(2, 2, 4, 2))

# 为每个Source绘制饼状图
for (source in c("oral", "gut", "skin")) {
  # 筛选当前Source的数据
  data_subset <- plotdata %>% filter(Source == source)
  
  # 计算Quality(MIMAG)的分布
  quality_counts <- table(data_subset$`Quality(MIMAG)`)
  
  # 准备饼状图数据
  labels <- names(quality_counts)
  values <- as.numeric(quality_counts)
  
  # 确保只有High和Medium类别
  valid_labels <- intersect(labels, c("High", "Medium"))
  valid_values <- values[labels %in% c("High", "Medium")]
  valid_colors <- sapply(valid_labels, function(x) colors[[source]][[x]])
  
  # 绘制3D饼状图
  pie3D(
    valid_values,
    labels = paste0(valid_labels, "\n(", round(valid_values/sum(valid_values)*100, 1), "%)"),
    col = valid_colors,
    theta = pi/4,  # 3D角度
    radius = 0.8,  # 饼图大小
    height = 0.06,  # 3D厚度
    main = paste("Quality Distribution for", source),
    labelcex = 1.2,  # 标签字体大小
    shade = 0.7,     # 阴影效果
    explode = 0.1,   # 扇形间隙
    start = pi/5     # 左右旋转饼图
  )
}

# 关闭PDF设备
dev.off()