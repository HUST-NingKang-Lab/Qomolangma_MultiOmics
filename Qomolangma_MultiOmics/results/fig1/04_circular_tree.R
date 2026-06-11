

library(ape)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(RColorBrewer)

# =========================
# 1. 读取数据
# =========================
tree <- read.tree("data/processed/user_bac120.unrooted.tree")
phylum_data <- read.csv("data/processed/01_01.Phylum.csv") %>% select(refMAG, Color) %>% mutate(refMAG = str_remove(refMAG, "\\.fa$"))
coverage_data <- read.csv("data/processed/01_02.coverage.csv", check.names = FALSE) %>% mutate(refMAG = str_remove(refMAG, "\\.fa$"))
MinT_MaxT_data <- read.csv("data/processed/01_04.coverage_MinT_MaxT.csv") %>% mutate(Genome = str_remove(Genome, "\\.fa$"))

# Rename columns
colnames(phylum_data) <- c("label", "phylum_color")
colnames(coverage_data) <- c("label", "gut_coverage", "skin_coverage", "oral_coverage")
colnames(MinT_MaxT_data) <- c("label", "oral_MinT", "skin_MinT", "gut_MinT", "oral_MaxT", "skin_MaxT", "gut_MaxT")

# 保留树上叶子
tree_tips <- tree$tip.label
all_data <- reduce(list(phylum_data, coverage_data, MinT_MaxT_data), full_join, by = "label") %>%
  filter(label %in% tree_tips)

meta_df <- all_data %>% column_to_rownames("label") %>% as.data.frame()
meta_df$tip <- rownames(meta_df)

# =========================
# 2. 给树的clade上色（Phylum）
# =========================
phylum_color_df <- phylum_data %>% filter(label %in% tree_tips)
tree <- groupOTU(tree, split(phylum_color_df$label, phylum_color_df$phylum_color))

# 修改角度：350°
p <- ggtree(tree, layout = "fan", open.angle = 10, branch.length = "none", aes(color = group)) +
  scale_color_identity() +
  theme_tree2() +
  theme(legend.position = "none") +
  geom_tiplab(size = 0)

# =========================
# 3. 添加annotation圈
# =========================

## 第1圈：Phylum颜色块
p <- p + geom_fruit(
  data = meta_df,
  geom = geom_tile,
  mapping = aes(y = tip, fill = phylum_color),
  width = 1,
  offset = 0.08,
  color = NA,
  alpha = 0.9
)

## 第2~4圈：覆盖度分层柱状图
coverage_long <- meta_df %>%
  select(tip, gut_coverage, skin_coverage, oral_coverage) %>%
  pivot_longer(cols = -tip, names_to = "Source", values_to = "Coverage")

coverage_colors <- c("oral_coverage" = "#E41A1C",   # 红
                     "skin_coverage" = "#377EB8",  # 蓝
                     "gut_coverage" = "#4DAF4A")  # 绿

for (i in seq_along(coverage_colors)) {
  src <- names(coverage_colors)[i]
  sub_df <- coverage_long %>% filter(Source == src)
  p <- p + geom_fruit(
    data = sub_df,
    geom = geom_bar,
    mapping = aes(y = tip, x = Coverage, fill = Source),
    orientation = "y",
    stat = "identity",
    width = 1.5,
    offset = 0.1,
    color = NA
  )
}

## 第5~10圈：MinT / MaxT 热图
time_levels <- c("T1","T2","Ta","T3","T4","Tb","T5","T6","T7","T8","T9","T10")

time_colors <- c("#EF5350", "#D2E3F1", "#B8D0E7", "#A0BDDB", "#84ABD1", "#6999C7", "#4F86BC", "#3574B1", "#1B61A7", "#014F9C" , "#E53935", "#B71C1C")
names(time_colors) <- time_levels

MinT_MaxT_long <- meta_df %>%
  select(tip, oral_MinT, skin_MinT, gut_MinT, oral_MaxT, skin_MaxT, gut_MaxT) %>%
  pivot_longer(cols = -tip, names_to = "Source", values_to = "TimePoint") %>%
  filter(TimePoint %in% time_levels)

for (i in 1:6) {
  src <- unique(MinT_MaxT_long$Source)[i]
  sub_df <- MinT_MaxT_long %>% filter(Source == src)
  p <- p + geom_fruit(
    data = sub_df,
    geom = geom_tile,
    mapping = aes(y = tip, fill = TimePoint),
    width = 1.6,       # 圈加宽
    offset = 0.12,
    color = NA
  )
}

# =========================
# 4. 统一颜色
# =========================
p <- p + scale_fill_manual(
  values = c(
    setNames(phylum_data$phylum_color, phylum_data$phylum_color),
    coverage_colors,
    time_colors
  ),
  na.value = "white"
)

# =========================
# 5. 美化主题并保存
# =========================
p <- p + theme(
  plot.margin = margin(15, 15, 15, 15),
  panel.background = element_blank(),   # 去掉背景色
  plot.background = element_blank(),    # 去掉背景色
  legend.position = "right"             # 添加图例
)

ggsave("annotated_tree.pdf", p, width = 20, height = 20, dpi = 300)



######### 打印从里到外的变量
for (i in 1:6) {
  src <- unique(MinT_MaxT_long$Source)[i]
  print(src)   # 打印当前层对应的变量
  sub_df <- MinT_MaxT_long %>% filter(Source == src)
  p <- p + geom_fruit(
    data = sub_df,
    geom = geom_tile,
    mapping = aes(y = tip, fill = TimePoint),
    width = 1.6,
    offset = 0.12,
    color = NA
  )
}

# [1] "oral_MinT"
# [1] "oral_MaxT"
# [1] "skin_MinT"
# [1] "skin_MaxT"
# [1] "gut_MinT"
# [1] "gut_MaxT"

for (i in seq_along(coverage_colors)) {
  src <- names(coverage_colors)[i]
  print(src)   # 打印当前圈对应的变量
  sub_df <- coverage_long %>% filter(Source == src)
  p <- p + geom_fruit(
    data = sub_df,
    geom = geom_bar,
    mapping = aes(y = tip, x = Coverage, fill = Source),
    orientation = "y",
    stat = "identity",
    width = 1.5,
    offset = 0.1,
    color = NA
  )
}

# [1] "oral_coverage"
# [1] "skin_coverage"
# [1] "gut_coverage"
