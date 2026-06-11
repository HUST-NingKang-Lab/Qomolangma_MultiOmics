
library(dplyr)
library(ggplot2)
library(ComplexUpset)

MAG_stat <- read.csv("data/processed/09.MAG_stat_all_20250821.csv", check.names = FALSE) %>% 
    select(refMAG, Source) %>% 
    unique()
    
# 加载自定义主题
source("code/00.theme.R")

# 转宽表
wide_data <- MAG_stat %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = Source, values_from = present, values_fill = FALSE) %>%
  select(refMAG, skin, gut, oral)

# 绘图
p <- upset(
  wide_data,
  intersect = c("skin", "gut", "oral"),
  name = "refMAG Overlap",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 4, family = "sans", color = "black")
    ) +
    labs(y = "Intersection Size") +
    theme_publish() +
    theme(
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
  ),
  set_sizes = upset_set_size() +
    labs(y = "Set Size") +
    theme_publish() +
    theme(
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    ),
  width_ratio = 0.25,
  sort_sets = FALSE,
  themes = theme_publish()
)


# 展示
print(p)

# 保存图形
ggsave("Upset图_不同部位的genome_species数量.pdf",
       plot = p, width = 9, height = 5.5, dpi = 300)


