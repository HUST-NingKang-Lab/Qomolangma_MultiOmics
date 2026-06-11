# 功能：读取过滤后的 abs → 按 time12 分组 → 筛选物种 → 计算 alpha 多样性 → 绘制箱线图
# 依赖：vegan, ggplot2, ggsignif, gridExtra, ggpubr, ggsci, dplyr, purrr, grid

library(vegan)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(ggpubr)
library(ggsci)
library(dplyr)
library(purrr)
library(grid)
source("code/00.theme.R")

dir_data <- "data/processed"
out_dir <- "data/processed"

# 读取 metadata ---------------------------------------------------------
meta <- read.csv("data/metadata/Matedata_Information.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time5", "time12")]
colnames(meta)[1] <- "sampleid"

# 读取所有过滤后的 TPM 并计算多样性 -------------------------------------
files <- list.files(dir_data, pattern = "^TPM_filtered_.+\\.rds$", full.names = TRUE)
diversity_all <- NULL

for (f in files) {
  src <- sub("TPM_filtered_(.+?)\\.rds$", "\\1", basename(f))
  TPM_fil <- readRDS(f)
  
  # 按 time12 分组
  time12_groups <- unique(meta$time12[meta$Source == src & meta$time12 != "Tc"])
  
  for (t in time12_groups) {
    # 获取当前 time12 组的样本
    samples <- meta$sampleid[meta$time12 == t & meta$Source == src]
    TPM_sub <- TPM_fil[, samples, drop = FALSE]
    
    if (nrow(TPM_sub) > 0) {
      abs_t <- t(TPM_sub)
      shan  <- diversity(abs_t, "shannon")
      simp  <- 1 - diversity(abs_t, "simpson")
      sr    <- specnumber(abs_t)
      piel  <- shan / log(sr)
      
      df <- data.frame(sampleid = rownames(abs_t),
                      Shannon.wiener = shan,
                      simpson = simp,
                      SR = sr,
                      Pielou = piel,
                      Source = src,
                      time12 = t)
      diversity_all <- bind_rows(diversity_all, df)
    }
  }
}

# 保存多样性结果 --------------------------------------------------------
write.csv(diversity_all,
          file.path(out_dir, "01.diversity_time12.csv"),
          row.names = FALSE)

# 合并 metadata ---------------------------------------------------------
plotdata <- inner_join(meta, diversity_all[, c("Shannon.wiener", "sampleid", "Source", "time12")],
                      by = c("sampleid", "Source", "time12")) %>%
  filter(time12 != "Tc") %>%
  mutate(time12 = factor(time12,
                        levels = c("T1", "T2", "Ta", "Tb", "T4", "T4", "T5", "T6", "T7", "T8", "T9", "T10")))

# 定义不同Source的颜色配置 ----------------------------------------------
region_colors <- list(
  "gut" = c("#e10000", "#e31515", "#e62a2a", "#e83e3f", "#eb5354", "#ed6869", 
            "#f29292", "#f29292", "#f5a7a7", "#f7bbbc", "#fad0d1", "#fce5e6"),
  "oral" = c("#0070ca", "#157cce", "#2a87d3", "#3e93d7", "#539fdb", "#68aadf", 
             "#7db6e4", "#92c1e8", "#a7cdec", "#bbd9f0", "#d0e4f5", "#e5f0f9"),
  "skin" = c("#f1b500", "#f2bb15", "#f4c12a", "#f5c73e", "#f6cd53", "#f7d368", 
             "#f9da7d", "#fae092", "#fbe6a7", "#fcecbb", "#fef2d0", "#fff8e5")
)

# 获取颜色函数 ----------------------------------------------------------
get_colors_for_source <- function(src) {
  time12_levels <- c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6", "T7", "T8", "T9", "T10")
  colors <- region_colors[[src]]
  if (is.null(colors)) {
    # 如果没有匹配的颜色，使用默认颜色
    colors <- rainbow(length(time12_levels))
  }
  # 确保颜色数量匹配时间点数量
  colors <- colors[1:length(time12_levels)]
  names(colors) <- time12_levels
  return(colors)
}

# 两两比较组合 ----------------------------------------------------------
time12_levels <- levels(plotdata$time12)
combn_matrix <- combn(time12_levels, 2)
compaired <- split(combn_matrix, col(combn_matrix)) |> lapply(as.character)

# 绘图函数 --------------------------------------------------------------
make_plot <- function(src) {
  sub <- filter(plotdata, Source == src)
  if (nrow(sub) == 0) return(NULL)
  
  # 获取当前source对应的颜色
  current_colors <- get_colors_for_source(src)

  # 计算显著性差异
  significant_comparisons <- lapply(compaired, function(comp) {
    group1 <- sub$Shannon.wiener[sub$time12 == comp[1]]
    group2 <- sub$Shannon.wiener[sub$time12 == comp[2]]
    if (length(group1) > 0 && length(group2) > 0) {
      p_value <- wilcox.test(group1, group2)$p.value
      if (p_value < 0.05) return(comp)
    }
    return(NULL)
  }) %>% compact()  # 移除NULL元素

  p <- ggplot(sub, aes(time12, Shannon.wiener, fill = time12)) +
    geom_boxplot() +
    theme_publish() +
    labs(title = paste("Alpha Diversity -", src),
         x = "Time Point", y = "Shannon Diversity") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = current_colors)

  # 仅当存在显著性差异时添加显著性条
  if (length(significant_comparisons) > 0) {
    p <- p + stat_compare_means(aes(group = time12),
                               method = "wilcox.test",
                               label = "p.signif",
                               comparisons = significant_comparisons)
  }

  return(p)
}

# 生成并保存图 ----------------------------------------------------------
# source_levels <- unique(plotdata$Source)
source_levels <- c("oral", "gut", "skin")
plot_list <- keep(source_levels, ~ !is.null(make_plot(.x))) %>%
  lapply(make_plot)

if (length(plot_list) > 0) {
  p_all <- grid.arrange(grobs = plot_list, ncol = length(plot_list))

  pdf("Figures/02.MAG_analysis/02.群落结构分析/02.计算alpha多样性/01.time12_shannon_only_sign.pdf",
      width = 10 * length(plot_list), height = 10)
  grid.draw(p_all)
  dev.off()
  cat("箱线图已保存！\n")
} else {
  cat("警告：没有足够数据绘图\n")
}

cat("多样性数据已保存至:", file.path(out_dir, "01.diversity_time12.csv"), "\n")