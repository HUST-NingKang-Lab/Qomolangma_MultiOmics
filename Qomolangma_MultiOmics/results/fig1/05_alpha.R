
library(vegan)
library(ggplot2)
library(dplyr)
library(purrr)
library(grid)
library(gridExtra)
library(ggpubr)
library(ggsci)
source("code/00.theme.R")

# 输入/输出目录 ---------------------------------------------------------
dir_data <- "data/processed"
out_dir <- "alpha多样性"

# 读取 metadata ---------------------------------------------------------
meta <- read.csv("data/metadata/Matedata_Information.csv", check.names = FALSE)
meta <- meta[, c("SampleID", "Source", "time12", "people", "altitude")]
colnames(meta)[1] <- "sampleid"
meta$time12 <- factor(meta$time12,
                      levels = c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5", "T6",
                                 "T7", "T8", "T9", "T10"))
meta$people <- factor(meta$people)

# 读取标准化后的 MAG 绝对丰度数据 ----------------------------------------
abs <- read.csv(file.path(dir_data, "02.MAG_abs_composition_normalized.csv"),
                check.names = FALSE, row.names = 1)

# 根据 Source 提取数据集（皮肤/肠道/口腔）
abs_list <- list(
  skin = abs[, colnames(abs) %in% meta$sampleid[meta$Source == "skin"]],
  gut  = abs[, colnames(abs) %in% meta$sampleid[meta$Source == "gut"]],
  oral = abs[, colnames(abs) %in% meta$sampleid[meta$Source == "oral"]]
)

# 过滤低丰度和低发生率的 MAG
filter_MAG <- function(df) {
  df[rowSums(df) > 0.001 & rowSums(df > 0) >= 0.2 * ncol(df), , drop = FALSE]
}
abs_list <- lapply(abs_list, filter_MAG)

# 计算 alpha 多样性 -----------------------------------------------------
diversity_all <- NULL

for (src in names(abs_list)) {
  abs_fil <- abs_list[[src]]
  if (ncol(abs_fil) == 0) next
  
  # 遍历 time12
  time12_groups <- unique(meta$time12[meta$Source == src & meta$time12 != "Tc"])
  
  for (t in time12_groups) {
    samples <- meta$sampleid[meta$time12 == t & meta$Source == src]
    abs_sub <- abs_fil[, colnames(abs_fil) %in% samples, drop = FALSE]
    
    # 再次筛选物种（避免 time12 内出现空表）
    abs_sub <- abs_sub[rowSums(abs_sub) > 0.001 & rowSums(abs_sub > 0) >= 3, , drop = FALSE]
    
    if (nrow(abs_sub) > 0) {
      abs_t <- t(abs_sub)
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
          file.path("data/processed", "02.绘制alpha多样性的个体变化_alpha_table.csv"),
          row.names = FALSE)

# 合并 metadata ---------------------------------------------------------
plotdata <- inner_join(meta, 
                       diversity_all[, c("Shannon.wiener", "sampleid", "Source", "time12")],
                       by = c("sampleid", "Source", "time12")) %>%
  filter(time12 != "Tc") %>%
  mutate(time12 = factor(time12,
                         levels = c("T1", "T2", "Ta", "T3", "Tb", "T4", "T5",
                                    "T6", "T7", "T8", "T9", "T10")))

# 调整 alpha 多样性值：每个人的时间点减去最早时间点 ------------------------
plotdata <- plotdata %>%
  group_by(Source, people) %>%
  arrange(time12) %>%
  mutate(earliest_time = first(time12[complete.cases(Shannon.wiener)]),
         baseline_shannon = Shannon.wiener[time12 == earliest_time & complete.cases(Shannon.wiener)][1],
         adjusted_shannon = if_else(complete.cases(Shannon.wiener),
                                    Shannon.wiener - baseline_shannon,
                                    NA_real_),
         adjusted_shannon = if_else(time12 == earliest_time & complete.cases(Shannon.wiener),
                                    0,
                                    adjusted_shannon)) %>%
  ungroup() %>%
  select(-baseline_shannon, -earliest_time)

# 计算每个时间点的调整后均值和altitude均值 -----------------------------
mean_data <- plotdata %>%
  group_by(Source, time12) %>%
  summarise(mean_shannon = mean(adjusted_shannon, na.rm = TRUE),
            mean_altitude = mean(altitude, na.rm = TRUE)) %>%
  ungroup()

# 绘图函数（修改版） ------------------------------------------------
make_plot <- function(src) {
  sub <- filter(plotdata, Source == src)
  if (nrow(sub) == 0) return(NULL)
  
  mean_sub <- filter(mean_data, Source == src)
  shannon_range <- range(sub$adjusted_shannon, na.rm = TRUE)
  
  # 设置右侧 y 轴范围
  altitude_min <- 0
  altitude_max <- 20000
  shannon_min <- min(shannon_range[1], na.rm = TRUE)  # 左侧 y 轴最小值根据数据自动设置
  shannon_span <- shannon_range[2] - shannon_min
  altitude_span <- altitude_max - altitude_min
  
  if (altitude_span > 0) {
    scale_factor <- shannon_span / altitude_span
    offset <- shannon_min - altitude_min * scale_factor
  } else {
    scale_factor <- 1
    offset <- 0
  }
  
  # 计算 Pearson 相关性
  cor_value <- cor(mean_sub$mean_shannon, mean_sub$mean_altitude * scale_factor + offset, 
                   method = "pearson", use = "complete.obs")
  cor_label <- sprintf("Pearson r = %.2f", cor_value)
  
  p <- ggplot(sub, aes(x = time12, y = adjusted_shannon, group = people)) +
    geom_smooth(aes(group = people), method = "loess", color = "grey",
                size = 0.5, se = TRUE, alpha = 0.1) +
    geom_line(data = mean_sub, aes(x = time12, y = mean_shannon, group = 1),
              color = "#5CB85C", size = 2) +
    geom_line(data = mean_sub, aes(x = time12, y = mean_altitude * scale_factor + offset, group = 1),
              color = "#357EBD", size = 2) +
    scale_y_continuous(
      name = "Adjusted Shannon Diversity",
      limits = c(shannon_min, shannon_range[2]),  # 左侧 y 轴最小值自动设置
      sec.axis = sec_axis(~ (. - offset) / scale_factor,
                          name = "Altitude (m)", 
                          breaks = seq(0, 20000, by = 2500))
    ) +
    theme_publish() +
    labs(title = paste("Adjusted Alpha Diversity -", src),
         x = "Time Point") +
    annotate("text", x = Inf, y = Inf, label = cor_label, hjust = 1.1, vjust = 1.1, size = 8) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y.right = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "blue"))
  
  return(p)
}

# 生成并保存图 ----------------------------------------------------------
# source_levels <- unique(plotdata$Source)
source_levels <- c("oral", "gut", "skin")
plot_list <- keep(source_levels, ~ !is.null(make_plot(.x))) %>% lapply(make_plot)

if (length(plot_list) > 0) {
  p_all <- grid.arrange(grobs = plot_list, ncol = length(plot_list))
  
  pdf("02.绘制alpha多样性的个体变化.pdf",
      width = 10 * length(plot_list), height = 10)
  grid.draw(p_all)
  dev.off()
  cat("折线图已保存！\n")
} else {
  cat("警告：没有足够数据绘图\n")
}

cat("多样性数据已保存至:", file.path(out_dir, "01.diversity_time12.csv"), "\n")