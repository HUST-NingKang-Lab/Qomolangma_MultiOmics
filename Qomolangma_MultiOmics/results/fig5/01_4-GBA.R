library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)

selected_metabolite <- "4-Guanidinobutanoic acid"
pheno_vars <- c("Hemoglobin(HGB)", "HCT", "Red Cell Distribution Width - Standard Deviation(RDW-SD)")
met_data_path <- "data/processed/Metabolin_plasma.csv"
pheno_data_path <- "/data/processed/physiology.csv"
OUTPUT_DIR <- getwd()

# 创建输出目录
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

met_data <- read.csv(met_data_path, check.names = FALSE)
pheno_data <- read.csv(pheno_data_path, check.names = FALSE)

# ================================
# step1: 统一处理 time 列
# ================================
met_data <- met_data %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time),
         subject_id = gsub("^S0([0-9])$", "S\\1", subject_id)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

pheno_data <- pheno_data %>%
  mutate(time = as.character(time),
         time = gsub("^T", "", time),
         time = as.integer(time)) %>%
  filter(subject_id %in% paste0("S", 1:13)) %>%
  filter(!is.na(time))

# Add age_group
met_data <- met_data %>%
  mutate(age_group = ifelse(subject_id %in% c("S1", "S10"), "Elderly", "Young"))

pheno_data <- pheno_data %>%
  mutate(age_group = ifelse(subject_id %in% c("S1", "S10"), "Elderly", "Young"))

# ================================
# (1) 代谢物随时间变化的箱线图（分年龄组）
# ================================
cat("\n=== 绘制代谢物箱线图 ===\n")

if (selected_metabolite %in% colnames(met_data)) {
  
  # 提取目标代谢物数据
  plot_data <- met_data %>%
    select(subject_id, time, age_group, !!selected_metabolite) %>%
    rename(metabolite_value = !!selected_metabolite) %>%
    filter(!is.na(metabolite_value))
  
  # 检查每个时间点每个年龄组的样本量
  sample_counts <- plot_data %>%
    group_by(time, age_group) %>%
    summarise(n = n(), .groups = "drop")
  
  cat("样本量统计:\n")
  print(sample_counts)
  
  # 计算中位数用于连线
  median_data <- plot_data %>%
    group_by(time, age_group) %>%
    summarise(median_value = median(metabolite_value, na.rm = TRUE), .groups = "drop")
  
  # 绘图
  p1 <- ggplot(plot_data, aes(x = factor(time), y = metabolite_value, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_line(data = median_data, 
              aes(x = factor(time), y = median_value, group = age_group, color = age_group),
              linewidth = 1) +
    geom_point(data = median_data, 
               aes(x = factor(time), y = median_value, color = age_group),
               size = 3, shape = 18) +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    scale_color_manual(values = c("Elderly" = "#D55E00", "Young" = "#0072B2")) +
    labs(title = paste("Temporal Changes of", selected_metabolite),
         x = "Time Point",
         y = "Metabolite Level",
         fill = "Age Group",
         color = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.position = "top")
  
  # 尝试添加统计检验
  tryCatch({
    # 只在每个时间点都有两组数据时进行检验
    time_points_with_both_groups <- sample_counts %>%
      group_by(time) %>%
      filter(n() == 2, all(n >= 2)) %>%
      pull(time) %>%
      unique()
    
    if (length(time_points_with_both_groups) > 0) {
      # 为每个时间点创建比较对
      comparisons <- lapply(time_points_with_both_groups, function(t) {
        list(c("Elderly", "Young"))
      })
      
      p1 <- p1 + stat_compare_means(
        aes(group = age_group),
        method = "t.test",
        label = "p.format",
        label.y.npc = 0.95,
        size = 3.5
      )
      cat("已添加统计检验\n")
    } else {
      cat("样本量不足，跳过统计检验\n")
    }
  }, error = function(e) {
    cat("统计检验失败:", e$message, "\n")
  })
  
  # 保存图像
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", selected_metabolite)
  ggsave(filename = file.path(OUTPUT_DIR, paste0("metabolite_", safe_name, ".pdf")),
         plot = p1, width = 10, height = 6)
  ggsave(filename = file.path(OUTPUT_DIR, paste0("metabolite_", safe_name, ".png")),
         plot = p1, width = 10, height = 6, dpi = 300)
  
  cat("代谢物箱线图已保存\n")
  
} else {
  cat("警告: 代谢物", selected_metabolite, "不存在于数据中\n")
}

# ================================
# (2) 每个生理表型随时间变化的箱线图（分年龄组）
# ================================
cat("\n=== 绘制表型箱线图 ===\n")

for (pheno_var in pheno_vars) {
  
  cat("\n处理表型:", pheno_var, "\n")
  
  if (!pheno_var %in% colnames(pheno_data)) {
    cat("警告: 表型", pheno_var, "不存在于数据中，跳过\n")
    next
  }
  
  # 提取表型数据
  plot_data_pheno <- pheno_data %>%
    select(subject_id, time, age_group, !!pheno_var) %>%
    rename(pheno_value = !!pheno_var) %>%
    filter(!is.na(pheno_value))
  
  # 检查样本量
  sample_counts_pheno <- plot_data_pheno %>%
    group_by(time, age_group) %>%
    summarise(n = n(), .groups = "drop")
  
  cat("样本量统计:\n")
  print(sample_counts_pheno)
  
  # 样本量过滤：至少有一个时间点有足够的数据
  if (nrow(plot_data_pheno) < 5) {
    cat("样本量不足，跳过该表型\n")
    next
  }
  
  # 计算中位数用于连线
  median_data_pheno <- plot_data_pheno %>%
    group_by(time, age_group) %>%
    summarise(median_value = median(pheno_value, na.rm = TRUE), .groups = "drop")
  
  # 绘图
  p2 <- ggplot(plot_data_pheno, aes(x = factor(time), y = pheno_value, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_line(data = median_data_pheno, 
              aes(x = factor(time), y = median_value, group = age_group, color = age_group),
              linewidth = 1) +
    geom_point(data = median_data_pheno, 
               aes(x = factor(time), y = median_value, color = age_group),
               size = 3, shape = 18) +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    scale_color_manual(values = c("Elderly" = "#D55E00", "Young" = "#0072B2")) +
    labs(title = paste("Temporal Changes of", pheno_var),
         x = "Time Point",
         y = "Phenotype Value",
         fill = "Age Group",
         color = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.position = "top")
  
  # 尝试添加统计检验
  tryCatch({
    time_points_with_both_groups <- sample_counts_pheno %>%
      group_by(time) %>%
      filter(n() == 2, all(n >= 2)) %>%
      pull(time) %>%
      unique()
    
    if (length(time_points_with_both_groups) > 0) {
      p2 <- p2 + stat_compare_means(
        aes(group = age_group),
        method = "t.test",
        label = "p.format",
        label.y.npc = 0.95,
        size = 3.5
      )
      cat("已添加统计检验\n")
    } else {
      cat("样本量不足，跳过统计检验\n")
    }
  }, error = function(e) {
    cat("统计检验失败:", e$message, "\n")
  })
  
  # 保存图像
  safe_pheno_name <- gsub("[^A-Za-z0-9_-]", "_", pheno_var)
  ggsave(filename = file.path(OUTPUT_DIR, paste0("phenotype_", safe_pheno_name, ".pdf")),
         plot = p2, width = 10, height = 6)
  ggsave(filename = file.path(OUTPUT_DIR, paste0("phenotype_", safe_pheno_name, ".png")),
         plot = p2, width = 10, height = 6, dpi = 300)
  
  cat("表型箱线图已保存:", safe_pheno_name, "\n")
}

# ================================
# (3) 代谢物 vs 表型的散点图（整体 Spearman 相关，分年龄组）
# ================================
cat("\n=== 绘制代谢物-表型散点图 ===\n")

if (selected_metabolite %in% colnames(met_data)) {
  
  # 合并代谢数据和表型数据
  merged_data <- met_data %>%
    select(subject_id, time, age_group, !!selected_metabolite) %>%
    rename(metabolite_value = !!selected_metabolite) %>%
    inner_join(pheno_data %>% select(subject_id, time, age_group, all_of(pheno_vars)),
               by = c("subject_id", "time", "age_group"))
  
  # 重塑为长格式
  long_data <- merged_data %>%
    pivot_longer(cols = all_of(pheno_vars),
                 names_to = "phenotype",
                 values_to = "pheno_value") %>%
    filter(!is.na(metabolite_value), !is.na(pheno_value))
  
  if (nrow(long_data) > 0) {
    
    # 为每个表型计算整体 Spearman 相关
    cor_stats <- long_data %>%
      group_by(phenotype) %>%
      summarise(
        cor = cor(metabolite_value, pheno_value, method = "spearman", use = "complete.obs"),
        p_value = cor.test(metabolite_value, pheno_value, method = "spearman", exact = FALSE)$p.value,
        .groups = "drop"
      ) %>%
      mutate(label = sprintf("ρ = %.3f\np = %.3g", cor, p_value))
    
    # 绘图
    p3 <- ggplot(long_data, aes(x = metabolite_value, y = pheno_value, color = age_group)) +
      geom_point(size = 2.5, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
      scale_color_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
      facet_wrap(~ phenotype, scales = "free", ncol = 2) +
      labs(title = paste("Correlation between", selected_metabolite, "and Phenotypes"),
         x = paste(selected_metabolite, "Level"),
         y = "Phenotype Value",
         color = "Age Group") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 11),
            legend.position = "top",
            strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_rect(fill = "grey90"))
    
    # 添加相关系数标签
    p3 <- p3 + geom_text(data = cor_stats,
                         aes(x = -Inf, y = Inf, label = label),
                         hjust = -0.1, vjust = 1.2,
                         size = 3.5, color = "black",
                         inherit.aes = FALSE)
    
    # 保存图像
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", selected_metabolite)
    ggsave(filename = file.path(OUTPUT_DIR, paste0("correlation_", safe_name, "_phenotypes.pdf")),
           plot = p3, width = 12, height = 8)
    ggsave(filename = file.path(OUTPUT_DIR, paste0("correlation_", safe_name, "_phenotypes.png")),
           plot = p3, width = 12, height = 8, dpi = 300)
    
    cat("散点图已保存\n")
    cat("\n相关性统计结果:\n")
    print(cor_stats)
    
  } else {
    cat("警告: 合并后无有效数据\n")
  }
  
} else {
  cat("警告: 代谢物", selected_metabolite, "不存在于数据中\n")
}

cat("\n=== 所有图表生成完成 ===\n")
cat("输出目录:", OUTPUT_DIR, "\n")





# ================================
# 方法一：同步变化率（Synchrony Rate）分析
# ================================
cat("\n=== 方法一：同步变化率分析 ===\n")

# 计算同步变化率的函数
calculate_synchrony <- function(metabolite_vec, phenotype_vec) {
  # 移除NA值
  valid_idx <- !is.na(metabolite_vec) & !is.na(phenotype_vec)
  met <- metabolite_vec[valid_idx]
  phen <- phenotype_vec[valid_idx]
  
  if (length(met) < 2) return(NA)
  
  # 计算变化量
  delta_met <- diff(met)
  delta_phen <- diff(phen)
  
  # 判断同步（同向变化）
  synchrony <- (delta_met > 0 & delta_phen > 0) | (delta_met < 0 & delta_phen < 0)
  
  # 计算同步比例
  sync_rate <- sum(synchrony) / length(synchrony)
  
  return(list(
    sync_rate = sync_rate,
    n_intervals = length(synchrony),
    n_sync = sum(synchrony)
  ))
}

# 为每个受试者计算同步率
if (selected_metabolite %in% colnames(met_data)) {
  
  # 合并数据
  sync_data <- met_data %>%
    select(subject_id, time, age_group, !!selected_metabolite) %>%
    rename(metabolite_value = !!selected_metabolite) %>%
    inner_join(pheno_data %>% select(subject_id, time, age_group, all_of(pheno_vars)),
               by = c("subject_id", "time", "age_group")) %>%
    arrange(subject_id, time)
  
  # 对每个表型分别计算
  synchrony_results_list <- list()
  
  for (pheno_var in pheno_vars) {
    cat("\n处理表型:", pheno_var, "\n")
    
    # 计算每个个体的同步率
    individual_sync <- sync_data %>%
      group_by(subject_id, age_group) %>%
      summarise(
        sync_result = list(calculate_synchrony(metabolite_value, .data[[pheno_var]])),
        .groups = "drop"
      ) %>%
      mutate(
        sync_rate = sapply(sync_result, function(x) ifelse(is.na(x), NA, x$sync_rate)),
        n_intervals = sapply(sync_result, function(x) ifelse(is.na(x), NA, x$n_intervals)),
        n_sync = sapply(sync_result, function(x) ifelse(is.na(x), NA, x$n_sync))
      ) %>%
      select(-sync_result) %>%
      filter(!is.na(sync_rate)) %>%
      mutate(phenotype = pheno_var)
    
    synchrony_results_list[[pheno_var]] <- individual_sync
    
    # 打印统计信息
    cat("\n个体同步率统计:\n")
    print(individual_sync %>% select(subject_id, age_group, sync_rate, n_intervals, n_sync))
  }
  
  # 合并所有表型的结果
  all_synchrony <- bind_rows(synchrony_results_list)
  
  # 计算组间统计
  group_stats <- all_synchrony %>%
    group_by(phenotype, age_group) %>%
    summarise(
      mean_sync = mean(sync_rate, na.rm = TRUE),
      sd_sync = sd(sync_rate, na.rm = TRUE),
      median_sync = median(sync_rate, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  cat("\n组间同步率统计:\n")
  print(group_stats)
  
  # 统计检验
  cat("\n统计检验结果:\n")
  for (pheno_var in pheno_vars) {
    pheno_sync <- all_synchrony %>% filter(phenotype == pheno_var)
    
    elderly_rates <- pheno_sync %>% filter(age_group == "Elderly") %>% pull(sync_rate)
    young_rates <- pheno_sync %>% filter(age_group == "Young") %>% pull(sync_rate)
    
    if (length(elderly_rates) >= 1 && length(young_rates) >= 1) {
      if (length(elderly_rates) >= 2 && length(young_rates) >= 2) {
        test_result <- t.test(young_rates, elderly_rates)
        cat(sprintf("\n%s: Young(%.1f%%) vs Elderly(%.1f%%), p=%.4f\n",
                    pheno_var,
                    mean(young_rates)*100,
                    mean(elderly_rates)*100,
                    test_result$p.value))
      } else {
        cat(sprintf("\n%s: 样本量不足，无法进行t检验\n", pheno_var))
      }
    }
  }
  
  # ===== 可视化1: 箱线图 =====
  p_sync_box <- ggplot(all_synchrony, aes(x = age_group, y = sync_rate * 100, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2.5) +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    facet_wrap(~ phenotype, ncol = 2) +
    labs(title = "Synchrony Rate: Metabolite-Phenotype Coupling",
         subtitle = paste("Metabolite:", selected_metabolite),
         x = "Age Group",
         y = "Synchrony Rate (%)",
         fill = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.position = "top",
          strip.text = element_text(size = 10, face = "bold"))
  
  # 添加统计检验标注
  tryCatch({
    p_sync_box <- p_sync_box + 
      stat_compare_means(aes(group = age_group),
                         method = "t.test",
                         label = "p.format",
                         size = 3.5)
  }, error = function(e) {
    cat("统计标注失败:", e$message, "\n")
  })
  
  ggsave(filename = file.path(OUTPUT_DIR, "synchrony_rate_boxplot.pdf"),
         plot = p_sync_box, width = 10, height = 8)
  ggsave(filename = file.path(OUTPUT_DIR, "synchrony_rate_boxplot.png"),
         plot = p_sync_box, width = 10, height = 8, dpi = 300)
  
  cat("\n同步率箱线图已保存\n")
  
  # ===== 可视化2: 柱状图（平均同步率） =====
  p_sync_bar <- ggplot(group_stats, aes(x = age_group, y = mean_sync * 100, fill = age_group)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_errorbar(aes(ymin = (mean_sync - sd_sync) * 100, 
                      ymax = (mean_sync + sd_sync) * 100),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", mean_sync * 100)),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    facet_wrap(~ phenotype, ncol = 2) +
    labs(title = "Average Synchrony Rate by Age Group",
         subtitle = paste("Metabolite:", selected_metabolite),
         x = "Age Group",
         y = "Mean Synchrony Rate (%)",
         fill = "Age Group") +
    ylim(0, max(group_stats$mean_sync * 100) * 1.2) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.position = "top",
          strip.text = element_text(size = 10, face = "bold"))
  
  ggsave(filename = file.path(OUTPUT_DIR, "synchrony_rate_barplot.pdf"),
         plot = p_sync_bar, width = 10, height = 8)
  ggsave(filename = file.path(OUTPUT_DIR, "synchrony_rate_barplot.png"),
         plot = p_sync_bar, width = 10, height = 8, dpi = 300)
  
  cat("同步率柱状图已保存\n")
  
  # ===== 可视化3: 个体轨迹示例图 =====
  # 选择一个年轻和一个老年个体展示
  example_young <- "S2"
  example_elderly <- "S1"
  
  example_data <- sync_data %>%
    filter(subject_id %in% c(example_young, example_elderly)) %>%
    select(subject_id, time, age_group, metabolite_value, all_of(pheno_vars[1])) %>%
    rename(phenotype_value = !!pheno_vars[1])
  
  if (nrow(example_data) > 0) {
    # 标准化数据用于展示
    example_data_norm <- example_data %>%
      group_by(subject_id) %>%
      mutate(
        metabolite_norm = scale(metabolite_value)[,1],
        phenotype_norm = scale(phenotype_value)[,1]
      ) %>%
      ungroup()
    
    # 转为长格式
    example_long <- example_data_norm %>%
      select(subject_id, time, age_group, metabolite_norm, phenotype_norm) %>%
      pivot_longer(cols = c(metabolite_norm, phenotype_norm),
                   names_to = "variable",
                   values_to = "normalized_value") %>%
      mutate(variable = recode(variable,
                               metabolite_norm = selected_metabolite,
                               phenotype_norm = pheno_vars[1]))
    
    p_trajectory <- ggplot(example_long, aes(x = time, y = normalized_value, 
                                              color = variable, shape = variable)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3.5) +
      facet_wrap(~ paste(subject_id, "-", age_group), ncol = 2, scales = "free_x") +
      scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
      scale_shape_manual(values = c(16, 17)) +
      labs(title = "Example Trajectories: Metabolite-Phenotype Coupling",
           subtitle = "Normalized values showing synchronous changes",
           x = "Time Point",
           y = "Normalized Value (Z-score)",
           color = "Variable",
           shape = "Variable") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 11),
            legend.position = "top",
            strip.text = element_text(size = 10, face = "bold"))
    
    ggsave(filename = file.path(OUTPUT_DIR, "synchrony_trajectory_examples.pdf"),
           plot = p_trajectory, width = 12, height = 6)
    ggsave(filename = file.path(OUTPUT_DIR, "synchrony_trajectory_examples.png"),
           plot = p_trajectory, width = 12, height = 6, dpi = 300)
    
    cat("轨迹示例图已保存\n")
  }
  
  # 保存同步率结果表
  write.csv(all_synchrony, 
            file.path(OUTPUT_DIR, "synchrony_rate_individual.csv"),
            row.names = FALSE)
  write.csv(group_stats,
            file.path(OUTPUT_DIR, "synchrony_rate_summary.csv"),
            row.names = FALSE)
}

# ================================
# 方法二：动态相关系数（Within-Subject Correlation）
# ================================
cat("\n=== 方法二：个体内时间序列相关分析 ===\n")

# 计算个体内相关系数的函数
calculate_within_subject_cor <- function(metabolite_vec, phenotype_vec, method = "spearman") {
  # 移除NA值
  valid_idx <- !is.na(metabolite_vec) & !is.na(phenotype_vec)
  met <- metabolite_vec[valid_idx]
  phen <- phenotype_vec[valid_idx]
  
  if (length(met) < 3) return(list(cor = NA, p_value = NA, n = length(met)))
  
  # 计算相关系数
  test_result <- cor.test(met, phen, method = method, exact = FALSE)
  
  return(list(
    cor = as.numeric(test_result$estimate),
    p_value = test_result$p.value,
    n = length(met)
  ))
}

if (selected_metabolite %in% colnames(met_data)) {
  
  # 合并数据
  cor_data <- met_data %>%
    select(subject_id, time, age_group, !!selected_metabolite) %>%
    rename(metabolite_value = !!selected_metabolite) %>%
    inner_join(pheno_data %>% select(subject_id, time, age_group, all_of(pheno_vars)),
               by = c("subject_id", "time", "age_group")) %>%
    arrange(subject_id, time)
  
  # 对每个表型分别计算
  correlation_results_list <- list()
  
  for (pheno_var in pheno_vars) {
    cat("\n处理表型:", pheno_var, "\n")
    
    # 计算每个个体的相关系数
    individual_cor <- cor_data %>%
      group_by(subject_id, age_group) %>%
      summarise(
        cor_result = list(calculate_within_subject_cor(metabolite_value, .data[[pheno_var]], method = "spearman")),
        .groups = "drop"
      ) %>%
      mutate(
        correlation = sapply(cor_result, function(x) x$cor),
        p_value = sapply(cor_result, function(x) x$p_value),
        n_points = sapply(cor_result, function(x) x$n)
      ) %>%
      select(-cor_result) %>%
      filter(!is.na(correlation)) %>%
      mutate(phenotype = pheno_var)
    
    correlation_results_list[[pheno_var]] <- individual_cor
    
    # 打印统计信息
    cat("\n个体相关系数:\n")
    print(individual_cor %>% select(subject_id, age_group, correlation, p_value, n_points))
  }
  
  # 合并所有表型的结果
  all_correlations <- bind_rows(correlation_results_list)
  
  # 计算组间统计
  cor_group_stats <- all_correlations %>%
    group_by(phenotype, age_group) %>%
    summarise(
      mean_cor = mean(correlation, na.rm = TRUE),
      sd_cor = sd(correlation, na.rm = TRUE),
      median_cor = median(correlation, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  cat("\n组间相关系数统计:\n")
  print(cor_group_stats)
  
  # 统计检验
  cat("\n统计检验结果:\n")
  for (pheno_var in pheno_vars) {
    pheno_cor <- all_correlations %>% filter(phenotype == pheno_var)
    
    elderly_cors <- pheno_cor %>% filter(age_group == "Elderly") %>% pull(correlation)
    young_cors <- pheno_cor %>% filter(age_group == "Young") %>% pull(correlation)
    
    if (length(elderly_cors) >= 1 && length(young_cors) >= 1) {
      if (length(elderly_cors) >= 2 && length(young_cors) >= 2) {
        test_result <- t.test(young_cors, elderly_cors)
        cat(sprintf("\n%s: Young(r=%.3f) vs Elderly(r=%.3f), p=%.4f\n",
                    pheno_var,
                    mean(young_cors),
                    mean(elderly_cors),
                    test_result$p.value))
      } else {
        cat(sprintf("\n%s: 样本量不足，无法进行t检验\n", pheno_var))
      }
    }
  }
  
  # ===== 可视化1: 箱线图（相关系数分布） =====
  p_cor_box <- ggplot(all_correlations, aes(x = age_group, y = correlation, fill = age_group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    facet_wrap(~ phenotype, ncol = 2) +
    labs(title = "Within-Subject Correlation: Metabolite-Phenotype Coupling",
         subtitle = paste("Metabolite:", selected_metabolite, "| Method: Spearman"),
         x = "Age Group",
         y = "Correlation Coefficient (ρ)",
         fill = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.position = "top",
          strip.text = element_text(size = 10, face = "bold"))
  
  # 添加统计检验标注
  tryCatch({
    p_cor_box <- p_cor_box + 
      stat_compare_means(aes(group = age_group),
                         method = "t.test",
                         label = "p.format",
                         size = 3.5)
  }, error = function(e) {
    cat("统计标注失败:", e$message, "\n")
  })
  
  ggsave(filename = file.path(OUTPUT_DIR, "within_subject_correlation_boxplot.pdf"),
         plot = p_cor_box, width = 10, height = 8)
  ggsave(filename = file.path(OUTPUT_DIR, "within_subject_correlation_boxplot.png"),
         plot = p_cor_box, width = 10, height = 8, dpi = 300)
  
  cat("\n相关系数箱线图已保存\n")
  
  # ===== 可视化2: 柱状图（平均相关系数） =====
  p_cor_bar <- ggplot(cor_group_stats, aes(x = age_group, y = mean_cor, fill = age_group)) +
    geom_col(alpha = 0.8, width = 0.6) +
    geom_errorbar(aes(ymin = mean_cor - sd_cor, ymax = mean_cor + sd_cor),
                  width = 0.2, linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", mean_cor)),
              vjust = ifelse(cor_group_stats$mean_cor > 0, -0.5, 1.5), 
              size = 4, fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    scale_fill_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    facet_wrap(~ phenotype, ncol = 2) +
    labs(title = "Average Within-Subject Correlation by Age Group",
         subtitle = paste("Metabolite:", selected_metabolite),
         x = "Age Group",
         y = "Mean Correlation Coefficient (ρ)",
         fill = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.position = "top",
          strip.text = element_text(size = 10, face = "bold"))
  
  ggsave(filename = file.path(OUTPUT_DIR, "within_subject_correlation_barplot.pdf"),
         plot = p_cor_bar, width = 10, height = 8)
  ggsave(filename = file.path(OUTPUT_DIR, "within_subject_correlation_barplot.png"),
         plot = p_cor_bar, width = 10, height = 8, dpi = 300)
  
  cat("相关系数柱状图已保存\n")
  
  # ===== 可视化3: 散点图矩阵（个体相关强度可视化） =====
  # 为每个个体绘制散点图，按相关系数强度排序
  top_subjects <- all_correlations %>%
    filter(phenotype == pheno_vars[1]) %>%
    arrange(desc(abs(correlation))) %>%
    head(6) %>%
    pull(subject_id)
  
  scatter_data <- cor_data %>%
    filter(subject_id %in% top_subjects) %>%
    select(subject_id, time, age_group, metabolite_value, !!pheno_vars[1]) %>%
    rename(phenotype_value = !!pheno_vars[1]) %>%
    left_join(
      all_correlations %>% 
        filter(phenotype == pheno_vars[1]) %>%
        select(subject_id, correlation),
      by = "subject_id"
    ) %>%
    mutate(label = sprintf("%s (%s)\nρ=%.3f", subject_id, age_group, correlation))
  
  if (nrow(scatter_data) > 0) {
    p_scatter_matrix <- ggplot(scatter_data, 
                                aes(x = metabolite_value, y = phenotype_value, color = age_group)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
      scale_color_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
      facet_wrap(~ label, scales = "free", ncol = 3) +
      labs(title = "Individual Coupling Strength Examples",
           subtitle = paste("Top 6 subjects ranked by |correlation|"),
           x = paste(selected_metabolite, "Level"),
           y = pheno_vars[1],
           color = "Age Group") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 11),
            axis.text = element_text(size = 9),
            axis.title = element_text(size = 10),
            legend.position = "top",
            strip.text = element_text(size = 9, face = "bold"))
    
    ggsave(filename = file.path(OUTPUT_DIR, "within_subject_scatter_matrix.pdf"),
           plot = p_scatter_matrix, width = 14, height = 10)
    ggsave(filename = file.path(OUTPUT_DIR, "within_subject_scatter_matrix.png"),
           plot = p_scatter_matrix, width = 14, height = 10, dpi = 300)
    
    cat("散点图矩阵已保存\n")
  }
  
  # 保存相关系数结果表
  write.csv(all_correlations, 
            file.path(OUTPUT_DIR, "within_subject_correlation_individual.csv"),
            row.names = FALSE)
  write.csv(cor_group_stats,
            file.path(OUTPUT_DIR, "within_subject_correlation_summary.csv"),
            row.names = FALSE)
}

# ================================
# 综合对比图
# ================================
cat("\n=== 生成综合对比图 ===\n")

if (exists("all_synchrony") && exists("all_correlations")) {
  
  # 合并两种方法的结果
  combined_results <- all_synchrony %>%
    select(subject_id, age_group, phenotype, sync_rate) %>%
    left_join(
      all_correlations %>% select(subject_id, age_group, phenotype, correlation),
      by = c("subject_id", "age_group", "phenotype")
    )
  
  # 绘制散点图：同步率 vs 相关系数
  p_combined <- ggplot(combined_results, 
                       aes(x = sync_rate * 100, y = correlation, color = age_group, shape = age_group)) +
    geom_point(size = 3.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1) +
    scale_color_manual(values = c("Elderly" = "#E69F00", "Young" = "#56B4E9")) +
    scale_shape_manual(values = c(16, 17)) +
    facet_wrap(~ phenotype, ncol = 2) +
    labs(title = "Comparison of Two Coupling Methods",
         subtitle = "Synchrony Rate vs Within-Subject Correlation",
         x = "Synchrony Rate (%)",
         y = "Correlation Coefficient (ρ)",
         color = "Age Group",
         shape = "Age Group") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 11),
          legend.position = "top",
          strip.text = element_text(size = 10, face = "bold"))
  
  ggsave(filename = file.path(OUTPUT_DIR, "methods_comparison.pdf"),
         plot = p_combined, width = 10, height = 8)
  ggsave(filename = file.path(OUTPUT_DIR, "methods_comparison.png"),
         plot = p_combined, width = 10, height = 8, dpi = 300)
  
  cat("综合对比图已保存\n")
  
  # 保存合并结果
  write.csv(combined_results,
            file.path(OUTPUT_DIR, "coupling_methods_combined.csv"),
            row.names = FALSE)
}

cat("\n=== 耦合分析完成 ===\n")
cat("所有结果已保存至:", OUTPUT_DIR, "\n")
cat("\n生成的图表包括:\n")
cat("1. 同步率箱线图 (synchrony_rate_boxplot)\n")
cat("2. 同步率柱状图 (synchrony_rate_barplot)\n")
cat("3. 轨迹示例图 (synchrony_trajectory_examples)\n")
cat("4. 相关系数箱线图 (within_subject_correlation_boxplot)\n")
cat("5. 相关系数柱状图 (within_subject_correlation_barplot)\n")
cat("6. 散点图矩阵 (within_subject_scatter_matrix)\n")
cat("7. 方法对比图 (methods_comparison)\n")