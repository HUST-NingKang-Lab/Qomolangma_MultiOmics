library(dplyr)
library(stringr)

MAG_stat <- read.csv("data/processed/02.MAG_stat_with_metadata_R226.csv")

MAG_checkM1 <- read.delim("data/processed/all_bins_checkM.tsv", check.names = FALSE)
MAG_checkM2 <- read.delim("data/processed/checkm_bu.tsv", check.names = FALSE)
MAG_checkM <- rbind(MAG_checkM1, MAG_checkM2)
colnames(MAG_checkM) <- c("MAG", "Marker_Lineage", "Marker_genomes", "Marker_genes", "Marker_genes_set", "miss_gene",
     "one_gene", "two_gene", "three_gene", "four_gene", "more_gene", "Completeness", "Contamination", "Strain_heterogeneity")

### 处理MAG的名字 ###
MAG_checkM$MAG <- gsub("-1_bin", "_bin", MAG_checkM$MAG)  # 补运行的16个样本，后缀多了个-1
MAG_checkM$MAG <- gsub("^D", "", MAG_checkM$MAG)  # 补运行的16个样本，前缀多了个D

### 提取SampleID列 ###
MAG_checkM <- MAG_checkM %>% mutate(SampleID = str_split_i(MAG, "_", 1))
MAG_checkM <- MAG_checkM %>% select(MAG, SampleID, everything())

### 给MAG添加后缀.fa ###
MAG_checkM <- MAG_checkM %>% mutate(MAG = paste0(MAG, ".fa"))

## 去除不要的样本
# QML6500-240426-14C, QML6500-240426-14D, QML7028-240502-07D
MAG_checkM <- MAG_checkM %>% filter(!SampleID %in% c("QML6500-240426-14C", "QML6500-240426-14D", "QML7028-240502-07D"))
print(length(unique(MAG_checkM$SampleID)))  # 正确的会有348个样本
print(length(unique(MAG_checkM$MAG)))  # 有12081个MAG
write.csv(MAG_checkM, "data/processed/07.MAG的checkM结果.csv", row.names = FALSE)



MAG_checkM <- read.csv("data/07.MAG_checkM_result.csv", check.names = FALSE)

MAG_checkM_inferior <- MAG_checkM %>% 
    filter(Completeness < 50 | Contamination > 10)

MAG_checkM_non_inferior <- MAG_checkM %>% 
    filter(!(Completeness < 50 | Contamination > 10))

print(length(MAG_checkM_non_inferior$MAG))

write.csv(MAG_checkM_inferior, "08.低质量的MAG的checkM结果.csv", row.names = FALSE)
write.csv(MAG_checkM_non_inferior, "08.非低质量的MAG的checkM结果.csv", row.names = FALSE)