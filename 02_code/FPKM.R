# 在线获取基因长度来计算FPKM
#test
# test2
# 加载包 ----
library(tidyverse)
library(openxlsx)
library(dplyr)
library(biomaRt)
library(readr)

# 设置工作路径 ----
setwd('.') 

# 设置输出目录 ----
dir.create("./01_data/Fpkm_data/")
dir <- "./01_data/Fpkm_data/"
# 下载基因长度 ----

# 连接到 Ensembl 数据库
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# 获取基因的起始和结束位置
gene_lengths <- getBM(attributes = c("ensembl_gene_id", 
                                     "start_position", 
                                     "end_position"),
                      mart = ensembl)

# 计算基因长度 (结束位置 - 起始位置 + 1)
gene_lengths$length <- gene_lengths$end_position - gene_lengths$start_position + 1

# 查看结果
head(gene_lengths)

## 读入counts数据 ----
expr_counts <- read.csv("./01_data/Counts_data/data_merge_adjusted.csv")

# 提取ensemble名
y <- expr_counts$X
ensemble <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][1]))
id <- data.frame(ensemble)

# 去除版本号(即去除小数点后的内容，ENSG00000005194.15)
id$ensemble <- gsub("\\.\\d+$", "", id$ensemble)  # 去除版本号

## 合并基因长度数据 ----
expr_counts <- cbind(id,expr_counts)
expr_counts <- merge(expr_counts, gene_lengths, by.x = "ensemble", 
                     by.y = "ensembl_gene_id", all.x = TRUE)

# 查看合并后的数据
head(expr_counts)

# 使用基因长度平均值填补缺失值
average_length <- mean(gene_lengths$length, na.rm = TRUE)
expr_counts$length[is.na(expr_counts$length)] <- average_length

# 计算每个样本的文库大小（总 read count）
library_sizes <- colSums(expr_counts[, 3:29], na.rm = TRUE)

# 初始化一个空的数据框来存储 FPKM 值
fpkm_data <- expr_counts[, 1:2] # 保留 ensemble ID 和原始 counts 数据的前两列

# 遍历每个样本，计算 FPKM
for (i in 3:29) {
  sample_name <- colnames(expr_counts)[i]
  fpkm_values <- (expr_counts[, i] / (expr_counts$length / 10^3)) / (library_sizes[i - 2] / 10^6)
  fpkm_data[[sample_name]] <- fpkm_values
}

# 查看 FPKM 数据
head(fpkm_data)
fpkm_data <- fpkm_data[,-1]

# 保存结果
write.csv(fpkm_data,file = paste0(dir ,"data_merge_adjusted_fpkm.csv"))
          
