# Introduction
此脚本使用sva包combat函数以去除批次效应

# 1. Library packages ----
library(readr)
library(readxl)
library(RCurl)
library(genefilter)
library(sva)
library(openxlsx)
library(dplyr)
library(edgeR)
source("./02_code/QC_boxplot.R")
source("./02_code/QC_heatmap.R")
source("./02_code/QC_PCA.R")
# 加载过程中显示版本过老可能不兼容不影响使用，只要不报错就行

# 2. Data input ----
## 2.1 expr input ----
data_batch1 <- read_csv("./01_data/Counts_data/gene_count_matrix(fuben).csv")
data_batch2 <- read_csv("./01_data/Counts_data/gene_count_matrix_1.csv")

# process data
colnames(data_batch1) <- gsub("-","_",colnames(data_batch1))
colnames(data_batch2) <- gsub("-","_",colnames(data_batch2))
data_merge <- merge(data_batch1,data_batch2,by = "gene_id")
data_merge <- as.data.frame(data_merge)
rownames(data_merge) <- data_merge$gene_id
data_merge <- subset(data_merge,select = -c(gene_id))
colnames(data_merge) <- gsub("cas9","",colnames(data_merge))      # 去除列名多余信息
colnames(data_merge) <- gsub("w","W",colnames(data_merge))

# cpm
cpm_data <- cpm(data_merge)
cpm_data <- as.data.frame(cpm_data)
combat_data

# filter
# approach 2
cpm_data_log2 <- log2(cpm_data+1)
filtered_data <- combat_data_log2[apply(combat_data_log2,1,mean)>1,] # 13387
filtered_data <- 2^filtered_data-1

## 2.2 group input ----
data_group_all <- read.xlsx("./01_data/Counts_data/Group_by_IC50.xlsx")
colnames(data_group_all)
data_group_1 <- subset(data_group_all,select = c(id,group))
colnames(data_group_1)[2] <- "group"
data_batch <- subset(data_group_all,select = c(id,batch))
colnames(data_batch)[2] <- "group"
data_condition <- subset(data_group_all,select = c(id,Condition))
colnames(data_condition)[2] <- "group"
table(data_group_all$group)

# 配色设置
value_colour_group <- c("High" = "#E64B35FF",# Experimental group
                        "Medium" = "#4DBBD5FA",# other group1
                        "Low" = "#F2A200"# other group2
)

value_colour_batch <- c("Batch1" = "#E64B35FF",# Experimental group
                        "Batch2" = "#4DBBD5FA"# other group1
)
value_colour_condition <- c("2W" = "#E64B35FF",# Experimental group
                            "6W" = "#4DBBD5FA",# other group1
                            "WT" = "#F2A200"
                            )
# 3. Bath detect ----
## 3.1 set output path ----

dir_qc <- "./03_result/QC/CPM/All/"
data_qc <- filtered_data

## 3.2 boxplot -----------------------------------------------------------------
# 函数有log2（x+1）！！！！！！！！！
pdf(file = paste0(dir_qc,"QC_boxplot.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_qc,data_group = data_group_1,
           value_colour = value_colour_group,
           title = NULL)
dev.off()

## 3.3 heatmap -----------------------------------------------------------------
# 进行heatmap图绘制时，需要去除标准差为零的基因，否则会出现缺失值
# 去除标准差不等于0后的结果赋值为data_before_1
data_qc_sd<-data_qc[apply(data_qc,1,sd)!=0,]

# group
# 函数有log2（x+1）！！！！！！！！！
pdf(file = paste0(dir_qc,"QC_heatmap_group.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_group_1,
           value_colour = value_colour_group)
dev.off()

# batch
pdf(file = paste0(dir_qc,"QC_heatmap_batch.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_batch,
           value_colour = value_colour_batch)
dev.off()

# condition
pdf(file = paste0(dir_qc,"QC_heatmap_condition.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_condition,
           value_colour = value_colour_condition)
dev.off()
## 3.4 PCA ---------------------------------------------------------------------
# 函数没有log2 ！！！！！！！！！
# group
pdf(file = paste0(dir_qc,"QC_pca_group.pdf"),
    width = 15,
    height = 15)
QC_PCA(data = log2(data_qc+1),
       data_group = data_group_1,
       value_colour = value_colour_group)
dev.off()

# batch
pdf(file = paste0(dir_qc,"QC_pca_batch.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = log2(data_qc+1),
       data_group = data_batch,
       value_colour = value_colour_batch)
dev.off()

# condition
pdf(file = paste0(dir_qc,"QC_pca_condition.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = log2(data_qc+1),
       data_group = data_condition,
       value_colour = value_colour_condition)
dev.off()

## 3.5 Cor ---------------------------------------------------------------------
library(corrplot)
# cor analysis
dir_cor <- "./03_result/QC/CPM_batch/Cor/"

  corr <- cor(filtered_data, method = 'pearson')      # cor函数计算两两样本（列与列）之间的相关系数
View(corr)	                                  # 查看样本之间的相关系数

pdf(paste0(dir_cor,'MV4_11_sample_correlation.pdf'), width = 10, height = 10)	
  corrplot(corr, type = 'upper',   # type='upper'：只显示右上角相关系数矩阵
           tl.col = 'black',       # tl.col='black'：字体颜色黑色
           order = 'hclust',       # order='hclust'：使用层次聚类算法
           tl.srt = 45,            # tl.srt = 45：x轴标签倾斜45度
           number.cex = 0.8,        # 相关性系数字体大小
           addCoef.col = 'white')	 # addCoef.col='white'：添加相关系数数值，颜色白色
dev.off() 

# 4. Sva::ComBat ----
## 4.1 design ----
data_group_all$biological <- paste0(data_group_all$Cell,"_",data_group_all$Condition)
data_group_all$biological <- factor(data_group_all$biological, levels = c("MOLM13_2W", "MOLM13_6W", "MOLM13_WT",
                                                                          "MV4_11_2W", "MV4_11_6W", "MV4_11_WT",

                                                                                                                                                    "OCI_AML2_2W", "OCI_AML2_6W", "OCI_AML2_WT"))
design <- model.matrix(~1 + as.factor(group), data=data_group_all) # 0 + 不选择截距项
rownames(design) <- data_group_all$id

## 4.2 adjust batch ----
combat_data <- ComBat(dat = cpm_data, 
                      batch = data_group_all$batch, 
                      mod = design, 
                      par.prior = TRUE)
                      #ref.batch = "Batch1")
combat_data <- as.data.frame(combat_data)
combat_data[combat_data < 0] <- 0
## 4.3 res output ----
write.csv(combat_data,file = "./01_data/Cpm_data/IC50_merged_adjusted.csv")


