sva包去除批次效应

# 加载包 ----
install.packages("BiocManager")
BiocManager::install("BiocParallel")
BiocManager::install("edgeR", version = "3.18")
install.packages("nlme")
BiocManager::install("sva")
library(readr)
library(readxl)
library(RCurl)
library(genefilter)
library(sva)
library(dplyr)
library(edgeR)
source("./02_code/QC_boxplot.R")
source("./02_code/QC_heatmap.R")
source("./02_code/QC_PCA.R")
# 加载过程中显示版本过老可能不兼容不影响使用，只要不报错就行

# input ----
## gene expression input ----
data_batch1 <- read_csv("./01_data/Counts_data/gene_count_matrix(fuben).csv")
data_batch2 <- read_csv("./01_data/Counts_data/gene_count_matrix_1.csv")
colnames(data_batch1) <- gsub("-","_",colnames(data_batch1))
colnames(data_batch2) <- gsub("-","_",colnames(data_batch2))
data_merge <- merge(data_batch1,data_batch2,by = "gene_id")
data_merge <- as.data.frame(data_merge)
rownames(data_merge) <- data_merge$gene_id
data_merge <- subset(data_merge,select = -c(gene_id))
data_merge <- data_merge[,-grep("4w", colnames(data_merge))]

## cpm ----
cpm_data <- cpm(data_merge)
cpm_data <- as.data.frame(cpm_data)
## filter ----
# approach 2
cpm_data_log2 <- log2(cpm_data+1)
filtered_data <- cpm_data_log2[apply(cpm_data_log2,1,mean)>1,] # 13454
filtered_data <- 2^filtered_data-1

## group input ----
data_group_all <- read.csv("./01_data/Sample_info/data_group_annotation.csv")
data_group_all <- data_group_all[-grep("4w",data_group_all$X),] # 删除4w样本
data_batch$id <- gsub("cas9","",data_batch$id)
colnames(data_group_all)
data_group_1 <- subset(data_group_all,select = c(id,Cell))
colnames(data_group_1)[2] <- "group"
data_batch <- subset(data_group_all,select = c(id,Batch))
colnames(data_batch)[2] <- "group"
data_condition <- subset(data_group_all,select = c(id,Condition))
colnames(data_condition)[2] <- "group"


# 配色设置 ----
# 配色设置
value_colour_group <- c("MOLM13" = "#E64B35FF",# Experimental group
                        "MV4_11" = "#4DBBD5FA",# other group1
                        "OCI_AML2" = "#F2A200"# other group2
)

value_colour_batch <- c("Batch1" = "#E64B35FF",# Experimental group
                        "Batch2" = "#4DBBD5FA"# other group1
)
value_colour_condition <- c("2W" = "#E64B35FF",# Experimental group
                            "6W" = "#4DBBD5FA",# other group1
                            "WT" = "#F2A200"
                            )

# 设置输出目录----
#dir.create("03_result")
dir <- "./03_result/QC/Count_batch/"

data_qc <- filtered_data

# boxplot -----------------------------------------------------------------
# 函数有log2（x+1）！！！！！！！！！
pdf(file = paste0(dir,"QC_boxplot.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_qc,data_group = data_group_1,
           value_colour = value_colour_group,
           title = " original_data")
dev.off()

# heatmap -----------------------------------------------------------------
# 进行heatmap图绘制时，需要去除标准差为零的基因，否则会出现缺失值
# 去除标准差不等于0后的结果赋值为data_before_1
data_qc_sd<-data_qc[apply(data_qc,1,sd)!=0,]

# group分组
# 函数有log2（x+1）！！！！！！！！！
pdf(file = paste0(dir,"QC_heatmap_group.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_group_1,
           value_colour = value_colour_group)
dev.off()

# batch分组

pdf(file = paste0(dir,"QC_heatmap_batch.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_batch,
           value_colour = value_colour_batch)
dev.off()

# condition分组

pdf(file = paste0(dir,"QC_heatmap_condition.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_qc_sd,data_group = data_condition,
           value_colour = value_colour_condition)
dev.off()
# PCA ---------------------------------------------------------------------
# 函数没有log2 ！！！！！！！！！
# group
pdf(file = paste0(dir,"QC_pca_group.pdf"),
    width = 15,
    height = 15)
QC_PCA(data = log2(data_qc+1),
       data_group = data_group_1,
       value_colour = value_colour_group)
dev.off()

# batch
pdf(file = paste0(dir,"QC_pca_batch.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = log2(data_qc+1),
       data_group = data_batch,
       value_colour = value_colour_batch)
dev.off()

# condition
pdf(file = paste0(dir,"QC_pca_condition.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = log2(data_qc+1),
       data_group = data_condition,
       value_colour = value_colour_condition)
dev.off()

# Cor ---------------------------------------------------------------------
library(corrplot)
# cor analysis

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

# Combat ----
# 实验设计矩阵（cell line和treatment） ----
data_group_all$biological <- paste0(data_group_all$Cell,"_",data_group_all$Condition)
data_group_all$biological <- factor(data_group_all$biological, levels = c("MOLM13_2W", "MOLM13_6W", "MOLM13_WT",
                                                                          "MV4_11_2W", "MV4_11_6W", "MV4_11_WT",
                                                                          "OCI_AML2_2W", "OCI_AML2_6W", "OCI_AML2_WT"))
design <- model.matrix(~as.factor(biological), data=data_group_all)
rownames(design) <- data_group_all$id

## 矫正 ----
combat_data <- ComBat(dat = cpm_data, 
                      batch = data_group_all$Batch, 
                      mod = design, 
                      par.prior = TRUE)
                      #ref.batch = "Batch1")
combat_data <- as.data.frame(combat_data)
combat_data[combat_data < 0] <- 0
write.csv(combat_data,file = "./01_data/Cpm_data/data_merged_adjusted.csv")


