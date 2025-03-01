# 1. Function&packages ----
source("./02_code/QC_boxplot.R")
source("./02_code/QC_heatmap.R")
source("./02_code/QC_PCA.R")
library(openxlsx)
library(readr)
library(readxl)
library(edgeR)
if (!requireNamespace('corrplot', quietly = TRUE))
  install.packages('corrplot')
library(corrplot)

# 2. Data input ----
## 2.1 Group input ----
# 导入分组信息
data_group <- read_excel("./01_data/Counts_data/aml_group_merge.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
data_group <- data_group[-grep("4w",data_group$id),]          # 删除4周样本
colnames(data_group)
data_group$id <- gsub("cas9","",data_group$id)                # 去除id列多余信息
table(data_group$group)
rownames(data_group) <- data_group$id
data_group <- data_group[-grep("MOLM13_WT_1",data_group$id),]

## 2.2 Expr input ----
exprSet_raw <- read.csv("./01_data/Counts_data/data_merge_adjusted.csv")
exprSet_raw <- as.data.frame(exprSet_raw)

# 提取gene_id列作为行名
colnames(exprSet_raw)
colnames(exprSet_raw)[1] <- "gene_id"
rownames(exprSet_raw) <- exprSet_raw$gene_id
exprSet_raw <- exprSet_raw[,-1]

# 去除样本名多余信息
colnames(exprSet_raw) <- gsub("cas9","",colnames(exprSet_raw))  
exprSet_raw <- exprSet_raw[,-grep("4w",colnames(exprSet_raw))] # 删除4周样本
exprSet_raw <- exprSet_raw[,-grep("MOLM13_WT_1",colnames(exprSet_raw))]

# Filter low genes 
# approach 1
min_sam <- 3
filtered_data <- exprSet_raw[apply(exprSet_raw,1,
                                function(x)sum(x>1)>min_sam),] 
# approach 2
exprSet_raw_log2 <- log2(exprSet_raw+1)
filtered_data <- exprSet_raw_log2[apply(exprSet_raw_log2,1,mean)>5,] # 13454
filtered_data <- 2^filtered_data-1

# 3. QC --------------------------------------------------------------------------
# 设置输出目录 
folder_path <- "./03_result/QC/Count/All"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

dir <- "./03_result/QC/CPM/MOLM13/"

qc_data <- filtered_data

## 3.0 set group ----
table(data_group$group)
targeted_group <- data_group

# VEN-WT 组别
targeted_group$group <- paste0(targeted_group$group, "_VEN")
targeted_group[grep("WT",targeted_group$id),2] <- gsub("_VEN","_WT",
                                                       targeted_group[grep("WT",targeted_group$id),2])
table(targeted_group$group)

# 2W/6W-WT 组别
targeted_group[grep("6w",targeted_group$id),2] <- gsub("_VEN","_6W",
                                                       targeted_group[grep("6w",targeted_group$id),2])
table(targeted_group$group)

targeted_group <- targeted_group[grep("MOLM13",targeted_group$id),]
# 配色设置
value_colour <- c("MOLM13" = "#E64B35FF",# Experimental group
                  "MV4_11" = "#4DBBD5FA",# other group1
                  "OCI_AML2" = "#F2A200")  # other group2


## 3.1 Boxplot -----------------------------------------------------------------
# 函数中有log2（x+1）

pdf(file = paste0(dir,"QC_boxplot.pdf"),
    width = 6,
    height = 4)
QC_boxplot(qc_data,
           data_group = targeted_group,
           value_colour = value_colour,
           title = "Normalization Data")
dev.off()

## 3.2 Heatmap -----------------------------------------------------------------
qc_data_1 <- qc_data[,grep("MOLM13",colnames(qc_data))]
qc_data_1_sd<-qc_data[apply(qc_data,1,sd)!=0,]

# 函数中进行了log2（x+1）
if(T){
pdf(file = paste0(dir,"QC_heatmap.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = qc_data_1_sd,
           data_group = targeted_group,
           value_colour = value_colour)
dev.off()
}
## 3.3 PCA ---------------------------------------------------------------------
# 函数中没有log2
pdf(file = paste0(dir,"QC_pca.pdf"),
    width = 7,
    height = 7)
QC_PCA(data = log2(qc_data+1),
       data_group = targeted_group,
       value_colour = value_colour)
dev.off()

## 3.4 Cor_heatmap -------------------------------------------------------------
#计算样本之间的相关性
dir_cor <- "./03_result/QC/CPM/Cor/"
gene_expr <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
rownames(gene_expr) <- gene_expr$X
gene_expr <- gene_expr[,-1]
colnames(gene_expr) <- gsub("cas9","",colnames(gene_expr))# 去除样本名中不必要信息

# Filter low genes 
filtered_data <- gene_expr[apply(gene_expr,1,
                                   function(x)sum(x>1)>0.5*ncol(gene_expr)),] 

# control 
gene_WT <- filtered_data[,grep("WT",colnames(filtered_data))]
sorted_colnames <- sort(colnames(gene_WT))            # 对列名进行排序
gene_WT <- gene_WT[, sorted_colnames]                 # 按照排序后的列名重新排列


# experimental
gene_2W <- filtered_data[,grep("2w",colnames(filtered_data))]
sorted_colnames <- sort(colnames(gene_2W))            # 对列名进行排序
gene_2W <- gene_2W[, sorted_colnames]                 # 按照排序后的列名重新排列

# Cell line
gene_OCI_AML2 <- filtered_data[,grep("OCI_AML2",colnames(filtered_data))]
sorted_colnames <- sort(colnames(gene_OCI_AML2))            # 对列名进行排序
gene_OCI_AML2 <- gene_OCI_AML2[, sorted_colnames]                 # 按照排序后的列名重新排列

# log2(expr+1)
gene_MV4_11 <-log2(gene_MV4_11+1)  

# cor analysis
corr <- cor(gene_MV4_11, method = 'pearson')      # cor函数计算两两样本（列与列）之间的相关系数
View(corr)	                                  # 查看样本之间的相关系数

if(T){
pdf(paste0(dir_cor,'MV4_11_sample_correlation.pdf'), width = 10, height = 10)	
corrplot(corr, type = 'upper',   # type='upper'：只显示右上角相关系数矩阵
         tl.col = 'black',       # tl.col='black'：字体颜色黑色
         order = 'hclust',       # order='hclust'：使用层次聚类算法
         tl.srt = 45,            # tl.srt = 45：x轴标签倾斜45度
         number.cex = 0.8,        # 相关性系数字体大小
         addCoef.col = 'white')	 # addCoef.col='white'：添加相关系数数值，颜色白色
 dev.off() 
}
