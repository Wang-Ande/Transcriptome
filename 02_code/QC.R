# 1. Function&packages ----
source("./02_code/QC_boxplot.R")
source("./02_code/QC_heatmap.R")
source("./02_code/QC_PCA.R")
library(openxlsx)
library(readr)
library(readxl)
library(edgeR)
library(corrplot)

# 2. Data input ----
## 2.1 Group input ----
data_group <- read_excel("./01_data/group_info.xlsx")
data_group <- as.data.frame(data_group)
data_group <- data_group[-21,]
data_group <- data_group[grep("OCI",data_group$id),-3]
colnames(data_group)[2] <- "group"
table(data_group$group)
rownames(data_group) <- data_group$id

## 2.2 Expr input ----
exprSet_raw <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)
exprSet_raw <- as.data.frame(exprSet_raw)
exprSet_raw <- exprSet_raw[,data_group$id]
## 2.3 Filter low expr  ----
# approach 1
min_sam <- 0.5*ncol(exprSet_raw)   # 一半样本中表达
filtered_data <- exprSet_raw[rowSums(exprSet_raw > 1) >= min_sam, ] # 14334

# approach 2 (需要画基因频率直方图)
exprSet_raw_log2 <- log2(exprSet_raw+1)
filtered_data <- exprSet_raw_log2[apply(exprSet_raw_log2,1,mean)>5,] 
filtered_data <- 2^filtered_data-1

# 3. QC --------------------------------------------------------------------------
# 设置输出目录 
folder_path <- "./03_result/01_QC/TPM/OCI/"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

dir_qc <- "./03_result/01_QC/TPM/OCI/"
## 3.0 set group ----
qc_data <- filtered_data
targeted_group <- data_group
table(targeted_group$group) 
# 配色设置
value_colour <- c("High" = "#E64B35FF",
                  "Low" = "#3C5488FF",
                  "WT" = "#F2A200")
                  #"OE_ZBTB7A" = "#3a9491" 
  
## 3.1 Boxplot -----------------------------------------------------------------
# 函数中有log2（x+1）
if(T){
pdf(file = paste0(dir_qc,"QC_boxplot.pdf"),
    width = 6,
    height = 4)
QC_boxplot(qc_data,
           data_group = targeted_group,
           value_colour = value_colour,
           title = NULL)
dev.off()
}
## 3.2 Heatmap -----------------------------------------------------------------
qc_data_sd<-qc_data[apply(qc_data,1,sd)!=0,]
# 函数中进行了log2（x+1）
if(T){
pdf(file = paste0(dir_qc,"QC_heatmap.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = qc_data,
           data_group = targeted_group,
           value_colour = value_colour)
dev.off()
}
## 3.3 PCA ---------------------------------------------------------------------
# 函数中没有log2
if(T){
pdf(file = paste0(dir_qc,"QC_pca.pdf"),
    width = 6,
    height = 6)
QC_PCA(data = log2(qc_data+1),
       data_group = targeted_group,
       value_colour = value_colour)
dev.off()
}
## 3.4 Cor_heatmap -------------------------------------------------------------
#计算样本之间的相关性
dir_cor <- "./03_result/01_QC/TPM/"
gene_expr <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)

# Filter low genes 
filtered_data <- gene_expr[apply(gene_expr,1,
                                   function(x)sum(x>1)>0.5*ncol(gene_expr)),] 
# control 
# data_run <- filtered_data[,grep("WT",colnames(filtered_data))]
sorted_colnames <- sort(colnames(filtered_data))             # 对列名进行排序
data_run <- filtered_data[, sorted_colnames]                 # 按照排序后的列名重新排列

# log2(expr+1)
data_run <-log2(data_run+1)  

# cor analysis
corr_matrix <- cor(data_run, method = 'pearson')      # cor函数计算两两样本（列与列）之间的相关系数
View(corr_matrix)	                                  # 查看样本之间的相关系数

# two visualization methods
# methods 1
if(T){
  cairo_pdf(paste0(dir_cor,'QC_corrplot.pdf'), width = 11, height = 11)	
  breaks <- seq(0.5, 1, length.out = 100)  # 强制颜色范围放大 0.9~1
  color_palette <- colorRampPalette(c("white", "#2166AC"))(99)
  corrplot(corr_matrix, type = 'upper',   # type='upper'：只显示右上角
           method = "color",      # ("circle", "square", "ellipse", "number", "shade", "color")
           col = color_palette,
           col.lim = c(0.78, 1),
           tl.col = 'black',       # tl.col='black'：字体颜色黑色
           order = 'hclust',       # order='hclust'：使用层次聚类算法
           tl.srt = 45,            # tl.srt = 45：x轴标签倾斜45度
           number.cex = 0.75,       # 相关性系数字体大小
           addCoef.col = 'white')	 # addCoef.col='white'：添加相关系数数值，颜色白色
  dev.off() 
}

# methods 2
# 转换相关性矩阵为长格式
library(reshape2)
library(ggplot2)
corr_long <- melt(corr_matrix)

# 使用层次聚类重新排序列名
hc <- hclust(dist(corr_matrix))
ordered_names <- rownames(corr_matrix)[hc$order]

# 重新调整数据框的因子顺序
corr_long$Var1 <- factor(corr_long$Var1, levels = ordered_names)
corr_long$Var2 <- factor(corr_long$Var2, levels = ordered_names)
min(corr_long$value)
# plot
p1 <- ggplot(corr_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#D1E5F0", "#2166AC"), 
                       limits = c(min(corr_matrix), 1),                  # 设定颜色映射范围
                       name = expression(R)) +            # 更改图例标题为 R²
  geom_text(aes(label = sprintf("%.2f", value)), size = 3, color = "white") +   # 显示相关系数，保留2位小数
  labs(title = "Pearson correlation between samples") +     # 添加标题
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签
        axis.title = element_blank(),                       # 隐藏 X 和 Y 轴标题
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中、加粗
        legend.position = "right")                          # 保持图例在右侧
print(p1)
ggsave(paste0(dir_cor,'QC_corrplot.pdf'), width = 9, height = 8)

