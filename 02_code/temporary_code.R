#计算转录和蛋白之间的相关性 ----
folder_path <- "../MultiOmics/Corr/"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}
dir_cor <- "../MultiOmics/Corr/"
# 蛋白矩阵 
intensity <- read.csv("../Proteome/01_Data/report.pg_matrix_fill_norma.csv", row.names = 1)
intensity <- intensity[,-23]
anno <- read.xlsx("../Proteome/01_Data/data_anno.xlsx")
anno <- anno[anno$Protein.Group%in%rownames(intensity),]
intensity <- intensity[,order(colnames(intensity))]
identical(anno$Protein.Group, rownames(intensity))
intensity$Genes <- anno$Genes
y <- intensity$Genes
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
intensity$gene <- gene1
a <- intensity[duplicated(intensity$gene),]
# 去重
library(dplyr)
data_clean <- intensity %>%
  drop_na() %>%                                      # 去除含 NA 的行
  mutate(mean_expr = rowMeans(across(1:26))) %>%    # 1~26 列算均值
  group_by(gene) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() 

# 把 gene 作为行名
intensity_matrix <- as.data.frame(data_clean)
rownames(intensity_matrix) <- intensity_matrix$gene
intensity_matrix <- intensity_matrix[, !colnames(intensity_matrix) %in% c("gene", "Genes", "mean_expr")]

# 转录矩阵
tpm <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)
tpm <- tpm[,-21]
y <- rownames(tpm)
gene1 <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
tpm$gene <- gene1
a <- tpm[duplicated(tpm$gene),]
# 去重
library(dplyr)
data_clean <- tpm %>%
  drop_na() %>%                                      # 去除含 NA 的行
  mutate(mean_expr = rowMeans(across(1:26))) %>%    # 1~26 列算均值
  group_by(gene) %>%
  slice_max(order_by = mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() 
# 把 gene 作为行名
tpm_matrix <- as.data.frame(data_clean)
rownames(tpm_matrix) <- tpm_matrix$gene
tpm_matrix <- tpm_matrix[, !colnames(tpm_matrix) %in% c("gene", "Genes", "mean_expr")]

# 转录和蛋白取交集
common_genes <- intersect(rownames(tpm_matrix), rownames(intensity_matrix))
tpm_matrix <- tpm_matrix[common_genes, , drop = FALSE]
intensity_matrix <- intensity_matrix[common_genes, , drop = FALSE]
identical(rownames(tpm_matrix), rownames(intensity_matrix))          # 检查行名是否相同, 包括顺序
# log2
tpm_matrix <- log2(tpm_matrix+1)
intensity_matrix <- log2(intensity_matrix+1)
# 散点图
df <- data.frame(
  x = as.vector(as.matrix(clean_tpm1_1)),
  y = as.vector(as.matrix(clean_tpm1_2))
)

ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", col = "red") +
  xlab("TPM Run1") +
  ylab("TPM Run2") +
  theme_bw()

# cor analysis
corr_matrix <- cor(intensity_matrix, tpm_matrix, method = "spearman")
View(corr_matrix)	                                  # 查看样本之间的相关系数
# 创建显示矩阵，只保留对角线
display_mat <- matrix("", nrow = nrow(corr_matrix), ncol = ncol(corr_matrix))
diag(display_mat) <- sprintf("%.2f", diag(corr_matrix))  # 对角线显示相关系数

# 热图
library(pheatmap)
p <- pheatmap(corr_matrix,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              display_numbers = display_mat,             # 显示数字
              number_format = "%.2f",             # 数字格式：保留两位小数
              number_color = "black",             # 数字颜色
              fontsize_number = 7)               # 数字大小
# 导出为图片
pdf("../MultiOmics/Corr/T&P_correlation_heatmap.pdf", width = 6.5, height = 5.5)
print(p)  # 打印对象到 pdf
dev.off()

# 只计算对应样本之间的相关性
sample_cor <- sapply(1:ncol(clean_tpm1), function(i){
  cor(clean_tpm1[,i], clean_tpm2[,i], method = "pearson")
})
names(sample_cor) <- colnames(clean_tpm1)
view(sample_cor)
# two visualization methods
# methods 1
if(T){
  cairo_pdf(paste0(dir_cor,'QC_corrplot.pdf'), width = 11, height = 11)	
  breaks <- seq(0.5, 1, length.out = 100)  # 强制颜色范围放大 0.9~1
  color_palette <- colorRampPalette(c("white", "#2166AC"))(99)
  corrplot(corr_matrix, type = 'upper',   # type='upper'：只显示右上角
           method = "color",      # ("circle", "square", "ellipse", "number", "shade", "color")
           col = color_palette,
           col.lim = c(0, 1),
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
