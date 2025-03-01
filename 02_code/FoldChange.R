# Foldchange计算
# single sample foldchange
# foldchange定义为对照组和实验组中基因平均表达值的比值，单样本直接计算比值即可
library(openxlsx)
library(readr)
library(ggplot2)
library(ggrepel)
# 一、single sample logfc ----
# 1 data input ----
expr <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
rownames(expr) <- expr$X
expr <- expr[,-1]

expr_anno <- read.xlsx("../Proteome/01_Data/data_anno.xlsx")
rownames(expr_anno) <- expr_anno$Protein.Group
expr_anno <- expr_anno[rownames(expr_anno)%in%rownames(expr),]
expr$gene <- expr_anno$Genes
colnames(expr)
# filter lowexpr genes
expr <- log2(expr+1)
expr_filter <- expr[apply(expr,1,mean)>1,] 

# 去除样本名多余信息
colnames(expr_filter) <- gsub("cas9","",colnames(expr_filter)) 
colnames(expr_filter) <- gsub("w","W",colnames(expr_filter))

# 确保所有的实验组和对照组配平
expr_filter$MOLM13_2W_1 <- 1
expr_filter <- expr_filter[,-grep("4W",colnames(expr_filter))]
colnames(expr_filter)

# 2. subset symbol names ----
y <- rownames(expr_filter)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
expr_filter$gene <- gene
expr_filter <- expr_filter[, order(names(expr_filter))]


# 若进行log2，先转置回原格式
expr_filter <- 2^expr_filter[,2:28]-1
expr_filter <- expr_filter+1
expr <- expr_filter
# 3. creat fc df ----
# 定义细胞系和对应的时间点
cell_lines <- c("MOLM13", "MV4_11", "OCI_AML2")
time_points <- c("WT", "2W", "6W")

# 创建一个空的 fold_change 数据框
fold_change <- data.frame(Gene = gene)

# 计算单样本 fold change
for (cell_line in cell_lines) {
  # 获取对应细胞系的列名
  wt_samples <- grep(paste0(cell_line, "_WT"), colnames(expr), value = TRUE)
  two_w_samples <- grep(paste0(cell_line, "_2W"), colnames(expr), value = TRUE)
  six_w_samples <- grep(paste0(cell_line, "_6W"), colnames(expr), value = TRUE)
  
  # 单样本对照组 (WT) 和实验组 (2W, 6W) 之间的 Fold Change
  for (i in 1:length(two_w_samples)) {
    fold_change[paste0(cell_line, "_2W_", i, "_vs_WT_", i)] <- expr[[two_w_samples[i]]] / expr[[wt_samples[i]]]
  }
  
  for (i in 1:length(six_w_samples)) {
    fold_change[paste0(cell_line, "_6W_", i, "_vs_WT_", i)] <- expr[[six_w_samples[i]]] / expr[[wt_samples[i]]]
  }
}

fold_change[,c(2:19)] <- log2(fold_change[,c(2:19)])
# 4. res output ----
write.xlsx(fold_change,file = "./01_Data/Cpm_data/Foldchange.xlsx")

# foldchange比对
# transcriptome
transcriptome_fc <- read.xlsx("./01_data/Cpm_data/Foldchange.xlsx")

# proteome
proteome_fc <- read.xlsx("../Proteome/01_Data/Foldchange.xlsx")

# common set
common_set <- intersect(transcriptome_fc$Gene, proteome_fc$Gene)

transcriptome_fc_select <- transcriptome_fc[transcriptome_fc$Gene%in%common_set,]
transcriptome_fc_select <- transcriptome_fc_select[!duplicated(transcriptome_fc_select$Gene),]
rownames(transcriptome_fc_select) <- transcriptome_fc_select$Gene
transcriptome_fc_select <- transcriptome_fc_select[,-1]
transcriptome_fc_select <- transcriptome_fc_select[order(rownames(transcriptome_fc_select)),]

proteome_fc_select <- proteome_fc[proteome_fc$Gene%in%common_set,]
proteome_fc_select <- proteome_fc_select[!duplicated(proteome_fc_select$Gene),]
rownames(proteome_fc_select) <- proteome_fc_select$Gene
proteome_fc_select <- proteome_fc_select[, -1]
proteome_fc_select <- proteome_fc_select[order(rownames(proteome_fc_select)), ]

# 合并
merged_df <- cbind(transcriptome_fc_select, proteome_fc_select)
colnames(merged_df) <- c(paste0("Transcript_", c(1:17)), 
                         paste0("Protein_", c(1:17)))

# 5. cor analysis ----
library(corrplot)
correlation_matrix <- cor(transcriptome_fc_select, proteome_fc_select)

# plot
pdf("../FC/fc_cor_cpm.pdf", width = 8,height = 6)
colnames(correlation_matrix) <- c(paste0("Transcript_", c(1:17)))
rownames(correlation_matrix) <- c(paste0("Protein_", c(1:17)))
cor <- corrplot(correlation_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.7, tl.cex = 0.7)
dev.off()

library(ggplot2)

# 6. pca analysis ----
pca_result <- prcomp(t(merged_df), center = TRUE, scale. = TRUE)

# 创建分组信息
sample_names <- rownames(pca_result$x)
group_labels <- rep(paste0("Group_", 1:17), times = 2,)  # 每对Transcript和Protein为一组

# 将分组信息添加到PCA结果数据框中
pca_data <- as.data.frame(pca_result$x)
pca_data$Group <- group_labels  # 添加Group列

# 5. 绘制PCA图并根据分组着色
library(ggplot2)
library(ggrepel)
library(ggsci)
colors <- brewer.pal(12, "Paired")  # "Paired" 提供高对比度颜色
colors <- rep(colors, length.out = 17)  # 扩展到 17 组
pca_data$Group <- factor(pca_data$Group, levels = paste0("Group_", 1:17))  # 确保按顺序排列
# 提取方差贡献率（转换为百分比）
explained_variance <- summary(pca_result)$importance[2, ] * 100

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +  # 绘制PCA散点图
  geom_label_repel(aes(label = rownames(pca_data)),
                   size = 3, alpha = 0.8, label.size = 0.2,
                   max.overlaps = 50,  # 允许更多标签
                   force = 5) +  # 标签排斥力
  labs(title = "PCA Plot of Transcript and Protein Groups",
       x = paste0("PC1 (", round(explained_variance[1], 2), "%)"),
       y = paste0("PC2 (", round(explained_variance[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = colors)  

ggsave(filename = "fc_pca.pdf",plot = pca,device = pdf,
       path = NULL,dpi = 300)

# 7.  boxplot ----
library(tidyr)

fc_data_long <- pivot_longer(merged_df, 
                             cols = everything(), 
                             names_to = "Sample", 
                             values_to = "FC")
fc_data_long$Group <- rep(paste0("Group_", 1:17), times = nrow(fc_data_long)/17)

# 绘制箱线图和散点图结合图
ggplot(fc_data_long, aes(x = Sample, y = FC, fill = Group)) + 
  geom_boxplot(outlier.shape = NA, color = "black") +  # 绘制箱线图，去掉异常值
  geom_point(aes(color = Sample), position = position_jitter(width = 0.2), size = 3) +  # 绘制散点图
  labs(x = "Sample", y = "FC") +
  theme_minimal() +
  scale_fill_manual(values = c("Group_1" = "#1D3557",  # 深海蓝
                               "Group_2" = "#457B9D",  # 雾霾蓝
                               "Group_3" = "#A8DADC",  # 浅海蓝
                               "Group_4" = "#F1FAEE",  # 暗米白
                               "Group_5" = "#E63946",  # 酱紫红
                               "Group_6" = "#F1A7A6",  # 温暖粉
                               "Group_7" = "#2A9D8F",  # 深青色
                               "Group_8" = "#264653",  # 深墨绿色
                               "Group_9" = "#D4A5A5",  # 辣椒红
                               "Group_10" = "#A4A9B3", # 暗灰蓝
                               "Group_11" = "#6A4C93", # 高贵紫
                               "Group_12" = "#1F2A44", # 墨蓝色
                               "Group_13" = "#0F4B5F", # 钢青色
                               "Group_14" = "#EF476F", # 明亮珊瑚红
                               "Group_15" = "#073B4C", # 青黑色
                               "Group_16" = "#118C94", # 天青色
                               "Group_17" = "#F77F00"))+ # 深橙色  
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # 将x轴标签旋转90度

# 8. scatter plot ----
library(ggplot2)
start_column <- 1
end_column <- start_column + 1     # 因为循环3次，所以结束列为 start_column + 2
dir_scatter <- "../FC/fc_scatter_cpm/"


for(i in start_column:end_column)  {
# 提取对应fc
scatter_df <- merged_df[,c(i, i+17)]                      # 每次取 i 和 i+17 作为一对
colnames(scatter_df) <- c("Transcript","Protein")
scatter_df$Transcript <- scale(scatter_df$Transcript)
scatter_df$Protein <- scale(scatter_df$Protein)
scatter_df <- as.data.frame(scatter_df)
scatter_df <- scatter_df[order(scatter_df$Transcript), ] 
cor_value <- cor(scatter_df$Transcript, scatter_df$Protein)

# 绘制散点图
# 创建 PDF 文件名和标题
pdf_filename <- paste0(dir_scatter,"fc_scatter_", i, ".pdf")       # 文件名：fc_scatter_1.pdf...
plot_title <- paste0("MOLM13_2W", i - start_column + 1, 
                     ".vs.WT_", i - start_column + 1)  # 图标题：MOLM13_6W_1.vs.WT_1

ggplot_plot <- ggplot(scatter_df, aes(x = Transcript, y = Protein)) +
  geom_point(shape = 1, size = 3) +  # x 点颜色
  labs(x = "Fc of Transcript", y = "Fc of Protein", title = plot_title) +
  theme_minimal() +
  theme(legend.position = "none")+  # 隐藏图例
  annotate("text", x = min(scatter_df$Transcript), y = max(scatter_df$Protein), 
           label = paste("Correlation: ", round(cor_value, 2)), hjust = 0, vjust = 1, size = 5, color = "black")
ggsave(filename = pdf_filename, plot = ggplot_plot, width = 8, height = 6 )
}

# 9. 热图 ----
library(pheatmap)
pheatmap_df <- merged_df[,c(6:8,23:25)]
pdf("fc_pheatmap.pdf", width = 8,height = 6)
pheatmap_plot <- pheatmap(pheatmap_df, 
         cluster_rows = TRUE,     # 对行进行聚类
         cluster_cols = TRUE,     # 对列进行聚类
         scale = "none",          # 对每行数据进行标准化
         color = colorRampPalette(c("blue", "white", "red"))(50),  # 自定义颜色
         main = "Heatmap of FC",  # 标题
         show_rownames = FALSE,   # 显示基因名
         show_colnames = TRUE     # 显示样本名
)
dev.off()  

# 表达量相关性
proteome_expr <- expr_filter
transcriptome_expr <- expr_filter
common_set <- intersect(transcriptome_expr$gene, proteome_expr$gene)

# 

proteome_expr_select <- proteome_expr[proteome_expr$gene%in%common_set,]
proteome_expr_select <- proteome_expr_select[!duplicated(proteome_expr_select$gene),]
rownames(proteome_expr_select) <- proteome_expr_select$gene
proteome_expr_select <- proteome_expr_select[,-27]
proteome_expr_select <- proteome_expr_select[order(rownames(proteome_expr_select)),]
proteome_expr_select <- proteome_expr_select[,order(colnames(proteome_expr_select))]

transcriptome_expr_select <- transcriptome_expr[transcriptome_expr$gene%in%common_set,]
transcriptome_expr_select <- transcriptome_expr_select[!duplicated(transcriptome_expr_select$gene),]
rownames(transcriptome_expr_select) <- transcriptome_expr_select$gene
transcriptome_expr_select <- transcriptome_expr_select[,-27]
transcriptome_expr_select <- transcriptome_expr_select[order(rownames(transcriptome_expr_select)),]
transcriptome_expr_select <- transcriptome_expr_select[,order(colnames(transcriptome_expr_select))]

proteome_scaled <- scale(proteome_expr_select)
transcriptome_scaled <- scale(transcriptome_expr_select)

# cor analysis 
correlation_matrix_1 <- cor(transcriptome_scaled,proteome_scaled)
colnames(correlation_matrix_1) <- c(paste0("Transcript_", c(1:26)))
rownames(correlation_matrix_1) <- c(paste0("Protein_", c(1:26)))
# plot
pdf("expr_cor.pdf", width = 8,height = 6)
cor <- corrplot(correlation_matrix_1, method = "color", type = "upper", 
                tl.col = "black", tl.srt = 45, addCoef.col = "black",
                number.cex = 0.7, tl.cex = 0.7)
dev.off()

# 二、mean sample logfc
# 1. data_input ----
transcriptome_expr_select <- read.csv()
proteome_expr <- read.csv("../FC/Proteome_selected_expr.csv")
colnames(proteome_expr)[1] <- "gene"

proteome_expr_select <- proteome_expr[proteome_expr$gene%in%common_set,]
proteome_expr_select <- proteome_expr_select[!duplicated(proteome_expr_select$gene),]
rownames(proteome_expr_select) <- proteome_expr_select$gene
proteome_expr_select <- proteome_expr_select[,-1]
proteome_expr_select <- proteome_expr_select[order(rownames(proteome_expr_select)),]
proteome_expr_select <- proteome_expr_select[,order(colnames(proteome_expr_select))]

transcriptome_expr_select <- transcriptome_expr[transcriptome_expr$gene%in%common_set,]
transcriptome_expr_select <- transcriptome_expr_select[!duplicated(transcriptome_expr_select$gene),]
rownames(transcriptome_expr_select) <- transcriptome_expr_select$gene
transcriptome_expr_select <- transcriptome_expr_select[,-1]
transcriptome_expr_select <- transcriptome_expr_select[order(rownames(transcriptome_expr_select)),]
transcriptome_expr_select <- transcriptome_expr_select[,order(colnames(transcriptome_expr_select))]

# 函数：计算logFC
calculate_logfc <- function(expr_data, group1_pattern, group2_pattern) {
  
  # log2
  expr_data <- log2(expr_data)
  
  # 获取group1（如2W）和group2（如WT）的列名
  group1_samples <- grep(group1_pattern, colnames(expr_data), value = TRUE)
  group2_samples <- grep(group2_pattern, colnames(expr_data), value = TRUE)
  
  # 计算组内每个样本的均值
  group1_mean <- rowMeans(expr_data[, group1_samples], na.rm = TRUE)
  group2_mean <- rowMeans(expr_data[, group2_samples], na.rm = TRUE)
  
  # 计算log2 Fold Change (logFC)
  logfc <- group1_mean - group2_mean
  
  return(logfc)
}

cell_lines <- c("MOLM13", "MV4_11", "OCI_AML2")
time_points <- c("2W", "6W")

# 结果保存的列表
logfc_results <- list()

# 迭代处理每个细胞系和时间点

for (cell in cell_lines) {
  for (time_point in time_points) {
    # 构造2W和WT样本名称的模式
    group1_pattern <- paste0(cell, "_", time_point)
    group2_pattern <- paste0(cell, "_WT")
    
    # 提取和计算logFC
    logfc_results[[paste0(cell, "_", time_point,".vs.WT")]] <- calculate_logfc(transcriptome_expr_select, group1_pattern, group2_pattern)
  }
}
transcriptome_logfc_results <- as.data.frame(logfc_results)
proteome_logfc_results <- as.data.frame(logfc_results)

# 合并
merged_expr_df <- cbind(transcriptome_logfc_results, proteome_logfc_results)
colnames(merged_expr_df) <- c(paste0("Trans_", colnames(transcriptome_logfc_results)), 
                         paste0("Prote_", colnames(proteome_logfc_results)))

# 2. cor analysis ----
library(corrplot)
expr_correlation_matrix <- cor(transcriptome_logfc_results, proteome_logfc_results)

# plot
pdf("../FC/fc_all_cpm_cor.pdf", width = 8,height = 6)
colnames(expr_correlation_matrix) <- c(paste0("Trans_", colnames(transcriptome_logfc_results)))
rownames(expr_correlation_matrix) <- c(paste0("Prote_", colnames(proteome_logfc_results)))
cor <- corrplot(expr_correlation_matrix, method = "color", type = "upper", 
                tl.col = "black", tl.srt = 45, addCoef.col = "black",
                number.cex = 0.7, tl.cex = 0.7)
dev.off()
