#' Analysis of Avariance, ANOVA 
#'

# Packages ----
library(readr)
library(openxlsx)
library(ggplot2)

# Data input ----
data <- read.csv("./01_data/gene_tpm_matrix.csv", row.names = 1)
group <- read.xlsx("./01_data/group_info.xlsx")
data <- data[,-grep("OCI_4W",colnames(data))]
group <- group[-grep("OCI_4W",group$id),]
head(data)
head(group)
rownames(group) <- group$id
data <- data[,-grep('MOLM13|MV4',colnames(data))]
group <- group[-grep('MOLM13|MV4',group$id),]

# 过滤低表达基因
source("../../Code_Box/Transcriptome/run_filter_coding_gene.R")
data_clean <- filter_coding_genes(data)    # 去除非编码蛋白基因
data_clean <- data_clean[,!colnames(data_clean)%in%c("ense","symbol")]
min_sample <- ceiling(length(colnames(data_clean)) / 2)  # 向上取整
exprSet_filtered <- data_clean[rowSums(data_clean > 1) >= min_sample,]  
dim(exprSet_filtered)

data_target <- exprSet_filtered
group_target <- group

# 创建分组因子
group_target$group <- factor(group_target$group, levels =  c("WT", "Low", "High"))  

# 转置表达矩阵：行为样本，列为蛋白
expr_matrix <- t(log2(data_target+1))        # log2 标准化处理

# 对比样本名和分组，需确保完全对应
identical(rownames(expr_matrix),rownames(group_target)) 

if(T){
# ----------------------------- ANOVA ------------------------------------------
anova_results <- apply(expr_matrix, 2, function(x) {
  fit <- aov(x ~ group_target$group)
  summary(fit)[[1]][["Pr(>F)"]][1]  # 提取 p 值
})
# ------------------------- Welch ANOVA -------------------------------------
welch_results <- apply(expr_matrix, 2, function(x) {
  oneway.test(x ~ group_target$group, var.equal = FALSE)$p.value
})
# ------------------------- Kruskal-Wallis -------------------------------------
kruskal_results <- apply(expr_matrix, 2, function(x) {
  kruskal.test(x ~ group_target$group)$p.value
})
# 汇总结果
stats_table <- data.frame(
  Protein_ID = colnames(expr_matrix),
  ANOVA_p = anova_results,
  Welch_p = welch_results,
  KW_p = kruskal_results
)

# 添加 FDR 校正
stats_table$ANOVA_FDR <- p.adjust(stats_table$ANOVA_p, method = "BH")
stats_table$Welch_FDR <- p.adjust(stats_table$Welch_p, method = "BH")
stats_table$KW_FDR <- p.adjust(stats_table$KW_p, method = "BH")

# 筛选显著蛋白（如 ANOVA 或 Kruskal-Wallis）
sig_anova <- subset(stats_table, ANOVA_p < 0.05)
sig_welch <- subset(stats_table, Welch_p < 0.05)
sig_kw    <- subset(stats_table, KW_p < 0.05)
}

# Res output ----
write.csv(stats_table , "./03_Result/05_ANOVA/Result_table.csv")
