
  title: "RNA-Seq_normalization.R"

# data input ----

## 加载packages和function ----
setwd("/download/RStudio/AML_project/Transcriptome/")
source("./02_code/QC_PCA.R")
source("./02_code/QC_boxplot.R")
source("./02_code/QC_heatmap.R")
source("./02_code/run_limma_DE.R")
source("./02_code/run_enrichment_analysis.R")
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(pathview)
library(ggplot2)
## group input ----

  # 导入分组信息
data_group <- read_xlsx("./01_data/Counts_data/aml_group_merge.xlsx")
table(data_group$group)
data_group <- as.data.frame(data_group)
head(data_group,n=6) #检验

# 添加bacth信息
data_group$batch <- c(rep("Batch1",23),rep("Batch2",4))

# 添加实验对照分组信息
data_group$condition <- c(rep("Experiment",5),"Control",
                          rep("Experiment",5),rep("Control",3),
                          rep("Experiment",6),rep("Control",5),
                          rep("Experiment",2))
# 配色设置
value_colour <- c("MOLM13" = "#E64B35FF",# Experimental group
                  "MV4_11" = "#4DBBD5FA",# other group1
                  "OCI_AML2" = "#F2A200"# other group2
                  )

value_colour_batch <- c("Batch1" = "#E64B35FF",# Experimental group
                        "Batch2" = "#4DBBD5FA"# other group1
                  )
value_colour_condition <- c("Experiment" = "#E64B35FF",# Experimental group
                    "Control" = "#4DBBD5FA"# other group1
)
rownames(data_group) <- data_group$id
## RNA-Seq matrix input ---- 

data_input <- read_csv("./01_data/Tpm_data/gene_tpm_matrix_merged.csv")
head(data_input,n=6) #检验是否读取正确
data_input <- as.data.frame(data_input)
names(data_input)[1] <- "gene_id"
rownames(data_input) <- data_input$gene_id
data_input <- data_input[,-1]
colnames(data_input) <- gsub("-","_",colnames(data_input))
# 处理RNA-Seq数据，将低表达的基因删除
data_before <- data_input
#data_before<-data_input[apply(data_input,1,sum)>15,] 

data_before<-data_input[apply(data_input,1,function(x)sum(x>1)>0.5*ncol(data_input)),] 

# 设置输出目录----
#dir.create("03_result")
dir <- "./03_result/QC/"
data_group_1 <- data_group[,-c(3,4)]
data_batch <- data_group[,-c(2,4)]
colnames(data_batch)[2] <- "group"
data_condition <- data_group[,-c(2,3)]
colnames(data_condition)[2] <- "group"

# boxplot -----------------------------------------------------------------
pdf(file = paste0(dir,"QC_boxplot.pdf"),
    width = 6,
    height = 4)
QC_boxplot(data_before,data_group = data_group,
           value_colour = value_colour,
           title = " original_data")
dev.off()

# heatmap -----------------------------------------------------------------
# 进行heatmap图绘制时，需要去除标准差为零的基因，否则会出现缺失值
# 去除标准差不等于0后的结果赋值为data_before_1
data_before_1<-data_before[apply(data_before,1,sd)!=0,]

# group分组

pdf(file = paste0(dir,"QC_heatmap_group.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_before_1,data_group = data_group_1,
           value_colour = value_colour)
dev.off()

# batch分组

pdf(file = paste0(dir,"QC_heatmap_batch.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_before_1,data_group = data_batch,
           value_colour = value_colour_batch)
dev.off()

# condition分组

pdf(file = paste0(dir,"QC_heatmap_condition.pdf"),
    width = 6,
    height = 6)
QC_heatmap(data = data_before_1,data_group = data_condition,
           value_colour = value_colour_condition)
dev.off()
# PCA ---------------------------------------------------------------------

# group
pdf(file = paste0(dir,"QC_pca_adjusted_1.pdf"),
    width = 15,
    height = 15)
QC_PCA(data = data_before,
       data_group = data_group_1,
       value_colour = value_colour)
dev.off()

# batch
pdf(file = paste0(dir,"QC_pca_batch.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = data_before,
       data_group = data_batch,
       value_colour = value_colour_batch)
dev.off()

# condition
pdf(file = paste0(dir,"QC_pca_condition.pdf"),
    width = 10,
    height = 10)
QC_PCA(data = data_before,
       data_group = data_condition,
       value_colour = value_colour_condition)
dev.off()
# DEseq2 -----------------------------------------------------------

## 加载包 ----
library(tidyverse)
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library(DESeq2)
library(openxlsx)

## count数据导入----
dat <- read.xlsx("./01_data/aml数据_time/OCI_AML2_counts_sva_WT_6W.xlsx")
# 检验数据格式
head(dat,n=8) 

# 提取geneid列作为行名
colnames(dat)[1] <- "gene_id"
rownames(dat) <- dat$gene_id
dat <- dat[,-1]

## group信息导入----
colData <- read.xlsx('./01_data/aml数据_time/OCI_AML2_group_WT_6W.xlsx')

# 检验数据格式
head(colData,n=8) 

# 显示分组信息
table(colData$group) 
colData$group <- factor(colData$group, levels = c("OCI_AML2_WT","OCI_AML2_D"))

# R会默认将因子向量的第一个水平作为参考（对照）组，然后将其他水平与这个参考组进行比较。OCI_AML2/MV4_11/MOLM13
table(colData$group)

## 添加批次信息 ----
#batch<-factor(c(rep(1,6),2,rep(1,3)))
#colData$batch<-batch

# 让表达矩阵的样本与分析信息表保持一致
dat <- dat[, colData$id] 

## 创建DESeqDataSet对象 ----
dds <- DESeqDataSetFromMatrix(countData = dat, 
                              colData = colData, 
                              design = ~group)

# 过滤一些低表达的genes
dds <- dds[rowSums(counts(dds)) > 15, ]
#dds <- dds[ rowSums(counts(dds[,1:5])) > 15, ]
#dds <- dds[ rowSums(counts(dds[,6:8])) > 15, ]

# 估计每个样本的大小因子，用于归一化原始的表达量数据
dds <- estimateSizeFactors(dds) 

## 差异表达分析 ----
dds <- DESeq(dds)

## 导出标准化结果 ----
meta <- rowSums(dat)
View(as.data.frame(table(meta)))
ggplot(as.data.frame(table(meta)),aes(x = log2(Freq))) + 
  geom_density()
exprSet <- dat
rld <- rlogTransformation(dds)## 得到经过DESeq2软件normlization的表达矩阵！
vsd <- vst(dds)
exprSet_new=assay(rld)

# 绘图 ----

## 绘制par图 ----
png("MOLM13.png",width = 800,height = 600)
par(cex.axis = 0.8, cex.lab = 0.7,mfrow=c(2,2),mar=c(6,4,4,2))
n.sample=ncol(exprSet_new)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(as.matrix(exprSet))
hist(exprSet_new)
dev.off()

## 绘制PCA图 ----
raw <- SummarizedExperiment(counts(dds, normalized=FALSE),
                            colData=colData(dds))
nor <- SummarizedExperiment(counts(dds, normalized=TRUE),
                            colData=colData(dds))
vsd <- vst(dds)
rld <- rlog(dds)
pdf("MOLM13_PCA.pdf")
plotPCA( DESeqTransform(raw), intgroup= "group" )
plotPCA( DESeqTransform(nor), intgroup= "group" )
plotPCA(vsd, intgroup= "group" )
plotPCA(rld, intgroup= "group" )
dev.off()

# 提取结果，results函数中contrast参数来指定比较时
# contrast <- c("condition", "treated", "control")
# R会默认将因子向量的第一个水平作为参考（对照）组，
# 然后将其他水平与这个参考组进行比较。
# 所以要将因子水平对调一下 rev(levels(group)))

## result差异表达分析 ----
res <- results(dds, contrast = c("group",rev(levels(colData$group))))
DE <- as.data.frame(res) 

## 删除NA值 ----
DE<-na.omit(DE)

## 保存删除NA值后的DE ----
write.csv(DE,file = "./03_result/DE/OCI_AML2_WT_vs_OCI_AML2_D/DE(NA.omit)_WT_6W.csv")

# volcano plot ------------------------------------------------------------
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
# group 1为实验组
# group 2为对照组
#OCI_AML2/MV4_11/MOLM13/Venetoclax
group_1 <- "OCI_AML2_WT"
group_2 <- "OCI_AML2_VEN"
res_data <- DE
data <- res_data
#data <- res_data[res_data[,"padj"] <= 1,]
# 颜色划分padj <0.05，且log2FC绝对值大于sd(tempOutput$logFC)*3
# 设定cutoff ,cutoff要以DE原始结果为准
# 计算阈值
cutoff <- 1 # sd(res_data$log2FoldChange) * 3

data$sig[data$pvalue >= 0.05 | abs(data$log2FoldChange < cutoff)] <- "Not"

data$sig[data$pvalue < 0.05 & data$log2FoldChange >= cutoff] <- "Up"

data$sig[data$pvalue < 0.05 & data$log2FoldChange <= -cutoff] <- "Down"

# table统计上下调基因数量
table(data$sig)

# 得到DEseq2结果
input <- data

## DEseq2结果导出 ----
# 从初始数据中提取差异基因
output<-dat[rownames(dat)%in%rownames(input),]

# 合并DEseq2数据
DEseq2_output<-cbind(output,input)

# 将结果写入相应文件名导出
write.csv(DEseq2_output,file = "./03_result/DE/OCI_AML2_WT_vs_OCI_AML2_D/DE(sig)_WT_6W.csv")
# MOLM13/MV4_11/OCI_AML2DE(sig)_2W_6W.csv

library(ggrepel)
library(ggplot2)
# 提取基因名
#y <- rownames(data)
#gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
#id<-data.frame(gene)
#data<-cbind(id,data)
# 筛选上调中显著性top10的基因
#top10_up <- filter(data, sig == 'Up') %>%  # 筛选出上调的基因
# distinct(gene, .keep_all = TRUE) %>%   # 去除重复的基因，保留第一个出现的行
#top_n(10, -log10(padj))               # 选择-P.Value值最大的前10个基因

# 筛选下调中显著性top10的基因
#top10_down <- filter(data, sig == 'Down') %>%  # 筛选出上调的基因
#distinct(gene, .keep_all = TRUE) %>%   # 去除重复的基因，保留第一个出现的行
# top_n(10, -log10(padj)) 
# 绘图
volc <- ggplot(data = data, aes(x = log2FoldChange,
                                y = -log10(pvalue),
                                color = sig)) +
  geom_point(alpha = 0.9) +  theme_classic() +
  theme(panel.grid = element_blank(),strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#4DBBD5","grey80","#E64B35")) +
  geom_vline(xintercept = c(-cutoff,cutoff),lty = 4,lwd = 0.6,alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05),lty = 4,lwd = 0.6,alpha = 0.8) +
  labs(title = paste0(group_1,"-",group_2))+
  xlim(c(-30, 25))#+
geom_label_repel(data = top10_up, 
                 aes(x = log2FoldChange, y = -log10(padj), label = gene),
                 show.legend = FALSE,size = 2,max.overlaps = Inf)+  # 添加上调基因的标签
  geom_label_repel(data = top10_down, 
                   aes(x = log2FoldChange, y = -log10(padj), label = gene),
                   show.legend = FALSE,size = 2,max.overlaps = Inf)  # 添加下调基因的标签
volc

# KEGG GO -----------------------------------------------------------------

# 设置输出目录----OCI_AML2/MV4_11/MOLM13
#dir.create("../GO+KEGG_enrichment/OCI_AML2_WT_vs_OCI_AML2_D_WT_6W")
dir <- "./03_result/GO&KEGG_DEseq2/MV4_11_WT_vs_MV4_11_6W"

# 差异分析结果输入 ----
res_data<-read.csv("./03_result/DE_DEseq2/MV4_11/DE(sig)_WT_6W.csv")

colnames(res_data) # 查看列名，根据需要进行下一项更改
colnames(res_data)[which(colnames(res_data) == "log2FoldChange")] <- "logFC"
colnames(res_data)
GeneSymbol <- subset(res_data,pvalue< 0.05)

## 提取基因名（symbol） ----
y <- GeneSymbol$X
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
id<-data.frame(gene)
GeneSymbol<-cbind(id,GeneSymbol)

# 检查有无重复项 
duplicates <- duplicated(GeneSymbol)
a <- print(GeneSymbol[duplicates,])

# 设置cutoff值
cutoff <- 1

## 提取down-regulated genes ----
down_genes <- subset(GeneSymbol, logFC < -cutoff)

## GO-down ----
### 设置数据库 ----
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库

### gene ID转换 ----
gene <- bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### GO富集分析----
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)# 设定q值阈值

## KEGG-down ----
### 设置数据库 ----
KEGG_database <- 'hsa'  # KEGG是hsa数据库

### gene ID转换 ----
gene <- bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### KEGG富集分析----
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
## GO、KEGG结果整合 ----
result <- list(enrichGO = kkd, enrichKEGG = kk)

## 结果标记为下调 ----
result_down <- result
kkd_down <- result_down$enrichGO
kk_down <- result_down$enrichKEGG

## 导出下调enrichGO数据为csv文件 ----
write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)

## 数据可视化并导出为PDF ----
pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
print(p1)
dev.off()

##导出下调enrichKEGG数据为csv文件 ----
write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)

## 数据可视化并导出为PDF ----
pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
p2 <- dotplot(kk_down,showCategory = 10)
print(p2)
dev.off()

## 提取up-regulated genes ----
up_genes <- subset(GeneSymbol, logFC > cutoff)

## GO-up ----
### 设置数据库 ----
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库

### gene ID转换 ----
gene <- bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### GO富集分析----
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, # 设定q值阈值
                readable = T)
## KEGG-up ----

### 设置数据库 ----
KEGG_database <- 'hsa'  # KEGG是hsa数据库

### gene ID转换 ----
gene <- bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### KEGG富集分析----
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
## GO、KEGG结果整合 ----
result <- list(enrichGO = kkd, enrichKEGG = kk)

# 结果标记为上调
result_up <- result
kkd_up <- result_up$enrichGO
kk_up <- result_up$enrichKEGG

## 导出上调enrichGO数据为csv文件 ----
write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)

## 数据可视化并导出为PDF ----
pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off() 

## 导出上调enrichKEGG数据为csv文件 ----
write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)

## 数据可视化并导出为PDF ----
pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
p4 <- dotplot(kk_up,showCategory = 10 )
print(p4)
dev.off()

# 统计上下调的通路的数量 ----
## 下调GO ----
table(kkd_down@result$p.adjust<0.05)
## 下调KEGG ----
table(kk_down@result$p.adjust<0.05)
## 上调GO ----
table(kkd_up@result$p.adjust<0.05)
## 上调KEGG ----
table(kk_up@result$p.adjust<0.05)
