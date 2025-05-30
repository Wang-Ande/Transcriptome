# 1. Library packages ----
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
library(openxlsx)
library(enrichplot)
library(patchwork)
library(gen)
source("02_code/run_limma_DE.R")                 # load前检查是否需要更改！
source("02_code/run_enrichment_analysis.R")

# 2. DE_limma ----

## 2.1 Expr input ----
exprSet_raw <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
exprSet_raw <- as.data.frame(exprSet_raw)

# 整理格式
colnames(exprSet_raw)
colnames(exprSet_raw)[1] <- "gene_id"
rownames(exprSet_raw) <- exprSet_raw$gene_id
exprSet_raw <- exprSet_raw[,-1]
colnames(exprSet_raw) <- gsub("cas9","", colnames(exprSet_raw))

## 2.2 log2(expr+1) ----

## 2.3 Filter low genes ----
# approach 1
exprSet_filtered <- exprSet_raw[apply(exprSet_raw,1,
                                   function(x)sum(x>1)>0.5*ncol(exprSet_raw)),] 

# approach 2
exprSet_raw_log2 <- log2(exprSet_raw+1)
exprSet_filtered <- exprSet_raw_log2[apply(exprSet_raw_log2,1,mean)>1,] # 13454
exprSet_filtered <- 2^exprSet_filtered-1

# 提取编码蛋白基因
#gtf_file <- "01_data/Homo_sapiens.GRCh38.113.gtf/Homo_sapiens.GRCh38.113.gtf"
#txdb <- makeTxDbFromGFF("01_data/Homo_sapiens.GRCh38.113.gtf/Homo_sapiens.GRCh38.113.gtf", format = "gtf")

#txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
#gene_info <- genes(txdb)

#gtf_data <- import(gtf_file)
#head(mcols(gtf_data))
#coding_genes <- gtf_data[gtf_data$gene_biotype == "protein_coding", ]
#coding_gene_ids <- unique(coding_genes$gene_id)

## 2.4 Group input ----
colData <- read_xlsx("./01_data/Counts_data/aml_group_merge.xlsx")
colData <- as.data.frame(colData)
colnames(colData)

# 整理样本
colData$id <- gsub("cas9","", colData$id)         # 去除样本名多余信息
colData <- colData[-grep("4w",colData$id),]       # 去除4w样本
rownames(colData) <- colData$id
table(colData$group) 

targeted_group <- colData

# VEN-WT 组别
if(T){
targeted_group$group <- paste0(targeted_group$group, "_VEN")
targeted_group[grep("WT",targeted_group$id),2] <- gsub("_VEN","_WT",
                                                       targeted_group[grep("WT",targeted_group$id),2])
table(targeted_group$group)
}

# 2W/6W-WT 组别
if(T){
targeted_group[grep("2w",targeted_group$id),2] <- gsub("_VEN","_2W",
                                                       targeted_group[grep("2w",targeted_group$id),2])
targeted_group[grep("6w",targeted_group$id),2] <- gsub("_VEN","_6W",
                                                       targeted_group[grep("6w",targeted_group$id),2])
table(targeted_group$group)
}

# 显示分组信息
#colData$group <- gsub("OCI_AML2_D", "OCI_AML2_2W", colData$group)
# MOLM13/MV4_11/OCI_AML2
## 2.5 Set group ----
group_1 <- "MOLM13_VEN"    # group 1为实验组
group_2 <- "MOLM13_WT"    # group 2为对照组

## 2.6 res output ---- 
source("02_code/run_limma_DE.R")                 # load前检查是否需要更改！
result_merge <- run_limma_DE(exprSet = exprSet_filtered,
                             colData = targeted_group,
                             group_1 = group_1,
                             group_2 = group_2,
                             log2 = TRUE,
                             logfc_threshold = 0.585,        # log2fc值
                             pvalue_threshold = 0.05, 
                             qvalue_threshold = NULL,
                             dir = "03_result/DE_limma")
table(result_merge$DE_res$Sig)

# 3. GO KEGG --------------------------------------------------------------
## 3.1 设置输出目录----
#OCI_AML2/MV4_11/MOLM13
# 指定文件夹路径
folder_path <- "./03_result/GO&KEGG_limma/MV4_11_WT_vs_MV4_11_VEN/"
dir <- "./03_result/GO&KEGG_limma/MOLM13_WT_vs_MOLM13_2W/"

# 检查文件夹是否存在
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
  print(paste("Folder", folder_path, "created."))
} else {
  print(paste("Folder", folder_path, "already exists."))
}

## 3.2 DE_res input ----
DE_result <- read.csv('./03_result/DE_limma/MOLM13_2W_vs_MOLM13_WT/DE.csv')

## 3.3 设置P.Value ----
Genesymbol <- subset(DE_result, P.Value < 0.05)

## 3.4 设置cutoff值 ----
cutoff <- 0.585

# 转换基因名 
y <- Genesymbol$Row.names
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
id <- data.frame(gene)
Genesymbol <- cbind(id,Genesymbol)

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

if(T){
## 3.5 down genes ----
down_genes <- subset(Genesymbol, logFC < -cutoff)

# gene ID转换 
gene <- bitr(down_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

### 3.5.1 GO ----
# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)# 设定q值阈值

### 3.5.2 KEGG ----
# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

## GO、KEGG结果整合 
result <- list(enrichGO = kkd, enrichKEGG = kk)

# 结果标记为下调 
result_down <- result
kkd_down <- result_down$enrichGO
kk_down <- result_down$enrichKEGG

### 3.5.3 下调res_output ----
# 导出下调enrichGO 
write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
print(p1)
dev.off()

# 导出下调enrichKEGG
write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
p2 <- dotplot(kk_down,showCategory = 10)
print(p2)
dev.off()
}

if(T){
## 3.6 up genes ----
up_genes <- subset(Genesymbol, logFC > cutoff)

### 3.6.1 GO-up ----
# gene ID转换 
gene <- bitr(up_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, # 设定q值阈值
                readable = T)
### 3.6.2 KEGG-up ----

# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

# GO、KEGG结果整合
result <- list(enrichGO = kkd, enrichKEGG = kk)

# 结果标记为上调
result_up <- result
kkd_up <- result_up$enrichGO
kk_up <- result_up$enrichKEGG

### 3.6.3 上调res_output ----

# 导出上调enrichGO
write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
print(p3)
dev.off()

# 导出上调enrichKEGG
write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)

# dotplot
pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
p4 <- dotplot(kk_up,showCategory = 10 )
print(p4)
dev.off()
}


## 3.7 统计上下调的通路的数量 ----
# 下调GO 
table(kkd_down@result$p.adjust<0.05)
# 下调KEGG 
table(kk_down@result$p.adjust<0.05)
# 上调GO 
table(kkd_up@result$p.adjust<0.05)
# 上调KEGG 
table(kk_up@result$p.adjust<0.05)
