# ClusterGVis

# 1. Packages ----
# Note: please update your ComplexHeatmap to the latest version!
library(pak)
pkg_install("junjunlab/ClusterGVis")
library(ClusterGVis)
library(grid)
library(ComplexHeatmap)
library(clusterProfiler)
library(TCseq)
library(e1071)
library(tkWidgets)
library(Mfuzz)
library(circlize)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(openxlsx)

# 2. Data input ----
# expr input(standardised)
expr <- read.csv("./01_data/gene_tpm_matrix.csv",row.names = 1)
colnames(expr)
expr <- expr[,-grep("OCI_4W",colnames(expr))]
expr <- expr[,-grep('MOLM13|MV4',colnames(expr))]

# expr <- expr[expr$Sig != "stable",c(1:11)]
# 用anova结果
DE_gene <- read.csv("./03_result/05_ANOVA/Result_table.csv", row.names = 1)
sig_welch <- subset(DE_gene, Welch_FDR < 0.05)
# top5000 <- sig_welch[order(sig_welch$Welch_FDR), ][1:5000, ]
expr <- expr[rownames(expr)%in%rownames(sig_welch),]

# 提取gene名
y <- rownames(expr)
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
expr$gene <- gene 
rownames(expr) <- expr$gene
expr <- expr[,!colnames(expr)%in%c("gene")]
# 若提示有重复，进行去重处理
# 1. 处理 NA 值
expr <- expr[!is.na(expr$gene), ]  # remove NA rows

# 2. 处理重复基因
expr <- aggregate(. ~ gene, data = expr, FUN = max)

# 3. 设定 rownames
rownames(expr) <- expr$gene
expr <- expr[,-1]

# group input
sample_group <- read.xlsx("./01_data/group_info.xlsx")
rownames(sample_group) <- sample_group$id
sample_group <- sample_group[-grep("OCI_4W",sample_group$id),]
sample_group <- sample_group[-grep('MOLM13|MV4',sample_group$id),]
table(sample_group$group)

# check
identical(rownames(sample_group),colnames(expr))

## 2.1 Data preprocess ----
# 要求输入归一化之后的数据（已进行TPM）
expr_log2 <- log2(expr+1)
colnames(expr_log2)

# 按细胞系分析（可选，根据目的是否需要这一步）
# source("./02_Code/extract_expr_by_cell.R")
# process_cell_line("OCI")  # 生成 expr_cellline

## 2.2 caculate average ----
### 2.2.1 approach 1 -----------------------------------------------------------
library(limma)
# limma::avereps：这是来自 limma 包的函数，avereps 用于对重复数据进行平均值计算。
# avereps 会根据指定的 ID 进行分组，并对相同 ID 的数据取平均值
colnames(expr_log2) <- sample_group$group    # 修改表达矩阵的列名
avereps_df <- t(limma::avereps(t(expr_log2) , ID = colnames(expr_log2))) #对用一组的表达值取平均
head(avereps_df)
avereps_df <- avereps_df[,order(colnames(avereps_df),decreasing = T)]    # 设置聚类方向
save(avereps_df,file = './03_result/07_C-Means/Control-Low-High/avereps_df.Rdata')

### 2.2.2 approach 2 -----------------------------------------------------------
# 处理列名，去掉 `.1`, `.2` 这种后缀
colnames(expr_molm13) <- sub("\\..*", "", colnames(expr_molm13))
selected_mean <- apply(expr_molm13[, c(2, 7, 8)], 1, mean, na.rm = TRUE)
selected_mean <- as.data.frame(selected_mean)
expr_molm13[,"34"] <- selected_mean$selected_mean
expr_molm13 <- expr_molm13[,-c(7,8)]

# 按列名数值大小排序
expr_molm13 <- expr_molm13[, order(as.numeric(colnames(expr_molm13)))]
avereps_df <- expr_molm13

# 3. Cluster Analysis ----
## 3.1 Set output path ----
dir_cl <- "./03_result/07_C-Means/Control-Low-High/"
exps <- as.matrix(avereps_df) 

# check optimal cluster numbers
pdf(paste0(dir_cl,"Elbow_plot.pdf"),width = 6,height = 5)
getClusters(obj = exps)
dev.off()
# choose cluster.num =  5

## 3.2 mfuzz for clustering ----
cm <- clusterData(obj = exps,
                  # scaleData = TRUE,
                  seed = 123,
                  cluster.method = "mfuzz",
                  cluster.num = 4)
save(cm,file = paste0(dir_cl,"C-Means_res.Rdata"))
# 输出excel格式，方便查看
write.xlsx(cm$wide.res, file = paste0(dir_cl,"C-Means_res.xlsx"))


## 3.3 kmeans for clustering ----
km <- clusterData(obj = as.data.frame(exps),
                  # scaleData = TRUE,
                  cluster.method = "kmeans",
                  cluster.num = 4,
                  seed = 123)
save(km,file = paste0(dir_cl,"K-Means_res.Rdata"))
# 输出excel格式，方便查看
write.xlsx(km$wide.res, file = paste0(dir_cl,"K-Means_res.xlsx"))

# 4. Plot ----
## 4.1 plot line only ----
visCluster(object = cm,
           plot.type = "line")

# change color
visCluster(object = cm,
           plot.type = "line",
           ms.col = c("green","orange","red"))

# remove meadian line
cairo_pdf(paste0(dir_cl,"Cluster_lineplot.pdf"),width = 8.5,height = 5)
visCluster(object = cm,
           plot.type = "line",
           add.mline = FALSE)
dev.off()

## 4.2 plot heatmap only ----
pdf(paste0(dir_cl,'Cluster_heatmap.pdf'),height = 10,width = 6)
visCluster(object = km,
           plot.type = "both",
           column_names_rot = 45)
dev.off()

# 5. Enrich analysis ----
library(org.Hs.eg.db)

##  5.1 enrich for clusters ----
enrich <- enrichCluster(object = cm,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        id.trans = TRUE,
                        fromType = "SYMBOL",
                        toType = c("ENTREZID"),
                        readable = TRUE,
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 123)
save(enrich,file = (paste0(dir_cl,"enrich.Rdata")))

# check
head(enrich,3)

# cluster num
cl_num <- 4   # ggsci::pal_d3()(cl_num) 

## 5.2 plot ----
cairo_pdf(paste0(dir_cl,'Cluster_enrichplot.pdf'),
          height = 6.5,width = 11,onefile = F)
visCluster(object = cm,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes.side = "left",
           genes.gp = c('italic',fontsize = 12,col = "black"),
           annoTerm.data = enrich,
           line.side = "left",
           go.col = rep(ggsci::pal_d3()(cl_num),each = 5), 
           go.size = "pval")
dev.off()

## 5.3 all Enrichpathway ----
# The above fun() only shows the top 5 pathways, so we need the next fun() to display all pathways.
# Gene prepare
# source("./02_Code/Extract_genes.R")   # c-means专用
load("./03_result/06_K-Means/Control-Low-High/K-Means_res.Rdata")

# subset target cluster gene
targeted_genes <- km$wide.res[km$wide.res$cluster == 3,]

# set database
GO_database <- 'org.Hs.eg.db' 
KEGG_database <- 'hsa'         

# gene ID转换 
gene <- clusterProfiler::bitr(targeted_genes$gene, fromType = 'SYMBOL', 
                              toType = 'ENTREZID', OrgDb = GO_database)

# GO 
# GO富集分析
go <- clusterProfiler::enrichGO(gene = gene$ENTREZID, 
                                 OrgDb = GO_database, 
                                 keyType = "ENTREZID", 
                                 ont = "ALL",          # (ALL,BP,CC,MF）
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 1,
                                 readable = T)   

# KEGG 
# KEGG富集分析
kegg <- clusterProfiler::enrichKEGG(gene = gene$ENTREZID,
                                  keyType = "kegg",
                                  organism = KEGG_database,
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 1,
                                  use_internal_data = T)
# 设置geneID可读格式
kegg <- setReadable(kegg, OrgDb='org.Hs.eg.db', keyType='ENTREZID')

#GO、KEGG结果整合 
result <- list(enrichGO = go, enrichKEGG = kegg)
GO_res <- result$enrichGO
KEGG_res <- result$enrichKEGG

# Res output
dir_enrich <- "./03_result/06_K-Means/Control-Low-High/Cluster_3/"
if (!dir.exists(dir_enrich)) {
  dir.create(dir_enrich)
  print(paste("Folder", dir_enrich, "created."))
} else {
  print(paste("Folder", dir_enrich, "already exists."))
}
if(T){
# output enrichGO 
write.xlsx(GO_res@result, file = paste0(dir_enrich, "GO.xlsx"))
pdf(file = paste0(dir_enrich, "GO.pdf"), width = 6, height = 8)
p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') +
  theme(axis.text.y = element_text(angle = 0, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))  # Show a maximum of 40 characters 
print(p1)
dev.off()

# output enrichKEGG
write.xlsx(KEGG_res@result, file = paste0(dir_enrich, "KEGG.xlsx"))
pdf(file = paste0(dir_enrich, "KEGG.pdf"), width = 6, height = 5)
p2 <- dotplot(KEGG_res,showCategory = 10)
print(p2)
dev.off()
}

# dotplot
# ggplot2画气泡图，scale_color_gradient设置蓝红配色
library(ggplot2)
library(stringr)

# Enrich_res input 
Enrich_res <- read.xlsx("./03_Result/C_means_cluster/")
Enrich_res <- Enrich_res[Enrich_res$pvalue<0.05,]
# 按 p 值升序排序后，取前 20 行（最显著的 20 个结果）
Enrich_res <- head(Enrich_res, 20)                    # 取前 20 行
Enrich_res <- Enrich_res[order(-Enrich_res$pvalue), ]  # 按 p 值排序

Enrich_res$Description <- factor(Enrich_res$Description, 
                                 levels = unique(Enrich_res$Description))
#Enrich_res <- Enrich_res[order(Enrich_res$Count),]

pdf(file =paste0(dir_enrich,"KEGG_dw_ms_0_selected.pdf") ,
    width = 5, height = 6 )
p2 <- ggplot(Enrich_res,aes(x=Count,y=Description))+
  geom_point(aes(size=Count,color= -log10(pvalue)))+
  theme_bw()+labs(y="",x="Count")+ 
  scale_color_gradient(low = "lightblue", high = "darkblue")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +  # 动态换行
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  
print(p2)      
dev.off()
#scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
#调整Y轴标签角度
