# 1. Library packages ----
library(readxl)
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(tidyr)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(pathview)
library(openxlsx)
library(enrichplot)
library(patchwork)
library(KEGG.db) # 本地富集必须
# library(gen)
source("02_code/run_limma_DE.R")                 # load前检查是否需要更改！
source("02_code/run_GoKegg_Enrich.R")

# 2. DE_limma ----

## 2.1 Expr input ----
exprSet <- read.csv("./01_data/gene_count_matrix.csv", row.names = 1)
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[,-grep('OCI|MOLM13',colnames(exprSet))]
## OCI 细胞系去除4w
exprSet <- exprSet[,-4]

## 2.3 Filter low genes ----
# 去除非编码蛋白
source("../../Code_Box/Transcriptome/run_filter_coding_gene.R")
exprSet <- filter_coding_genes(exprSet)
exprSet <- exprSet[,!colnames(exprSet)%in%c("ense","symbol")]

# approach 1
min_sample <- ceiling(length(colnames(exprSet)) / 2)  # 向上取整
exprSet_filtered <- exprSet[rowSums(exprSet > 10) >= min_sample,]  
dim(exprSet_filtered)

# approach 2
exprSet_log2 <- log2(exprSet+1)
exprSet_filtered <- exprSet_log2[apply(exprSet_log2,1,mean)>1,] # 13454
exprSet_filtered <- 2^exprSet_filtered-1

## 2.4 Group input ----
colData <- read_xlsx("./01_data/group_info.xlsx")
colData <- as.data.frame(colData)
colnames(colData)
colData <- colData[-grep('OCI|MOLM13',colData$id),-3]
## OCI 细胞系去除4w
colData <- colData[-4,]

# 整理样本
rownames(colData) <- colData$id
table(colData$group) 

## 2.5 Group design ----
group_1 <- "High"    # group 1为实验组
group_2 <- "WT"     # group 2为对照组
targeted_group <- colData
targeted_data <- as.data.frame(exprSet_filtered) 

## 2.6 Res output ---- 
source("02_code/run_limma_DE.R")                 # load前检查是否需要更改！
result_merge <- run_limma_DE(exprSet = targeted_data,
                             colData = targeted_group,
                             group_1 = group_1,
                             group_2 = group_2,
                             log2 = FALSE,       # 若使用Voom方法，则不对count进行log2标准化
                             logfc_threshold = 1,        # "1"、"0.585"、"0.263"
                             pvalue_threshold = 0.05, 
                             qvalue_threshold = NULL,
                             dir = "03_result/02_DE/MV4")
table(result_merge$DE_res$Sig)
# 导出gene symbol
DE_Genes <- read.csv('./03_result/02_DE/MV4/High_vs_WT/DE.csv',row.names = 1)
y <- DE_Genes$Row.names
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
DE_Genes$gene <- gene
DE_Genes <- DE_Genes[,c("gene","logFC","P.Value","adj.P.Val","Sig")]
write.xlsx(DE_Genes, file = "./03_result/02_DE/MV4/High_vs_WT/DE_Gene_Names.xlsx")

# 3. GO KEGG --------------------------------------------------------------

## 3.1 Set path ----
DE_dir <- './03_result/02_DE/MOLM13/'           # 差异分析结果目录
Enrich_dir <- './03_result/03_Enrichment/MOLM13/'  # 富集结果存放目录
log_file <- file.path(Enrich_dir, "enrichment_fail_log.txt")  # 失败日志文件

GO_database <- 'org.Hs.eg.db'  # GO数据库
KEGG_database <- 'hsa'         # KEGG数据库

## 3.2 DE_RES  ---- 
DE_files <- list.files(DE_dir, pattern = "DE\\.csv$", 
                       full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

## 3.3 Circulation for all DE_res ----
for (de_file in DE_files) {
  
  # 样本对名称，用子文件夹名作为标签
  sample_tag <- basename(dirname(de_file))
  
  # 创建输出子文件夹
  out_dir <- file.path(Enrich_dir, sample_tag)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  ## 读取 DE 文件
  DE_result <- read.csv(de_file)
  
  # 提取 gene symbol（假设 Row.names 列有 ENSEMBL|symbol）
  DE_result$gene <- sapply(strsplit(DE_result$Row.names, "\\|"), `[`, 2)
  
  # 构建 genelist
  genelist <- DE_result[, c("gene", "logFC", "P.Value", "adj.P.Val")]
  
  source("./02_code/run_GoKegg_Enrich.R")
  # 使用 tryCatch 运行富集分析，并记录失败日志
  tryCatch({
    Enrich_Res <- run_GoKegg_Enrich(Genelist = genelist,
                                    pvalue = 0.05,
                                    logfc = 0.5,
                                    go_db = GO_database,
                                    kegg_db = KEGG_database,
                                    Output_dir = out_dir)
    message(sample_tag, " enrichment analysis finished.")
  }, error = function(e) {
    # 打印错误信息
    message(sample_tag, " enrichment analysis failed: ", e$message)
    # 将失败信息写入日志文件（追加模式）
    write(paste(Sys.time(), sample_tag, e$message, sep = "\t"), 
          file = log_file, append = TRUE)
  })
}
