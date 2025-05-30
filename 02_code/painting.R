rm(ls())
setwd("")

# Packages
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)
library(openxlsx)

# 1. Enrich plot ---- 
# DE_gene input 
DE_gene <- read.csv("./03_Result/DE/OCI_AML2/OCI_VEN_vs_OCI_WT/result_DE.csv")
DE_gene <- DE_gene[DE_gene$Sig != "stable" ,]

# 转换基因名 
y <- DE_gene$Genes
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
DE_gene$gene <- gene

# 设置数据库 
GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
KEGG_database <- 'hsa'         # KEGG是hsa数据库

# gene ID转换 
gene <- bitr(DE_gene$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)

# GO 
# GO富集分析
kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                keyType = "ENTREZID", # 设定读取的gene ID类型
                ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2)# 设定q值阈值

# KEGG 
# KEGG富集分析
kk <- enrichKEGG(gene = gene$ENTREZID,
                 keyType = "kegg",
                 organism = KEGG_database,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

#GO、KEGG结果整合 
result <- list(enrichGO = kkd, enrichKEGG = kk)
GO_res <- result$enrichGO
KEGG_res <- result$enrichKEGG

# res_output 
dir_enrich <- "./03_Result/GO&KEGG/OCI_AML2/6W_vs_WT/"

# 导出enrichGO 
write.csv(GO_res@result, file = paste0(dir_enrich, "/GO_all.csv"), quote = F, row.names = F)

# 导出enrichKEGG
write.csv(KEGG_res@result, file = paste0(dir_enrich, "/KEGG_all.csv"), quote = F, row.names = F)

## 1.1 dotplot ----
# ggplot2画气泡图，scale_color_gradient设置蓝红配色
library(ggplot2)
library(stringr)

# Enrich_res input 
Enrich_res <- read.xlsx("./03_Result/GO&KEGG/OCI_AML2/6W_vs_WT/KEGG_all_selectd.xlsx")
Enrich_res <- Enrich_res[Enrich_res$p.adjust<0.05,]
Enrich_res <- Enrich_res[order(Enrich_res$Count),]

pdf(file ="./03_Result/GO&KEGG/OCI_AML2/6W_vs_WT/KEGG_all.pdf",
    width = 5, height = 4 )
p2 <- ggplot(Enrich_res,aes(x=Count,y=Description))+
      geom_point(aes(size=Count,color= -log10(p.adjust)))+
      theme_bw()+labs(y="",x="Count")+ 
      scale_color_gradient(low = "lightblue", high = "darkblue")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +  # 动态换行
      theme(axis.text.y = element_text(angle = 0, hjust = 1))  
print(p2)      
dev.off()
             #scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
             #调整Y轴标签角度

## 1.2 ridge plot ----
rm(list = ls()) 

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggridges)
library(openxlsx)  
library(org.Hs.eg.db)
library(clusterProfiler)

# output category 
dir_gsea <- "./03_Result/GSEA/OCI_M2/VEN_VS_WT/"

# 加载数据
# 包含logfc的genelist 
data_gene <- read.csv("./03_Result/DE/OCI_AML2/OCI_VEN_vs_OCI_WT/result_DE.csv")

if(T){
# 提取基因名（symbol） ----
y <- data_gene$Genes
SYMBOL <- unlist(lapply(y,
                        function(y) strsplit(as.character(y),";")[[1]][1]))
id <- data.frame(SYMBOL)
diff_gene <- cbind(id,data_gene)
# 检查有无重复symbol
duplicates <- duplicated(SYMBOL)
a <- print(diff_gene[duplicates,])

#如果有则进行以下处理
#计算行平均值，按降序排列!!!
matrix0 <- order(diff_gene$AveExpr,decreasing = T)

#调整表达谱的基因顺序
expr_ordered=diff_gene[matrix0,]

#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(expr_ordered$SYMBOL)#$Name是指基因名那一列

#得到最后处理之后的表达谱矩阵
diff_gene=expr_ordered[keep,]
}

# genelist 
diff_gene_sort <- diff_gene %>% 
  arrange(desc(logFC))
genelist_df <- diff_gene_sort[,c("SYMBOL","logFC")]

# 输入gsea结果
gsearesult <- read.xlsx("./03_Result/GSEA/OCI_M2/VEN_VS_WT/KEGG_results_selected.xlsx")

# 按p.adjust排列，只显示前20个
gsearesult <- gsearesult %>%
              dplyr::arrange(p.adjust) %>%
              head(26)

# 将 gene_df 和 gsearesult 合并
gsearesult <- gsearesult %>%
  dplyr::mutate(core_enrichment = as.character(core_enrichment)) %>%
  separate_rows(core_enrichment, sep = "/")  # 如果 core_enrichment 中是以 / 分隔的基因集合

# 合并 logFC 信息
gsearesult <- gsearesult %>%
  left_join(genelist_df, by = c("core_enrichment" = "SYMBOL"))


# 处理 GSEA 结果
gsearesult <- gsearesult %>%
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  dplyr::mutate(log10_adj_p = -log10(p.adjust)) %>%
  separate_rows(core_enrichment, sep = "/")

# 合并基因信息
gsearesult <- gsearesult %>%
  left_join(genelist_df, by = c("core_enrichment" = "SYMBOL"))

# 绘制 GSEA 图
ridge_plot <- ggplot(gsearesult, aes(x = logFC.x, y = Description, fill = log10_adj_p)) +
              geom_density_ridges(alpha = 0.7, scale = 1.5) +                       
              geom_point(aes(size = abs(NES), x = min(logFC.x)*1.2, color = NES)) +  
              # 加入点图，按 NES 的大小进行点的缩放,并将点的位置稍微外扩，放在与密度曲线重叠
              scale_size_continuous(range = c(1, 6), name = 'NES') +                   # 设定点大小的范围
              labs(x = "Log2 Fold Change", y = "Pathway", title = "GSEA Enrichment Analysis") +  # 设置标题和轴标签
              scale_fill_distiller(palette = 'Spectral', name = "-Log10(p.adj)") +    # 色带按 -log10(p.adjust) 调色
              scale_color_distiller(palette = 'RdBu', name = "NES", direction = -1, limits = c(-3, 3)) +  # NES 的颜色调节
              theme_minimal() +  # 使用简洁的主题
              theme(
              panel.background = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
              axis.text = element_text(color = "black", size = 12),
              axis.ticks = element_blank(),
              axis.text.y = element_text(color = "black", size = 11),
              axis.title = element_text(size = 12),
              axis.title.x = element_text(color = "black", size = 14),
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5)         # 标题加粗居中
                   ) +
              guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))# 设置图例

ggsave(filename = paste0(dir_gsea, "GSEA_Top20.pdf"), plot = ridge_plot,
       device = "pdf",  # 指定保存为 pdf 格式
       dpi = 300)
