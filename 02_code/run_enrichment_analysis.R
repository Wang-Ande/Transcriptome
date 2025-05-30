run_enrichment_analysis <- function(data, OrgDb, dir) {
  library(clusterProfiler)
  library(pathview)
  library(org.Hs.eg.db)
  library(ggplot2)
  
  # 定义辅助函数
  function1 <- function(OrgDb, genesymbol) {
    if (OrgDb == "Hs") {
      print("OrgDb: Human")
      gene.entrez.eg <- id2eg(ids = genesymbol$gene, category = 'SYMBOL', 
                              org = 'Hs', na.rm = F)
    } 
    
    gene.entrez.eg <- as.data.frame(gene.entrez.eg)
    gene.entrez.eg <- na.omit(gene.entrez.eg)
    # 检查是否有有效的Entrez ID
    if (is.null(gene.entrez.eg) || nrow(gene.entrez.eg) == 0) {
      stop("No valid Entrez IDs found.")
    }
    
    # GO KEGG分析
    if (OrgDb == "Hs") {
      kkd <- tryCatch({
        enrichGO(gene = gene.entrez.eg$ENTREZID,
                      OrgDb = "org.Hs.eg.db",
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)
      }, error = function(e) {
        message("Error in enrichGO: ", e$message)
        return(NULL)
      })
      kk <- tryCatch({
        enrichKEGG(gene = gene.entrez.eg$ENTREZID, 
                   organism = "hsa", 
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)
      }, error = function(e) {
        message("Error in enrichKEGG: ", e$message)
        return(NULL)
      })
    } 
    result <- list(enrichGO = kkd, enrichKEGG = kk)
    return(result)
  }
  
  # 创建目录
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  # down-regulated genes
  down_genes <- subset(data, logFC < -1)
  result_down <- function1(OrgDb, down_genes)
  kkd_down <- result_down$enrichGO
  kk_down <- result_down$enrichKEGG
  
  write.csv(kkd_down@result, file = paste0(dir, "/GO_down.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/GO_down.pdf"), width = 6, height = 7)
  p1 <- dotplot(kkd_down, showCategory = 5, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
  print(p1)
  dev.off()
  
  write.csv(kk_down@result, file = paste0(dir, "/KEGG_down.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/KEGG_down.pdf"), width = 6, height = 5)
  p2 <- dotplot(kk_down,showCategory = 10)
  print(p2)
  dev.off()
  
  # up-regulated genes
  up_genes <- subset(data, logFC > 1)
  result_up <- function1(OrgDb, up_genes)
  kkd_up <- result_up$enrichGO
  kk_up <- result_up$enrichKEGG
  
  write.csv(kkd_up@result, file = paste0(dir, "/GO_up.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/GO_up.pdf"), width = 6, height = 7)
  p3 <- dotplot(kkd_up, showCategory = 5, split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free')
  print(p3)
  dev.off()
  
  write.csv(kk_up@result, file = paste0(dir, "/KEGG_up.csv"), quote = F, row.names = F)
  pdf(file = paste0(dir, "/KEGG_up.pdf"), width = 6, height = 5)
  p4 <- dotplot(kk_up,showCategory = 10)
  print(p4)
  dev.off()
  
  # 统计上下调的通路的数量 ----
  ## 下调GO ----
  down_GO <- table(kkd_down@result$p.adjust<0.05)
  ## 下调KEGG ----
  down_KEGG <- table(kk_down@result$p.adjust<0.05)
  ## 上调GO ----
  up_GO <- table(kkd_up@result$p.adjust<0.05)
  ## 上调KEGG ----
  up_KEGG <- table(kk_up@result$p.adjust<0.05)
  print(list(down_GO,down_KEGG,up_GO,up_KEGG))
  }
