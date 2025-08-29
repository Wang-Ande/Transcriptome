run_GoKegg_Enrich <- function(Genelist,
                              Output_dir,
                              pvalue = 0.05,
                              logfc = 1,
                              go_db = "org.Hs.eg.db" ,
                              kegg_db = 'hsa') {
  library(clusterProfiler)
  library(ggplot2)
  library(openxlsx)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(KEGG.db)
  
  # 检查列名
  required_cols <- c("gene", "logFC", "P.Value")
  if (!all(required_cols %in% colnames(Genelist))) {
    stop("Input Genelist must contain columns: gene, logFC, P.Value")
  }
  
  # 创建输出目录
  dir.create(Output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 按上下调分组
  message("Filtering gene sets...")
  up_genes <- subset(Genelist, logFC > logfc & P.Value < pvalue)
  down_genes <- subset(Genelist, logFC < -logfc & P.Value < pvalue)
  all_genes <- subset(Genelist, abs(logFC) > logfc & P.Value < pvalue)
  
  # 内部函数
  do_enrichment <- function(df, tag) {
    gene_entrez <- tryCatch({
      bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = go_db)
    }, error = function(e) return(NULL))
    
    if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {
      warning(paste0("No ENTREZID found for ", tag))
      return(NULL)
    }
    
    enrich_go <- enrichGO(gene = gene_entrez$ENTREZID,
                          OrgDb = go_db,
                          keyType = "ENTREZID",
                          ont = "ALL",
                          pvalueCutoff = 0.5,
                          qvalueCutoff = 1,
                          readable = TRUE)
    
    enrich_kegg <- enrichKEGG(gene = gene_entrez$ENTREZID,
                              organism = kegg_db,
                              keyType = "kegg",
                              pvalueCutoff = 0.5,
                              qvalueCutoff = 1,
                              use_internal_data = TRUE)
    enrich_kegg <- setReadable(enrich_kegg, OrgDb = go_db, keyType = "ENTREZID")
    
    # 写出 GO 结果
    write.xlsx(enrich_go@result, file = file.path(Output_dir, paste0("GO_", tag, ".xlsx")), rowNames = FALSE)
    pdf(file = file.path(Output_dir, paste0("GO_", tag, ".pdf")), width = 5.5, height = 7)
    print(dotplot(enrich_go, showCategory = 5, color = "pvalue", split = "ONTOLOGY") +
            facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free'))
    dev.off()
    
    # 写出 KEGG 结果
    write.xlsx(enrich_kegg@result, file = file.path(Output_dir, paste0("KEGG_", tag, ".xlsx")), rowNames = FALSE)
    pdf(file = file.path(Output_dir, paste0("KEGG_", tag, ".pdf")), width = 6, height = 5)
    print(dotplot(enrich_kegg, showCategory = 10, color = "pvalue"))
    dev.off()
    
    return(list(GO = enrich_go, KEGG = enrich_kegg))
  }
  
  message("Running GO/KEGG for downregulated genes...")
  result_down <- do_enrichment(down_genes, "down")
  
  message("Running GO/KEGG for upregulated genes...")
  result_up <- do_enrichment(up_genes, "up")
  
  message("Running GO/KEGG for all DE genes...")
  result_all <- do_enrichment(all_genes, "all_DE")
  
  return(list(
    down = result_down,
    up = result_up,
    all = result_all
  ))
}
