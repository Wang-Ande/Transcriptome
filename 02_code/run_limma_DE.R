run_limma_DE <- function(exprSet, colData, group_1, group_2, log2,
                         logfc_threshold, pvalue_threshold, qvalue_threshold = NULL,
                         dir = getwd()) {
  library(limma)
  library(edgeR)
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(cowplot)
  
  # 创建输出目录
  output_dir <- paste0(dir, "/", group_1, "_vs_", group_2)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 提取group_1和group_2的样本
  group_1_samples <- colData$id[colData$group == group_1]
  group_2_samples <- colData$id[colData$group == group_2]
  target_group <- colData[c(group_1_samples, group_2_samples),]
  
  # 检查是否为整数型（counts）
  if (!all(exprSet == floor(exprSet))) {
    warning("Warning: The exprSet matrix contains non-integer values. Make sure you are inputting raw counts for voom.")
  }
  # 确保样本顺序与exprSet列顺序一致
  target_data <- exprSet[, c(group_1_samples, group_2_samples)]
  
  # check data: 检查exprSet和colData的列名与行名是否匹配
  if (!identical(sort(colnames(target_data)), sort(rownames(target_group)))) {
    stop("Error: The col names of 'target_data' and the row names of 'target_group' do not match.")
  }
  
  # 分组因子化：确保分组因子正确
  target_group$group <- factor(target_group$group, levels = c(group_1, group_2))
  
  
  # 可选log2转换
  if (log2) {
    target_data <- log2(target_data + 1)  # log2处理，避免0值
  }
  
  # 设计矩阵
  design <- model.matrix(~0 + factor(target_group$group))
  colnames(design) = levels(factor(target_group$group))
  rownames(design) = colnames(target_data)
  
  # dge对象创建 (若log2为true，则更改！)
  dge <- DGEList(counts = target_data, group = target_group$group)
  
  # filterByExpr()去除low genes 
  # keep <- filterByExpr(dge)
  # dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # 归一化 
  dge <- calcNormFactors(dge)
  
  # voom 方法标准化 
  pdf(file = paste0(output_dir, "/voom_plot.pdf"), width = 8, height = 6)
  v <- voom(dge, design, plot = TRUE, normalize = "none") # 若噪声过大，则设置normalize = "quantile" （一般不用防止矫枉过正）
  dev.off()
  # 如果没有原始count数据，是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化！！！！！
  # 可直接从这里开始进行差异分析
  #v <- target_data
  
  # 使用线性模型进行拟合
  fit <- lmFit(v, design)
  
  # 指定对比组别 
  con <- paste(levels(target_group$group), collapse = "-")
  # [1] "tumor-normal"
  
  # 创建对比矩阵 
  cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2,trend = TRUE)
  
  # 获取差异结果 
  tempOutput <- topTable(fit2, coef = con, n = Inf)
  DE_res <- na.omit(tempOutput)
  
  # 加change列 ----
  # 标记上下调基因，可根据需求设定阈值
  logFC = logfc_threshold
  P.Value = pvalue_threshold
  k1 <- (DE_res$P.Value < P.Value) & (DE_res$logFC < -logFC)
  k2 <- (DE_res$P.Value < P.Value) & (DE_res$logFC > logFC)
  DE_res <- mutate(DE_res, 
                   Sig = ifelse(k1, "down", 
                                ifelse(k2, "up", "stable")))
  table(DE_res$Sig)
  # down stable     up 
  #  749  25664    378 
  result_merge <- merge(target_data, DE_res, by = "row.names", all.y = TRUE)
  
  write.csv(result_merge,file = paste0(output_dir,"/DE.csv"))
  
  # 火山图 ----
  # change列因子化 ----
  DE_res$Sig <- factor(
    DE_res$Sig,
    levels = c("up", "down", "stable"))  # 强制按此顺序排列
  
  # 计算对称的横坐标范围（基于logFC绝对值的最大值）
  max_abs_logfc <- max(abs(DE_res$logFC), na.rm = TRUE)
  
  # 扩展5%的余量，避免点紧贴坐标轴边缘
  x_limit <- max_abs_logfc * 1.05
  
  # plot ----
  p <- ggplot(data = DE_res, 
              aes(x = logFC, 
                  y = -log10(P.Value))) +
    geom_point(alpha = 0.5, size = 2, 
               aes(color = Sig)) +
    ylab("-log10(Pvalue)")+
    # 按因子顺序指定颜色
    scale_color_manual(
      name = "Change",                              # 图例标题
      values = c(
        "up" = "#B30000",      # 红
        "stable" = "grey",     # 灰
        "down" = "#003366"),    # 蓝,
      labels = c(                                  # 显示上下调的个数
        "up" = paste0("Up: ", sum(DE_res$Sig == "up")),
        "stable" = paste0("Stable: ", sum(DE_res$Sig == "stable")),
        "down" = paste0("Down: ", sum(DE_res$Sig == "down"))))+
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), lty = 4, 
               col = "black", lwd = 0.8,alpha=0.4) +
    geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, 
               col = "black", lwd = 0.8,alpha=0.4) +
    labs(title = paste0(group_1,"-",group_2)) +
    xlim(-x_limit, x_limit)+
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),  # 标题居中
          aspect.ratio = 1.2)  # 设置纵横比，调整为更高
  
  ggsave(filename = paste0(output_dir,"/volc.pdf"),
         plot = p, device = "pdf", 
         width = 6, height = 5)
  
  
  # 热图 
  return(list(result_merge = result_merge, DE_res = DE_res))
}
