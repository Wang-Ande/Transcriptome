# library paks
rm(list = ls())
library(RcisTarget)
library(dplyr)
library(openxlsx)

# output dir
dir <- "./03_result/08_RcisTarget/OCI/"

# Load gene sets to analyze. e.g.:
DE_res <- read.xlsx("./03_result/02_DE/OCI/High_vs_WT/DE_Gene_Names.xlsx")
head(DE_res)

# Select motif database to use (i.e. organism and distance around TSS)
motifRankings <- importRankings("D:/R/R-Project/tools/Human/database/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
motifRankings 
head(motifRankings@rankings)
# motif annotation
data(motifAnnotations_hgnc)
motifAnnotations[1:3,1:6]

# level gene name
library(org.Hs.eg.db)
library(dplyr)

# 获取所有 protein-coding gene SYMBOL
protein_coding <- keys(org.Hs.eg.db, keytype = "SYMBOL") %>%
  AnnotationDbi::select(org.Hs.eg.db, ., columns = "GENETYPE", keytype = "SYMBOL") %>%
  filter(GENETYPE == "protein-coding") %>%
  pull(SYMBOL)
length(protein_coding)
# 只保留 protein-coding gene
DE_genes_pc <- DE_res$gene[DE_res$gene %in% protein_coding]
# 去掉 NA / 空值
DE_genes_pc <- DE_genes_pc[!is.na(DE_genes_pc)]
DE_res <- DE_res[DE_res$gene%in%DE_genes_pc,]

 # genelist
table(DE_res$Sig)
upgenes <- DE_res %>%
  dplyr::filter(Sig == "up") %>%
  dplyr::arrange(desc(logFC)) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(gene)
downgenes <- DE_res %>%
  dplyr::filter(Sig == "down") %>%
  dplyr::arrange(logFC) %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::pull(gene)
genelists <- list("upgenes"= upgenes,"downgenes"=downgenes)
genelists

# step by step ----
## Calculate enrichment ----
genelists <- c("MTHFD2","PHGDH","PSAT1","SHMT2","SLC1A4","SLC1A5","SLC6A9")
motifs_AUC <- calcAUC(genelists, motifRankings, nCores=8)
motifs_AUC
# AUC for 2 gene-sets and 5876 motifs.
# 
# AUC matrix preview:
#           bergman__Su_H_ bergman__croc bergman__pho bergman__tll c2h2_zfs__M0369
# upgenes       0.04908867     0.0759253            0            0      0.00000000
# downgenes     0.00196802     0.0000000            0            0      0.04194342

# 识别显著富集的motif：AUC值高于阈值的 motif 通常被认为是显著的
par(mfrow = c(1,2))
pdf(paste(dir, "auc.pdf"),width = 9,height = 7)
for(i in c("upgenes",'downgenes')){
  auc <- getAUC(motifs_AUC)[i,]
  hist(auc, main=i, xlab="AUC histogram",
       breaks=100, col = "#0000ff50", border = "darkblue") 
  nes3 <- (3*sd(auc)) + mean(auc) # 计算阈值,阈值设置为均值加上3倍标准差
  abline(v=nes3, col="tomato")
}
dev.off()

# motif-TF annotation ----
# Select significant motifs and/or annotate to TFs
# The motifs are considered significantly enriched if they pass the the Normalized Enrichment Score (NES) threshold.
# 显著性motif的选择是基于归一化富集评分( Normalized Enrichment Score,NES)进行的。
# 每个motif的NES是根据基因集中所有基序的AUC分布[(x-mean)/sd]计算。那些通过给定阈值(默认为3.0)的motifs 被 认为是重要的。
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations,
                                           highlightTFs = NULL)
# highlightTFs=list(upgenes="RFX5",downgenes='ZNF274'))
head(motifEnrichmentTable[,-"TF_lowConf", with=FALSE])
#    geneSet                                                  motif   NES   AUC highlightedTFs TFinDB                         TF_highConf
#     <char>                                                 <char> <num> <num>         <char> <char>                              <char>
# 1: upgenes                                      metacluster_174.5  5.74 0.183           RFX5        ZNF671; ZNF671 (directAnnotation). 
# 2: upgenes                                     metacluster_155.27  5.46 0.175           RFX5           ZFX; ZNF667 (directAnnotation). 
# 3: upgenes            dbtfbs__ZNF660_HEK293_ENCSR283DOU_merged_N1  5.15 0.166           RFX5                ZNF660 (directAnnotation). 
# 4: upgenes                                   transfac_pro__M06089  5.08 0.164           RFX5                 ZNF26 (directAnnotation). 
# 5: upgenes taipale_tf_pairs__TEAD4_CLOCK_NCACGTGNNNNNNNCATWCC_CAP  5.07 0.163           RFX5          CLOCK; TEAD4 (directAnnotation). 
# 6: upgenes                                   transfac_pro__M05767  5.04 0.163           RFX5                 ZNF77 (directAnnotation). 

# Identify the genes with the best enrichment for each Motif ----
motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=genelists)
motifEnrichmentTable_wGenes[1:4,1:4]
#    geneSet                                       motif   NES   AUC
#     <char>                                      <char> <num> <num>
# 1: upgenes                           metacluster_174.5  5.74 0.183
# 2: upgenes                          metacluster_155.27  5.46 0.175
# 3: upgenes dbtfbs__ZNF660_HEK293_ENCSR283DOU_merged_N1  5.15 0.166
# 4: upgenes                        transfac_pro__M06089  5.08 0.164

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })
anotatedTfs

motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
write.csv(
  motifEnrichmentTable_wGenes_wLogo,
  file = paste0(dir, "motif_enrichment_table.csv"),
)
# 可以把后面的[1:10]去掉
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

library(DT)
library(htmlwidgets)
dt <- datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))
saveWidget(dt, file = paste0(dir, "motif_enrichment_table.html"), selfcontained = TRUE)

# ----network
signifMotifNames <- motifEnrichmentTable$motif[1:3]
par(mfrow=c(1,3))
incidenceMatrix <- getSignificantGenes(genelists$upgenes,
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000,
                                       genesFormat="incidMatrix",
                                       method="iCisTarget")$incidMatrix
# 网络图
library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),
                    label=c(motifs, genes),
                    title=c(motifs, genes), # tooltip
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))

visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE,
                                        nodesIdSelection = TRUE) %>%
  visExport(type = "pdf", name = "network_image", label = "Export as PDF") # 添加导出按钮
