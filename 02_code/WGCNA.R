# WGCNA

# 加载包 ----
if(!require(biomaRt ,quietly = TRUE)){BiocManager::install("biomaRt")}
library(biomaRt)
library(readr)
library(openxlsx)
library(readxl)
library(dplyr)
library(scales)
library(tidyr)
library(WGCNA)
library(forcats)
library(ggplot2)
library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)

# set output category ----
dir_WGCNA <- "./03_result/WGCNA/CPM_All_gene/"

# 1.1 在线获取基因长度方式 ----
#若无法从源文件获取TPM，则通过下载在线基因长度数据来计算TPM

## 1.1.1 数据预处理 ----

## 1.1.2 标准化 TPM ----
library(biomaRt)

## 1.1.3 获取基因长度 ----
# 连接到 Ensembl 数据库 (可以根据需要修改物种)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") 

# 获取基因的起始和结束位置
gene_lengths <- getBM(attributes = c("ensembl_gene_id", 
                                     "start_position", 
                                     "end_position"),
                      mart = ensembl)

# 计算基因长度 (结束位置 - 起始位置 + 1)
gene_lengths$length <- gene_lengths$end_position - gene_lengths$start_position + 1

# 查看结果
head(gene_lengths)

## 1.1.4 读入counts数据 ----
expr_counts <- read.csv("./01_data/Counts_data/data_merge_adjusted.csv")

# 提取ensemble名
y <- expr_counts$X
ensemble <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][1]))
id <- data.frame(ensemble)

# 去除版本号(即去除小数点后的内容，ENSG00000005194.15)
id$ensemble <- gsub("\\.\\d+$", "", id$ensemble)  # 去除版本号

## 1.1.5 合并基因长度数据 ----
expr_counts <- cbind(id,expr_counts)
expr_counts <- merge(expr_counts, gene_lengths, by.x = "ensemble", 
                     by.y = "ensembl_gene_id", all.x = TRUE)

# 查看合并后的数据
head(expr_counts)

# 使用基因长度平均值填补缺失值
average_length <- mean(gene_lengths$length, na.rm = TRUE)
expr_counts$length[is.na(expr_counts$length)] <- average_length

# 计算每个基因的 RPK
rpk <- apply(expr_counts[,3:29], 2, 
             function(x) x /(expr_counts$length/1000) )  
rpk <- as.data.frame(rpk)

# 计算所有基因的 RPK 总和
total_rpk <- apply(rpk,2,sum)

# 计算每个基因的 TPM
tpm <- apply(rpk, 1, function(x) (x / total_rpk) * 1e6)
tpm <- t(tpm)

# 合并行名
gene_id <- expr_counts$X
gene_id <- as.data.frame(gene_id)
datExpr <- cbind(gene_id,tpm)

## 1.1.6 保留tpm结果
write.csv(datExpr,file = "./01_data/Counts_data/data_merge_adjusted_tpm.csv")

# 1.2 上游分析获取基因长度方式 ----
# 差异基因取并集
DEG1 <- read.csv("./03_result/DE_limma/MOLM13_2W_vs_MOLM13_WT/DE.csv")
DEG2 <- read.csv("./03_result/DE_limma/MOLM13_6W_vs_MOLM13_WT/DE.csv")
DEG3 <- read.csv("./03_result/DE_limma/MV4_11_2W_vs_MV4_11_WT/DE.csv")
DEG4 <- read.csv("./03_result/DE_limma/MV4_11_6W_vs_MV4_11_WT/DE.csv")
DEG5 <- read.csv("./03_result/DE_limma/OCI_AML2_2W_vs_OCI_AML2_WT/DE.csv")
DEG6 <- read.csv("./03_result/DE_limma/OCI_AML2_6W_vs_OCI_AML2_WT/DE.csv")

# 将它们存储在一个列表中
DEG_list <- list(DEG1, DEG2, DEG3, DEG4, DEG5, DEG6)  

# 初始化一个空的基因 ID 列表
union_genes <- character(0)

# 循环遍历 DEG_list 中的每个 DEG 数据框
for (DEG_df in DEG_list) {
    # 提取当前 DEG 数据框中的基因 ID 列，并排除 'stable' 样本
    DEG_genes <- DEG_df[DEG_df$Sig != "stable", 2]
    
    # 将当前基因 ID 列与 union_genes 的并集赋值给 union_genes
    union_genes <- union(union_genes, DEG_genes)
}


# 若从上游测序跑TPM，则从此步开始 
## 1.2.2 数据预处理 ----
# 添加偏移量，消除负值
# expr_input
count_raw <- read.csv("./01_data/Cpm_data/data_merged_adjusted.csv")
rownames(count_raw) <- count_raw$X
count_raw <- count_raw[,-1]
# 去除样本名多余信息
colnames(count_raw) <- gsub("cas9","",colnames(count_raw))
colnames(count_raw) <- gsub("w","W",colnames(count_raw))
count_DE <- count_raw[union_genes,]
min(count_DE)

datExpr <- count_raw

# 去除4w样本
datExpr <- datExpr[,-grep("4w",colnames(datExpr))]

## 1.2.3 标准化 log2+1 ----
datExpr <- log2(datExpr + 1)

## 1.2.4 filter ---- 
# WGCNA可以不过滤，WGCNA包中有自己的过滤方案
# approach 1
# var计算每个基因方差，筛选基因变化较大的基因，此次选取前75%的基因
vars_res <- apply(datExpr, 1, var)

# 计算百分位数截止值
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) # quantile生成分位数
per_res

upperGene <- datExpr[ which(vars_res > per_res[4]),]  # 选取方差位于前75%的基因
dim(upperGene)
datExpr <- data.matrix(upperGene)

# approach 2
datExpr <- datExpr[apply(datExpr,1,mean)>1,] # 13454

# 检查是否有目的基因被过滤掉
targetgene_rows <- grep("BAK1", rownames(datExpr))
targetgene <- datExpr[targetgene_rows,]
# TP53\BCL2\MCL1\BAX\BAK1基因都存在

# 2 t转置 ----
# 为了执行 WGCNA，通常需要转置数据，使其符合 WGCNA 的标准格式
datExpr <- t(datExpr)                       # 让每一列代表一个基因

# 3. 检查数据好坏 ----
gsg <- goodSamplesGenes(datExpr, verbose = 3)

# 查看标记结果
gsg$allOK
#[1] TRUE

# 如果为false则运行下段
if (!gsg$allOK) {
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", 
                         paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", 
                         paste(names(datExpr)[!gsg$goodSamples], collapse = ", ")))
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

ncol(datExpr)

# 对样本聚类查看有无异常样本
sampleTree = hclust(dist(datExpr), method = "average")
png(paste0(dir_WGCNA , "sample_clustering.png"), 
    width = 1500, 
    height = 1000,
    res = 300)
par(cex = 0.6)
par(mar = c(0,4,2,0)) 
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5, 
     cex.axis = 1.5, 
     cex.main = 2)
#abline(h = 150, col = 'red')
dev.off()

# 剔除离群样本，无则跳过
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10) 
table(clust)               #显示剔除多少样本， 0代表剔除的个数，1代表保留的个数
datExpr = datExpr[clust == 1, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 4. 选择软阈值 ----
# 定义1:30的power值，用于构建不同的基因共表达网络
powers <- c(c(1:10), 
            seq(from = 12, 
                to = 30,
                by = 2))

sft <-  pickSoftThreshold(datExpr,
                          powerVector = powers,
                          verbose = 5,
                          networkType = 'unsigned')
sft$powerEstimate
# 绘制拟合度图来选择合适的软阈值
pdf(file = paste0(dir_WGCNA , "sft_par.pdf"),width = 10,height = 6.5)
#sizeGrWindow(9, 5) 
par(mfrow = c(1,2))
#cex1 = 0.85
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     col="steelblue")
abline(h=0.8,col="red")
plot(sft$fitIndices[,1],  sft$fitIndices[,5], 
     type="n", 
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,  col="steelblue")
dev.off()

# 使用合适的软阈值（没有合适就选6）
power <- 20

# 5. 构建共表达矩阵 ----
net = blockwiseModules(datExpr,
                       power = power,           # 软阀值选为6
                       TOMType = "unsigned",    # 构建无尺度网络
                       minModuleSize = 30,      # 最小模块基因数为30
                       reassignThreshold = 0,   # 重分配阈值          
                       mergeCutHeight = 0.25,   # 模块合并阀值
                       numericLabels = F,       # 返回字符标签（如模块颜色名称）
                       pamRespectsDendro = FALSE,         
                       saveTOMs = FALSE,         # 不保存TOM矩阵
                       verbose = 3,
                       maxBlockSize = 20000)    # 可处理的最大模块基因数
# 显示所有模块个数和各个模块中基因数量
table(net$colors)

## 5.1 聚类分析 ----
# 使用层次聚类
geneTree = net$dendrograms[[1]] 
geneTree

png(paste0(dir_WGCNA ,"Module_colors_clusterdendrogram.png"),
    width = 2000,
    height = 1500,
    res = 300,
    type = "cairo")
moduleColors = net$colors
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",             # 标题
                    dendroLabels = FALSE,        # 不显示基因标签
                    hang = 0.03,                 # 树枝悬挂高度（相对于根部）
                    addGuide = TRUE,             # 颜色引导条
                    guideHang = 0.05)            # 颜色引导条与树状图悬挂距离
dev.off()

## 5.2 计算模块特征基因并绘制模块特征基因的层次聚类树 ----
# 计算模块特征基因（Module Eigengenes）
MEList = moduleEigengenes(datExpr, colors = net$colors)
MEs_1 = MEList$eigengenes

# 计算模块特征基因之间的相似性
MEDiss = 1 - cor(MEs_1)

# 使用层次聚类计算模块特征基因的相似性
METree = hclust(as.dist(MEDiss), method = "average")

# 绘制模块特征基因的树状图
pdf(file = paste0(dir_WGCNA, "Module_Eigengene_Clustering.pdf"), width = 7, height = 6)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
# 添加剪切高度
MEDissThres = 0.35  # 剪切高度
abline(h = MEDissThres, col = "red")
dev.off()

# 这里的MEDissThres是剪切模块的高度值，决定了哪些模块应该合并。

## 5.3 合并模块 ----
# 合并模块
merge = mergeCloseModules(datExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)

# 获取合并后的模块颜色
mergedColors = merge$colors

# 获取合并后的模块特征基因
mergedMEs = merge$newMEs

# 显示合并后模块的数量和颜色
table(mergedColors)

# 绘制合并后的模块树状图
pdf(file = paste0(dir_WGCNA,"merged_dynamic.pdf"), width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(moduleColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 6. 模块与表型的关系 ----

### 6.1 表型因子化 ----

trait <- read.csv("./01_data/Sample_info/trait.csv")
rownames(trait) <- trait$X
trait <- trait[,-1]
rownames(trait) <- gsub("cas9","",rownames(trait))
rownames(trait) <- gsub("w","W",rownames(trait))

if(T){
# 模块特征
MEs0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  # 计算模块特征向量
MEs = orderMEs(MEs0)                                       # 对模块特征向量排序

# 计算模块特征向量与表型的相关系数矩阵
moduleTraitCor <- cor(MEs,trait,use = "p")  

# 计算相关系数矩阵的p值
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  

# 构建绘图时用的文本矩阵
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")  

# 修改文本矩阵的维度，与相关系数矩阵相同
dim(textMatrix)=dim(moduleTraitCor)  
}

### 6.2 相关性热图 ----
if(T){
pdf(file = paste0(dir_WGCNA , "Module_Trait labeledHeatmap.pdf"),
    width = 12,height = 8)

# mar（）分别代表图形边距的底部、左侧、顶部和右侧的边距
par(mar = c(7, 7, 2, 2)) 
labeledHeatmap(Matrix = moduleTraitCor,        # 绘制带标签的热图
               xLabels = colnames(trait),     # x轴标签
               yLabels = names(MEs),           # y轴标签
               ySymbols = names(MEs),          # y轴符号
               colorLabels = FALSE,            # 不显示颜色标签
               colors = blueWhiteRed(50),      # 颜色范围
               textMatrix = textMatrix,        # 显示文本矩阵
               setStdMargins = FALSE,          # 不设置标准边距
               cex.text = 0.5,                 # 文本大小
               cex.lab.x = 0.7,                # X轴文本大小
               zlim = c(-1,1),                 # 颜色映射范围
               main = paste("Module-trait relationships"))  # 绘图标题
dev.off()
}

### 6.3 基因与性状和重要模块的关系：基因重要性和模块成员 ----
IC50 = as.data.frame(trait$IC50)
names(IC50) = "IC50"
modNames = substring(names(MEs), 3)                              # 提取模块名称
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")      # 修改列名
names(MMPvalue) = paste("p.MM", modNames, sep="")                # 修改列名

# 基因显著性
geneTraitSignificance = as.data.frame(cor(datExpr, IC50, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", names(IC50), sep="")
names(GSPvalue) = paste("p.GS.", names(IC50), sep="")

### 6.4 模内分析：鉴定具有高GS和MM的基因 ----
module = "greenyellow"                      # 选择模块
column = match(module, modNames)     # 匹配模块名称
moduleGenes = moduleColors==module   # 提取模块内基因
table(moduleGenes)

pdf(paste0(dir_WGCNA, "Greenyellow_MM.vs.GS.pdf"),
    width = 6,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for IC50",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")
dev.off()

### 6.5 提取指定模块的基因做PPI网络分析 ----

adjacency = adjacency(datExpr, power =6)  # 计算邻接矩阵
TOM = TOMsimilarity(adjacency)             # 计算拓扑重叠矩阵（TOM）

module="royalblue"                               # 选择要导出的模块
probes = colnames(datExpr)                  # 获取基因名称
inModule = (moduleColors==module)           # 找到属于当前模块的基因
modProbes=probes[inModule]                  # 提取属于当前模块的基因名称
head(modProbes)                             # 显示基因名称前几行
modTOM = TOM[inModule,inModule]             # 提取属于当前模块的基因之间的TOM值
dimnames(modTOM)=list(modProbes,modProbes)  # 修改维度名称

# 绘制模块基因表达量箱线图
library(reshape2)
library(ggplot2)

datExpr_royalblue <- datExpr[,modProbes]                   #提取基因模块内基因 
datExpr_royalblue <- t(datExpr_royalblue)

# 添加分组信息
group <- traitData[,-3]
colnames(group)[2] <- "group"
colnames(group)[1] <- "id"
group$id <- rownames(group)
# 计算起始和结束位置
start = nchar(group$id) - 3  # 倒数第四的位置
end = nchar(group$id) - 2    # 倒数第二的位置
# 提取子字符串
result = substr(group$id, start, end)
print(result)
group$group <- result

table(group$group)
value_colour <- c("WT" = "#E64B35FF",# Experimental group
                  "2w" = "#4DBBD5FA",# other group1
                  "6w" = "#F2A200")# other group2
# 提取不同细胞系的时间分组
group_molm13 <- group[grepl("MOLM13",group$id),]
group_mv4_11 <- group[grepl("MV4_11",group$id),]
group_oci <- group[grepl("OCI_AML2",group$id),]

pdf(file = "./03_result/WGCNA/Count_ALLgene/QC_boxplot_oci.pdf",
    width = 6,
    height = 4)
QC_boxplot(2^datExpr_royalblue,data_group = group_oci,
           value_colour = value_colour,title = paste("Module",module))
dev.off()


# 这里只选了top100的基因
nTop=165                                       # 设置要选择的基因数目
IMConn = softConnectivity(datExpr[,modProbes])  # 计算当前模块中基因之间的相似性
top=(rank(-IMConn)<=nTop)                      # 找到相似性排名前nTop的基因
filterTOM=modTOM[top,top]                      # 提取相似性排名前nTop的基因之间的TOM值
# for visANT
vis = exportNetworkToVisANT(filterTOM,
                            file = paste("visANTinput-",module,".txt",sep = ""),
                            weighted = T,threshold = 0)  

# for cytoscape
cyt = exportNetworkToCytoscape(filterTOM,
                               edgeFile = paste("./03_result/WGCNA/Count/CytoscapeInput-edges-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               nodeFile = paste("./03_result/WGCNA/Count/CytoscapeInput-nodes-", 
                                                paste(module, collapse="-"), 
                                                ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               altNodeNames = 
                               nodeNames = modProbes[top], 
                               nodeAttr = moduleColors[inModule][top])  


# 7. KEGG GO -------------------------------------------------------------------
## 7.1 Set output catagory----
#OCI_AML2/MV4_11/MOLM13
# 指定文件夹路径
dir.create("./03_Result/GO&KEGG/MOLM13/VEN_vs_WT/")
dir_WGCNA <- "./03_result/WGCNA/CPM_All_gene/"

df <- data.frame(gene = colnames(datExpr), module = moduleColors)
write.csv(df,file = paste0(dir_WGCNA,"module_gene.csv"))
table(df$module)

## 7.2 WGCNA_res input ----
WGCNA_gene <- read.csv('./03_result/WGCNA/CPM_All_gene/module_gene.csv')

# 转换基因名 
y <- WGCNA_gene$gene
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
WGCNA_gene$gene <- gene

## 7.3 Module genes ----
Module_genes <- subset(WGCNA_gene, WGCNA_gene$module == "greenyellow")

if(T){
    # 设置数据库 
    GO_database <- 'org.Hs.eg.db'  # GO是org.Hs.eg.db数据库
    KEGG_database <- 'hsa'         # KEGG是hsa数据库
    
    # gene ID转换 
    gene <- bitr(Module_genes$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = GO_database)
    
    ## 7.4 GO ----
    # GO富集分析
    kkd <- enrichGO(gene = gene$ENTREZID, # 导入基因的ENTREZID编号
                    OrgDb = GO_database, # 用到的数据库（人类是：org.Hs.eg.db）
                    keyType = "ENTREZID", # 设定读取的gene ID类型
                    ont = "ALL", # (ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)# 设定q值阈值
    
    ## 7.5 KEGG ----
    # KEGG富集分析
    kk <- enrichKEGG(gene = gene$ENTREZID,
                     keyType = "kegg",
                     organism = KEGG_database,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
    
    ## GO、KEGG结果整合 
    result <- list(enrichGO = kkd, enrichKEGG = kk)
    GO_res <- result$enrichGO
    KEGG_res <- result$enrichKEGG
    
    ## 7.6  res_output ----
    # 导出enrichGO 
    write.xlsx(GO_res@result, file = paste0(dir_WGCNA, "/greenyellowModule_GO_res.xlsx"))
    
    # dotplot
    pdf(file = paste0(dir_WGCNA, "/greenyellowModule_GO.pdf"), width = 6, height = 7)
    p1 <- dotplot(GO_res, showCategory = 5, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free', space = 'free') 
    print(p1)
    dev.off()
    
    # 导出enrichKEGG
    write.xlsx(KEGG_res@result, file = paste0(dir_WGCNA, "/greenyellowModule_KEGG_res.xlsx"))
    
    # dotplot
    pdf(file = paste0(dir_WGCNA, "/greenyellowModule_KEGG.pdf"), width = 7, height = 10)
    p2 <- dotplot(KEGG_res,showCategory = Inf,label_format=50)
    print(p2)
    dev.off()
}

# ggplot2
kegg <- KEGG_res@result[KEGG_res@result$pvalue<0.05,]
p2 <- ggplot(kegg,aes(x=GeneRatio,y=Description))+
    geom_point(aes(size=Count,color= -log10(pvalue)))+
    theme_bw()+labs(y="",x="GeneRatio")+ 
    scale_color_gradient(low="blue",high="red")+
    #scale_size_continuous(range = c(3, 12))+  # 调整气泡的大小范围
    theme(axis.text.y = element_text(angle = 0, hjust = 1))  # 调整Y轴标签角度

ggsave(plot = p2,filename = paste0(dir_WGCNA, "greenyellowModule_kegg.pdf"),height = 7,width = 6)

## 7.7 统计通路的数量 ----
# GO 
table(GO_res@result$p.adjust<0.05)
# KEGG 
table(KEGG_res@result$p.adjust<0.05)
