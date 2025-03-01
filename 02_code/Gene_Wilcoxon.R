# 删除工作空间中所有的对象
rm(list = ls()) 

# 加载包 ----
library(rstatix)
library(tidyverse)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggpubr) # 添加显著性

# 设置工作路径 ----
setwd('.') 
WT_2W
# 创建输出目录 ----
if(!dir.exists('./03_result/Wilcoxon/FPKM/OCI_AML2')){
  dir.create('./03_result/Wilcoxon/FPKM/MOLM13/WT_2W')
} 

# 设置输出路径
dir <- "./03_result/Wilcoxon/FPKM/MV4_11/WT_6W/"

# 下方代码为增加运行速度，添加if语句，
#需要注意更改表达量和分组数据导入即可 !!!!!!!!!!!!!

# 导入测序数据 ----
if(T){ 
  TrainRawData <- read.xlsx("./01_data/fpkm_data_by_time/MV4_11_fpkm_sva_WT_6W.xlsx") 

#提取symbol名 ----
y <- TrainRawData$X1
gene <- unlist(lapply(y,function(y) strsplit(as.character(y),"\\|")[[1]][2]))
id <- data.frame(gene)
TrainRawData <- cbind(id,TrainRawData)

# 检查重复基因
#duplicate_gene <- TrainRawData[duplicated(TrainRawData$gene) ,] 

# 根据gene去除重复行，保留第一条记录
TrainRawData <- TrainRawData %>%
  distinct(gene, .keep_all = TRUE)  

# 设置symbol为行名
rownames(TrainRawData) <- TrainRawData$gene
TrainRawData <- TrainRawData[,-c(1,2)]

# 导入分组信息 ----
TrainGroup <- read.xlsx("./01_data/fpkm_data_by_time/MV4_11_group_WT_6W.xlsx") 
colnames(TrainGroup) <- c('sample', 'group')

# 导入筛选基因 ----
HubGene <- data.frame(symbol = c("BAX","BCL2","MCL1","TP53","BAK1"))

# 从raw数据中筛选出hubgene的表达量
TrainData <- TrainRawData[HubGene$symbol, ] %>% t() %>% as.data.frame()
TrainData <- merge(TrainGroup, TrainData, by.x = "sample", by.y = 'row.names')

# 将数据转换为长格式 
TrainData <- TrainData %>% 
  pivot_longer(
    cols = -c("sample", "group"),
    names_to = "symbol",
    values_to = "Expression"
  )
TrainData$Expression <- as.numeric(TrainData$Expression)
# log2转换(根据情况)
#TrainData$Expression <- log2(TrainData$Expression)

# 基因Wilcoxon"配对"秩和检验 ----
WilcoxonResults <- TrainData%>%
  group_by(symbol)%>%                    # 按symbol对每个数据分组
  wilcox_test(Expression ~ group,
              paired = T)%>%          # 比较表达量在不同组别间的差异
  adjust_pvalue(method = 'fdr')          # 使用FDR方法来调整p值

# 箱线图可视化 ----
pdf(paste0(dir , "boxplot_Wilcoxon.pdf"), width = 8, 
    height = 6)
ggplot(TrainData, aes(x = symbol, y = Expression, fill = group)) +
  stat_boxplot(geom = "errorbar",
               width = 0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(aes(x = symbol, y = Expression, fill = group),
               width = 0.2,
               position = position_dodge(0.9), 
               outlier.shape = NA, 
               outlier.colour = NA)+ 
  scale_fill_manual(values = c('#355783', "gold"), name = "Group")+
  labs(title = "", x = "", y = "Expression", size = 20) +
  stat_compare_means(data = TrainData,
                     mapping = aes(group = group),#按照group变量来进行分组比较
                     label = "p.signif",         #显著性标记（ ***, **, *, ns）
                     method = 'wilcox.test',
                     paired = T) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust=1, colour = "black", face = "bold", size = 10), 
        axis.text.y = element_text(hjust = 0.5, colour ="black", face="bold", size=12), 
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

}


