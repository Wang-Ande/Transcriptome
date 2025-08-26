exprSet <- read.csv("./01_data/gene_tpm_transcript_matrix.csv", row.names = 1)
exprSet <- as.data.frame(exprSet)
data_group <- read_excel("./01_data/group_info.xlsx")
colnames(exprSet) <- gsub("w","W",colnames(exprSet))
colnames(exprSet)
data_group$id <- colnames(exprSet)
data_group <- data_group[,-3]
data_group[grep("WT",data_group$id),2] <- "WT"
table(data_group$group)
library(stringr)
data_group$cell <- gsub("_WT_1|_WT_2|_WT_3|_2W_1|_2W_2|_2W_3|_4W_2|_6W_1|_6W_2|_6W_3|","",data_group$id)
write.xlsx(data_group , file = "./01_data/group_info.xlsx")
write.csv(exprSet, file = "./01_data/gene_tpm_transcript_matrix.csv")


colnames(exprSet) <- gsub("w", "W", colnames(exprSet))
colnames(exprSet) <- gsub("MV", "MV4", colnames(exprSet))
colnames(exprSet) <- gsub("\\.", "_", colnames(exprSet))
