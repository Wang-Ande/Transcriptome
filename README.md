---
output:
  html_document: default
  pdf_document: default
---
# Transcriptome
 此项目为转录组下游数据处理

## 功能实现

此脚本通过读入RNA-seq counts数据，QC处理（boxplot，PCA plot和heatmap），
差异分析(Limma/DEseq2)、GO/KEGG通路富集分析，以及GSEA富集分析。

## 目录架构
|-- 01_data 【数据目录】
|-- 02_code 【代码目录】
|-- 03_result 【结果目录】
