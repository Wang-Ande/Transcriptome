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
|   |-- Counts_data 【测序数据】
|   |   |-- MOLM13_genecounts.xlsx
|   |   |-- MOLM13_genecounts_sva.xlsx
|   |   |-- MOLM13_group.xlsx
|   |   |-- MV4__11_genecounts.xlsx
|   |   |-- MV4__11_genecounts_sva.xlsx
|   |   |-- MV4__11_group.xlsx
|   |   |-- OCI-AML2_genecounts.xlsx
|   |   |-- OCI-AML2_genecounts_sva.xlsx
|   |   |-- OCI-AML2_group.xlsx
|   |   |-- aml_group.xlsx 【第一批测序数据分组】
|   |   |-- aml_group_1.xlsx 【第二批测序数据分组】
|   |   |-- aml_group_merge.xlsx 【测序数据分组汇总】
|   |   |-- data_merge_adjusted.csv 【批次矫正后counts数据】
|   |   |-- gene_count_matrix.csv 【第一批测序数据】
|   |   |-- gene_count_matrix_1.csv 【第二批测序数据】
|   |   `-- gene_count_matrix_merge.csv 【批次矫正前counts数据】
|   |-- Sample_info 【分组样本信息注释】
|   |   `-- data_group_annotation.csv 
|   `-- Counts_data_by_time 【测序数据按时间维度分组】
|       |-- MOLM13_counts_sva_2W_6W.xlsx
|       |-- MOLM13_counts_sva_WT_2W.xlsx
|       |-- MOLM13_counts_sva_WT_6W.xlsx
|       |-- MOLM13_group_2W_6W.xlsx
|       |-- MOLM13_group_WT_2W.xlsx
|       |-- MOLM13_group_WT_6W.xlsx
|       |-- MV4_11_counts_sva_2W_6W.xlsx
|       |-- MV4_11_counts_sva_WT_2W.xlsx
|       |-- MV4_11_counts_sva_WT_6W.xlsx
|       |-- MV4_11_group_2W_6W.xlsx
|       |-- MV4_11_group_WT_2W.xlsx
|       |-- MV4_11_group_WT_6W.xlsx
|       |-- OCI_AML2_counts_sva_2W_6W.xlsx
|       |-- OCI_AML2_counts_sva_WT_2W.xlsx
|       |-- OCI_AML2_counts_sva_WT_6W.xlsx
|       |-- OCI_AML2_group_2W_6W.xlsx
|       |-- OCI_AML2_group_WT_2W.xlsx
|       `-- OCI_AML2_group_WT_6W.xlsx
|-- 02_code 【代码目录】
|   |-- Combat_seq.R
|   |-- GSEA.R
|   |-- QC_PCA.R
|   |-- QC_boxplot.R
|   |-- QC_heatmap.R
|   |-- Transcriptome_DE.limma_KEGG.GO.R
|   |-- Transcriptome_QC_DE.DEseq2_KEGG.GO.R
|   |-- run_enrichment_analysis.R
|   `-- run_limma_DE.R
|-- 03_result 【结果目录】
|   |-- DE_DEseq2
|   |   |-- MOLM13
|   |   |   |-- 2W_6W_DEseq2.png
|   |   |   |-- DE(NA.omit).csv
|   |   |   |-- DE(NA.omit)_1.csv
|   |   |   |-- DE(NA.omit)_2W_6W.csv
|   |   |   |-- DE(NA.omit)_WT_2W.csv
|   |   |   |-- DE(NA.omit)_WT_6W.csv
|   |   |   |-- DE(sig).csv
|   |   |   |-- DE(sig)_1.csv
|   |   |   |-- DE(sig)_2W_6W.csv
|   |   |   |-- DE(sig)_WT_2W.csv
|   |   |   |-- DE(sig)_WT_6W.csv
|   |   |   |-- Mean-trend-variance.png
|   |   |   |-- WT_2W_DEseq2.png
|   |   |   |-- WT_6W_DEseq2.png
|   |   |   |-- WT_VEN_DEseq2.png
|   |   |   |-- WT_vs_VEN.limma.csv
|   |   |   |-- WT_vs_VEN.limma.pdf
|   |   |   `-- WT_vs_VEN_heatmap.limma.pdf
|   |   |-- MV4_11
|   |   |   |-- 2W_6W_DEseq2.png
|   |   |   |-- DE(NA.omit).csv
|   |   |   |-- DE(NA.omit)_1.csv
|   |   |   |-- DE(NA.omit)_2W_6W.csv
|   |   |   |-- DE(NA.omit)_WT_2W.csv
|   |   |   |-- DE(NA.omit)_WT_6W.csv
|   |   |   |-- DE(sig).csv
|   |   |   |-- DE(sig)_1.csv
|   |   |   |-- DE(sig)_2W_6W.csv
|   |   |   |-- DE(sig)_WT_2W.csv
|   |   |   |-- DE(sig)_WT_6W.csv
|   |   |   |-- DEgene_WT_2W.xlsx
|   |   |   |-- Mean-variance-trend.png
|   |   |   |-- WT_2W_DEseq2.png
|   |   |   |-- WT_6W_DEseq2.png
|   |   |   |-- WT_VEN_DEseq2.png
|   |   |   |-- WT_VEN_limma.csv
|   |   |   |-- WT_vs_VEN.limma.pdf
|   |   |   `-- WT_vs_VEN_heatmap.limma.pdf
|   |   |-- OCI_AML2
|   |   |   |-- 2W_6W_DEseq2.png
|   |   |   |-- DE(NA.omit).csv
|   |   |   |-- DE(NA.omit)_1.csv
|   |   |   |-- DE(NA.omit)_2W_6W.csv
|   |   |   |-- DE(NA.omit)_WT_2W.csv
|   |   |   |-- DE(NA.omit)_WT_6W.csv
|   |   |   |-- DE(sig).csv
|   |   |   |-- DE(sig)_1.csv
|   |   |   |-- DE(sig)_2W_6W.csv
|   |   |   |-- DE(sig)_WT_2W.csv
|   |   |   |-- DE(sig)_WT_6W.csv
|   |   |   |-- Mean-variance-trend.png
|   |   |   |-- WT_2W_DEseq2.png
|   |   |   |-- WT_6W_DEseq2.png
|   |   |   |-- WT_VEN_DEseq2.png
|   |   |   |-- WT_VEN_limma.csv
|   |   |   |-- WT_vs_VEN.limma.pdf
|   |   |   `-- WT_vs_VEN_heatmap.limma.pdf
|   |   `-- README.Rmd
|   |-- DE_limma
|   |   |-- MOLM13_2w_vs_MOLM13_6W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_2W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_6W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_D
|   |   |   |-- DE.csv
|   |   |   |-- Genes.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MV4_11_2w_vs_MV4_11_6w
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_2W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_6W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_D
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- OCI_AML2_2w_vs_OCI_AML2_6w
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- OCI_AML2_WT_vs_OCI_AML2_2W
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   |-- OCI_AML2_WT_vs_OCI_AML2_6w
|   |   |   |-- DE.csv
|   |   |   |-- heatmap.pdf
|   |   |   |-- volc.pdf
|   |   |   `-- voom_plot.pdf
|   |   `-- OCI_AML2_WT_vs_OCI_AML2_D
|   |       |-- DE.csv
|   |       |-- heatmap.pdf
|   |       |-- volc.pdf
|   |       `-- voom_plot.pdf
|   |-- DEseq2_nordata_output
|   |   |-- MOLM13_DEseq2batch.png
|   |   |-- MOLM13_DEseq2batch_rld.csv
|   |   |-- MOLM13_DEseq2batch_vsd.csv
|   |   |-- MOLM13_PCA_DEseq2batch.pdf
|   |   |-- MOLM13_PCA_svabatch.pdf
|   |   |-- MOLM13_heatmap_DEseq2batch.pdf
|   |   |-- MOLM13_heatmap_svabatch.pdf
|   |   |-- MOLM13_svabatch.png
|   |   |-- MOLM13_svabatch_rld.csv
|   |   |-- MOLM13_svabatch_vsd.csv
|   |   |-- MV4_11_DEseq2batch.png
|   |   |-- MV4_11_DEseq2batch_rld.csv
|   |   |-- MV4_11_DEseq2batch_vsd.csv
|   |   |-- MV4_11_PCA_DEseq2batch.pdf
|   |   |-- MV4_11_PCA_svabatch.pdf
|   |   |-- MV4_11_heatmap_DEseq2batch.pdf
|   |   |-- MV4_11_heatmap_svabatch.pdf
|   |   |-- MV4_11_svabatch.png
|   |   |-- MV4_11_svabatch_rld.csv
|   |   |-- MV4_11_svabatch_vsd.csv
|   |   |-- OCI_AML2_DEseq2batch.png
|   |   |-- OCI_AML2_DEseq2batch_rld.csv
|   |   |-- OCI_AML2_DEseq2batch_vsd.csv
|   |   |-- OCI_AML2_PCA_DEseq2batch.pdf
|   |   |-- OCI_AML2_PCA_svabatch.pdf
|   |   |-- OCI_AML2_heatmap_DEseq2batch.pdf
|   |   |-- OCI_AML2_heatmap_svabatch.pdf
|   |   |-- OCI_AML2_svabatch.png
|   |   |-- OCI_AML2_svabatch_rld.csv
|   |   `-- OCI_AML2_svabatch_vsd.csv
|   |-- GO&KEGG_DEseq2
|   |   |-- GO_KEGG_analysis.R
|   |   |-- MOLM13_2w_vs_MOLM13_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_2W_vs_MV4_11_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- OCI_AML2_2W_vs_OCI_AML2_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- OCI_AML2_WT_vs_OCI_AML2_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   `-- OCI_AML2_WT_vs_OCI_AML2_6W
|   |       |-- GO_down.csv
|   |       |-- GO_down.pdf
|   |       |-- GO_up.csv
|   |       |-- GO_up.pdf
|   |       |-- KEGG_down.csv
|   |       |-- KEGG_down.pdf
|   |       |-- KEGG_up.csv
|   |       `-- KEGG_up.pdf
|   |-- GO&KEGG_limma
|   |   |-- MOLM13_2W_vs_MOLM13_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MOLM13_WT_vs_MOLM13_D
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_2w_vs_MV4_11_6w
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_6W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- MV4_11_WT_vs_MV4_11_D
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- OCI_AML2_2w_vs_OCI_AML2_6w
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- OCI_AML2_WT_vs_OCI_AML2_2W
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   |-- OCI_AML2_WT_vs_OCI_AML2_6w
|   |   |   |-- GO_down.csv
|   |   |   |-- GO_down.pdf
|   |   |   |-- GO_up.csv
|   |   |   |-- GO_up.pdf
|   |   |   |-- KEGG_down.csv
|   |   |   |-- KEGG_down.pdf
|   |   |   |-- KEGG_up.csv
|   |   |   `-- KEGG_up.pdf
|   |   `-- OCI_AML2_WT_vs_OCI_AML2_D
|   |       |-- GO_down.csv
|   |       |-- GO_down.pdf
|   |       |-- GO_up.csv
|   |       |-- GO_up.pdf
|   |       |-- KEGG_down.csv
|   |       |-- KEGG_down.pdf
|   |       |-- KEGG_up.csv
|   |       `-- KEGG_up.pdf
|   |-- GSEA
|   |   |-- KEGG_GSEA
|   |   |   |-- Ridgeplot.png
|   |   |   |-- dotplot.png
|   |   |   `-- gseaplot.png
|   |   |-- MOLM13
|   |   |   `-- MsigDBR_GSEA
|   |   |       |-- 2W_6W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       |-- WT_2W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_GSEA_WT_2W.png
|   |   |       |   |-- H_GSEplot_WT_2W.png
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- H_ridgeplot_WT_2W.png
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       |-- WT_6W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       `-- WT_Ven
|   |   |           |-- KEGG_results.csv
|   |   |           `-- gseaplot.pdf
|   |   |-- MV4_11
|   |   |   `-- MsigDBR_GSEA
|   |   |       |-- 2W_6W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       |-- WT_2W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       |-- WT_6W
|   |   |       |   |-- C1_results.csv
|   |   |       |   |-- C2_results.csv
|   |   |       |   |-- C3_results.csv
|   |   |       |   |-- C4_results.csv
|   |   |       |   |-- C5_results.csv
|   |   |       |   |-- C6_results.csv
|   |   |       |   |-- C7_results.csv
|   |   |       |   |-- H_results.csv
|   |   |       |   |-- KEGG_results.csv
|   |   |       |   `-- gseaplot.pdf
|   |   |       `-- WT_Ven
|   |   |           |-- KEGG_results.csv
|   |   |           `-- gseaplot.pdf
|   |   `-- OCI_AML2
|   |       `-- MsigDBR_GSEA
|   |           |-- 2W_6W
|   |           |   |-- C1_results.csv
|   |           |   |-- C2_results.csv
|   |           |   |-- C3_results.csv
|   |           |   |-- C4_results.csv
|   |           |   |-- C5_results.csv
|   |           |   |-- C6_results.csv
|   |           |   |-- C7_results.csv
|   |           |   |-- H_results.csv
|   |           |   |-- KEGG_results.csv
|   |           |   `-- gseaplot.pdf
|   |           |-- WT_2W
|   |           |   |-- C1_results.csv
|   |           |   |-- C2_results.csv
|   |           |   |-- C3_results.csv
|   |           |   |-- C4_results.csv
|   |           |   |-- C5_results.csv
|   |           |   |-- C6_results.csv
|   |           |   |-- C7_results.csv
|   |           |   |-- H_results.csv
|   |           |   |-- KEGG_results.csv
|   |           |   `-- gseaplot.pdf
|   |           |-- WT_6W
|   |           |   |-- C1_results.csv
|   |           |   |-- C2_results.csv
|   |           |   |-- C3_results.csv
|   |           |   |-- C4_results.csv
|   |           |   |-- C5_results.csv
|   |           |   |-- C6_results.csv
|   |           |   |-- C7_results.csv
|   |           |   |-- H_results.csv
|   |           |   |-- KEGG_results.csv
|   |           |   `-- gseaplot.pdf
|   |           `-- WT_Ven
|   |               |-- KEGG_results.csv
|   |               `-- gseaplot.pdf
|   |-- QC
|   |   |-- QC_boxplot.pdf
|   |   |-- QC_boxplot_adjusted.pdf
|   |   |-- QC_heatmap_1.pdf
|   |   |-- QC_heatmap_adjusted_1.pdf
|   |   |-- QC_heatmap_batch.pdf
|   |   |-- QC_heatmap_batch_adjusted.pdf
|   |   |-- QC_heatmap_group.pdf
|   |   |-- QC_pca.pdf
|   |   |-- QC_pca_1_adjusted_1.pdf
|   |   |-- QC_pca_adjusted.pdf
|   |   |-- QC_pca_adjusted_1.pdf
|   |   |-- QC_pca_batch.pdf
|   |   |-- QC_pca_batch_adjusted.pdf
|   |   |-- QC_pca_condition.pdf
|   |   |-- QC_pca_condition_adjusted.pdf
|   |   |-- boxplot.png
|   |   |-- heatmap_adjusted.png
|   |   `-- pca.png
|   `-- Veen
|       |-- DownEVenn(1_MV4_11).csv
|       |-- DownEVenn(1_MV4_11).png
|       |-- DownEVenn(1_MV4_11).svg
|       |-- DownGeneSet(1).xlsx
|       |-- DownGeneSet(1_2W_6W).xlsx
|       |-- DownGeneSet(1_MOLM13_limma).xlsx
|       |-- DownGeneSet(1_MV4_11).xlsx
|       |-- DownGeneSet(1_MV4_11_limma).xlsx
|       |-- DownGeneSet(1_OCI_AML2_limma).xlsx
|       |-- DownGeneSet(1_WT_2W).xlsx
|       |-- DownGeneSet(1_WT_6W).xlsx
|       |-- DownGeneSet(1_WT_VEN_DEseq2).xlsx
|       |-- DownGeneSet(1_WT_VEN_limma).xlsx
|       |-- DownGeneSet(3_sd).xlsx
|       |-- Down_EVenn(2W_6W).csv
|       |-- Down_EVenn(2W_6W).png
|       |-- Down_EVenn(2W_6W).svg
|       |-- Down_EVenn(WT_2W).csv
|       |-- Down_EVenn(WT_2W).png
|       |-- Down_EVenn(WT_2W).svg
|       |-- Down_EVenn(WT_6W).csv
|       |-- Down_EVenn(WT_6W).png
|       |-- Down_EVenn(WT_6W).svg
|       |-- Downgene(1)_EVenn.csv
|       |-- Downgene(1)_EVenn.png
|       |-- Downgene(1)_EVenn.svg
|       |-- Downgene(3_sd)_EVenn.csv
|       |-- Downgene(3_sd)_EVenn.png
|       |-- Downgene(3_sd)_EVenn.svg
|       |-- UP_EVenn(1_WT_2W).csv
|       |-- UpEVenn(1_MV4_11).csv
|       |-- UpEVenn(1_MV4_11).png
|       |-- UpEVenn(1_MV4_11).svg
|       |-- UpGeneSet(1).xlsx
|       |-- UpGeneSet(1_2W_6W).xlsx
|       |-- UpGeneSet(1_MOLM13_limma).xlsx
|       |-- UpGeneSet(1_MV4_11).xlsx
|       |-- UpGeneSet(1_MV4_11_limma).xlsx
|       |-- UpGeneSet(1_OCI_AML2_limma).xlsx
|       |-- UpGeneSet(1_WT_2W).xlsx
|       |-- UpGeneSet(1_WT_6W).xlsx
|       |-- UpGeneSet(1_WT_VEN_DEseq2).xlsx
|       |-- UpGeneSet(1_WT_VEN_limma).xlsx
|       |-- UpGeneSet(3_sd).xlsx
|       |-- Up_EVenn(1_2W_6W).csv
|       |-- Up_EVenn(1_2W_6W).png
|       |-- Up_EVenn(1_2W_6W).svg
|       |-- Up_EVenn(1_WT_2W).png
|       |-- Up_EVenn(1_WT_2W).svg
|       |-- Up_EVenn(1_WT_6W).csv
|       |-- Up_EVenn(1_WT_6W).png
|       |-- Up_EVenn(1_WT_6W).svg
|       |-- Upgene(1)_EVenn.csv
|       |-- Upgene(1)_EVenn.png
|       |-- Upgene(1)_EVenn.svg
|       |-- Upgene(3_sd)_EVenn.csv
|       |-- Upgene(3_sd)_EVenn.png
|       `-- Upgene(3_sd)_EVenn.svg
|-- README.md
`-- Transcriptome.Rproj
