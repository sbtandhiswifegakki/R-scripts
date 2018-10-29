rm(list = ls())
setwd("/home/shebentao/Desktop")
getwd()

deseq_results_significant <- read.csv("DESeq2_result_1_significant.csv", row.names = "X")
significant_genes <- rownames(deseq_results_significant)
fpkm_gtf <- read.table("1fpkm.gtf", stringsAsFactors = FALSE)
rownames(fpkm_gtf) <- fpkm_gtf$V1
str(fpkm_gtf)
fpkm_gtf <- fpkm_gtf[-which(fpkm_gtf$V1 == "-" | fpkm_gtf$V1 == "."), ]
str(fpkm_gtf)
fpkm_gtf
fpkm_significant_genes <- fpkm_gtf[significant_genes, 8:13]

colnames(fpkm_significant_genes) <- c("WT_Al1", "WT_Al2", "WT_Al3", "WT_CK1", "WT_CK2", "WT_CK3")

fpkm_significant_genes <- na.omit(fpkm_significant_genes)
nrow(fpkm_significant_genes)

library(pheatmap)
pheatmap(log2(fpkm_significant_genes + 1), show_rownames = FALSE)
pheatmap(log2(t(fpkm_significant_genes[1:30, ] + 1)), show_colnames = TRUE)
pheatmap(log2(fpkm_significant_genes[1:30, ] + 1))



#have problems