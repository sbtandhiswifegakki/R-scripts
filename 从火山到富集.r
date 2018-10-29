#volcano_plot
rm(list=ls())
library(DESeq2)
getwd()
setwd("/home/shebentao/Desktop")
#at first, you should delete the first row which includes geneid,...
countTable1 <- read.table("2count.txt", sep = "\t", header = T, stringsAsFactors = FALSE, colClasses = c("character","numeric","numeric","numeric","numeric","numeric","numeric"))
rownames(countTable1) <- countTable1$V1
countTable1 <- countTable1[, -1]
colnames(countTable1) <- c("WT_Al1", "WT_Al2", "WT_Al3", "WT_CK1", "WT_CK2", "WT_CK3")
countTable1 <- countTable1[-which(rowSums(countTable1) <= 4), ]
nrow(countTable1)


#DESeq count log2FC&-log10(padj)
colData <- data.frame(row.names = colnames(countTable1), condition = c("Al", "Al", "Al", "CK", "CK", "CK"))
dds <- DESeqDataSetFromMatrix(countData = countTable1, colData = colData, design = ~ condition)
dds
dds$condition <- relevel(dds$condition, "CK")
dds$condition
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$padj), ]
write.csv(resOrdered, "DESeq2_result_all_2.csv")
deg <- subset(resOrdered, padj <= 0.01 & abs(log2FoldChange) >= 1)
summary(deg)
write.csv(deg, "DESeq2_result_significant_2.csv")


data1 <- read.csv("DESeq2_result_all_2.csv", stringsAsFactors = FALSE)
data1 <- na.omit(data1)
data1$change <- as.factor(ifelse(data1$padj<0.01 & abs(data1$log2FoldChange) >1, 
                                 ifelse(data1$log2FoldChange >1, 'UP:804', 'DOWN:391'), 'NOT'))
table(data1$change)
391 + 804
library(ggplot2)
p <- ggplot(data = data1, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha = 0.8, size = 1) + theme_classic(base_size = 15) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  scale_color_manual(name = "DEG(1195)", values = c("red", "green", "black"), limits = c("UP:804", "DOWN:391", "NOT")) 
p
rownames(data1) <- data1$X
data1$sign <- ifelse(-log10(data1$padj) > 10 & abs(data1$log2FoldChange) > 3, rownames(data1), NA)
library(ggrepel)
b <- p + geom_text_repel(aes(label=data1$sign)) + labs(title = "78_Al VS 78_CK")
b
library(plotly)
b <- ggplotly(b)
b

# heatmap
rm(list = ls())
setwd("/home/shebentao/Desktop")
getwd()

deseq_results_significant <- read.csv("DESeq2_result_significant_2.csv", row.names = "X")
significant_genes <- rownames(deseq_results_significant)
fpkm_gtf <- read.table("2fpkm.gtf", stringsAsFactors = FALSE)
rownames(fpkm_gtf) <- fpkm_gtf$V1
fpkm_significant_genes <- fpkm_gtf[significant_genes, 8:13]
colnames(fpkm_significant_genes) <- c("78_Al1", "78_Al2", "78_Al3", "78_CK1", "78_CK2", "78_CK3")
fpkm_significant_genes <- na.omit(fpkm_significant_genes)
nrow(fpkm_significant_genes)
library(pheatmap)
pheatmap(log2(fpkm_significant_genes + 1), show_rownames = FALSE, main = '78_Al VS 78_CK')
pheatmap(log2(fpkm_significant_genes[1:50,] + 1), show_rownames = TRUE, main = '78_Al VS 78_CK top50')

#Enrichment
setwd("/home/shebentao/Desktop")
library(clusterProfiler)
library(org.At.tair.db)
table <- read.csv("DESeq2_result_2_significant.csv", row.names = "X")
test <- bitr((table), fromType = "TAIR", toType = "ENTREZID", OrgDb = "org.At.tair.db")

geneID_file <- read.table("2_geneid.txt", quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_id <- geneID_file$To
go_results <- enrichGO(gene_id, 'org.At.tair.db', ont = "BP")
barplot(go_results, showCategory = 20, title = "GO_Enrichment_BP")
write.csv(as.data.frame(go_results), "deg_1_goEnrichment_BP(1947).csv")

go_results_MF <- enrichGO(gene_id, 'org.At.tair.db', ont = "MF")
barplot(go_results_MF, showCategory = 20, title = "GO_Enrichment_MF")
write.csv(as.data.frame(go_results_MF), "2_GO_enriched_MF_results.csv")

go_results_CC <- enrichGO(gene_id, 'org.At.tair.db', ont = "CC")
barplot(go_results_CC, showCategory = 20, title = "GO_Enrichment_CC")
write.csv(as.data.frame(go_results_CC), "2_GO_enriched_CC_results.csv")

kegg_results <- enrichKEGG(gene_id, organism = "ath", keyType = "ncbi-geneid")
dotplot(kegg_results, x = "GeneRatio", color = "p.adjust", split = NULL, title = "2_KEGG_dotplot")
write.csv(as.data.frame(kegg_results), "2_KEGG_enriched_results.csv")

a <- read.table("78_CK VS 78_Al_GO_special-up-enriched_BP_results.csv", sep = "," ,header = TRUE, stringsAsFactor = FALSE)
B <- a[1:2,]
c <- as.data.frame(B)
write.csv(c, "78 response to fungus&bacterium.csv",row.names = FALSE)
d <- c$geneID
d <- as.data.frame(d)
write.csv(d, "123.CSV", sep = ",",row.names = FALSE)
e <- as.vector(d[1,])
f <- as.vector(d[2,])
g <- union(e,f)
