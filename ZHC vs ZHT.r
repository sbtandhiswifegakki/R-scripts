getwd()
setwd("/home/shebentao/Desktop")
getwd()
ls()
library(DESeq2)
countTable1 <- read.table("count.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
tail(countTable1)
tail(countTable1)
rownames(countTable1) <- countTable1$V1
tail(countTable1)
countTable1 <- countTable1[, -1]
tail(countTable1)
colnames(countTable1) <- c("ZHC", "ZHT")
tail(countTable1)
head(countTable1)
countTable1 <- countTable1[-1, ]
head(countTable1)
tail(countTable1)
str(countTable1)

#filter
countTable1$ZHC = as.numeric(countTable1$ZHC)
countTable1$ZHT = as.numeric(countTable1$ZHT)
countTable1$WT_Al3 = as.numeric(countTable1$WT_Al3)
countTable1$WT_CK1 = as.numeric(countTable1$WT_CK1)
countTable1$WT_CK2 = as.numeric(countTable1$WT_CK2)
countTable1$WT_CK3 = as.numeric(countTable1$WT_CK3)
str(countTable1)

countTable1 <- countTable1[-which(rowSums(countTable1) < 4), ]
nrow(countTable1)
countTable1

#DESeq
colData <- data.frame(row.names = colnames(countTable1), condition = c("CK", "Al"))
dds <- DESeqDataSetFromMatrix(countData = countTable1, colData = colData, design = ~ condition)
dds
dds$condition <- relevel(dds$condition, "CK")
dds$condition
dds <- DESeq(dds)
res <- results(dds)
summary(res)
resOrdered <- res[order(res$padj), ]
str(resOrdered)
write.csv(resOrdered, "DESeq2_results_all.csv")
deg <- subset(resOrdered, padj <= 0.01 & abs(log2FoldChange) >= 2)
summary(deg)
write.csv(deg, "DESeq2_results_significant.csv")


#volcano plot
library(dplyr)
install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
volcano_data <- read.csv("DESeq2_results_all.csv", row.names = "X")
volcano_data <- na.omit(volcano_data)
threshold <- as.factor(abs(volcano_data$log2FoldChange) > 2 & volcano_data$padj <= 0.01)
p = ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(shape = threshold, color = threshold)) +
  scale_y_continuous(limits = c(0,200), expand = c(0,0)) + scale_shape_discrete(labels = c("not significant", "significant")) +
  scale_color_discrete(labels = c("not significant", "significant")) + labs(title = "WT_Al VS WT_CK(A.th)") 
p + geom_text_repel(aes(log2FoldChange, -log10(padj), label = rownames(volcano_data)))                                                                                                                        
