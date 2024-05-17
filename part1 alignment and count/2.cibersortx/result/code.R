rm(list=ls())
setwd("D:\\Thesis\\data\\result_tpm")

data <- read.table("CIBERSORTx_Adjusted.txt",header = T)
rownames(data) <- data$Mixture
data <- data[,-1]
library(pheatmap)
color_gradient <- colorRampPalette(c("white", "red"))(100)
colnames_string <- paste(colnames(data[1:50, 1:20]), collapse = "")
pdf("heatmap.pdf", width = 10, height = 100) # 根据需要调整PDF的宽度和高度

par(mar = c(5, 5, 5, 5))
pheatmap(data[,1:20], 
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         display_numbers=TRUE, 
         color=color_gradient,
         annotation_legend = FALSE,
         show_colnames = TRUE,
         )

dev.off()





