#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#The aim is to normalize read count and identify biology confounders



rm(list=ls())
setwd("/panfs/compbio/users/zzhu244/volumes")
matrix1 <- read.table("matrix_counts1.txt",header = T)
matrix2 <- read.table("matrix_counts2.txt",header = T)
matrix3 <- read.table("matrix_counts3.txt",header = T)
matrix4 <- read.table("matrix_counts4.txt",header = T)

matrix <- merge(matrix1,matrix2,by="gene_id")
matrix <- merge(matrix,matrix3,by="gene_id")
matrix <- merge(matrix,matrix4,by="gene_id")

#conduct transform
mRNAmatrix2<- matrix[,-1]
mRNAmatrix3<-as.matrix(mRNAmatrix2)
rownames(mRNAmatrix3)<-matrix$gene_id
exprSet <- mRNAmatrix3

#read clinical data
load("cluster_index.RData")
spname <- colnames(matrix) 
spname <- spname[-1]
spname <- as.data.frame(spname)
group <- merge(spname,data_index,by.x="spname",by.y="id")
name <- group$spname
rownames(group) <- group$spname
group$spname <- as.factor(group$spname)
group$cluster <- as.factor(group$cluster)
group$site <- as.factor(group$site)


#delete duplicate
exprSet <- exprSet[,colnames(exprSet)%in%name]
exprSet <- exprSet[,rownames(group)]
all(rownames(group) %in% colnames(exprSet))
all(rownames(group) == colnames(exprSet))





#Gene-level counts normalize

library("limma")
library("edgeR")
library("variancePartition")

# identify genes that pass expression cutoff
#remove 15815 genes
isexpr <- rowSums(cpm(exprSet) > 1) >= 0.5 * ncol(exprSet)

# create data structure with only expressed genes
gExpr <- DGEList(counts = exprSet[isexpr, ])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
design <- model.matrix(~1, data=data.frame(group = rep(1,ncol(exprSet))))

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design)
vobjGenes_matrix <- vobjGenes$E
vobjGenes_matrix <- t(vobjGenes_matrix)


#conduct PCA
library(FactoMineR)
library(factoextra)

result <- PCA(vobjGenes_matrix,scale.unit = T,graph = F)

eig.val <- get_eigenvalue(result)
head(eig.val)

pdf("elbow plot.pdf")

fviz_screeplot(result)

dev.off()

head(scale(result$ind$coord))
dimension <- scale(result$ind$coord)

dataset <- read.table("data_no_multicol.txt",header=T,sep = "\t")

all_dataset <- cbind(dataset,dimension)


#######regression on num

pvalue <- c()
for (i in 1:27)
{
reg <- c()
reg <- lm(all_dataset$Dim.1~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[1:27]

pvalue <- c()
for (i in 1:27)
{
reg <- c()
reg <- lm(all_dataset$Dim.2~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")



pvalue <- c()
for (i in 1:27)
{
reg <- c()
reg <- lm(all_dataset$Dim.3~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")


###############regression on category
pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$Dim.1~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$Dim.2~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")



pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$Dim.3~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")



#####PCA plot
pdf("Ppreg4 plot.pdf")
fviz_cluster(list(data=dimension[,1:2], cluster = as.factor(all_dataset$Ppreg4)),geom = c("point"))
dev.off()


pdf("Agegroup plot.pdf")
fviz_cluster(list(data=dimension[,1:2], cluster = as.factor(all_dataset$Agegroup)),geom = c("point"))
dev.off()


pdf("gender plot.pdf")
fviz_cluster(list(data=dimension[,1:2], cluster = as.factor(all_dataset$gender)),geom = c("point"))
dev.off()




############Variance partitioning analysis
all(rownames(dataset) %in% colnames(exprSet))
all(rownames(dataset) == colnames(exprSet))

dataset$Ppreg4 <- as.factor(dataset$Ppreg4)
dataset$Agegroup <- as.factor(dataset$Agegroup)
dataset$gender <- as.factor(dataset$gender)

scale_norm <- function(data)
{
	normalized_data <- (data - min(data)) / (max(data) - min(data))
	return(normalized_data)
}

dataset$CD14_positive_monocyte <- scale_norm(dataset$CD14_positive_monocyte)
dataset$CD16_positive_monocyte <- scale_norm(dataset$CD16_positive_monocyte)
dataset$CD4_positive_alpha_beta_T_cell <- scale_norm(dataset$CD4_positive_alpha_beta_T_cell)
dataset$CD56bright_natural_killer_cell <- scale_norm(dataset$CD56bright_natural_killer_cell)
dataset$eosinophil <- scale_norm(dataset$eosinophil)
dataset$gamma_delta_T_cell <- scale_norm(dataset$gamma_delta_T_cell)
dataset$macrophage_m2 <- scale_norm(dataset$macrophage_m2)
dataset$memory_B_cell <- scale_norm(dataset$memory_B_cell)
dataset$naive_B_cell <- scale_norm(dataset$naive_B_cell)
dataset$neutrophil <- scale_norm(dataset$neutrophil)
dataset$plasma_cell <- scale_norm(dataset$plasma_cell)
dataset$plasmacytoid_dendritic_cell <- scale_norm(dataset$plasmacytoid_dendritic_cell)
dataset$CD8_positive_alpha_beta_T_cell <- scale_norm(dataset$CD8_positive_alpha_beta_T_cell)
dataset$basophil <- scale_norm(dataset$basophil)
dataset$macrophage_m0 <- scale_norm(dataset$macrophage_m0)
dataset$macrophage_m1 <- scale_norm(dataset$macrophage_m1)
dataset$myeloid_dendritic_cell <- scale_norm(dataset$myeloid_dendritic_cell)







# Define formula
form <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD16_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ eosinophil+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ neutrophil+ plasma_cell+ plasmacytoid_dendritic_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell                                                                                                                  
# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
#form1 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric
#varPart1 <- fitExtractVarPartModel(vobjGenes, form1, dataset)
varPart <- fitExtractVarPartModel(vobjGenes, form, dataset)


vp <- sortCols(varPart)

pdf("bar plot_deseq2.pdf")

plotPercentBars(vp[1:10, ])

dev.off()

pdf("violin plot_deseq2.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()



##############remove neutrophil

# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil, dataset,
  fxn = residuals
)


# convert list to matrix
residMatrix <- do.call(rbind, residList)

form1 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD16_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ eosinophil+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ plasma_cell+ plasmacytoid_dendritic_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form1, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r1.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()


##############remove eosinophil
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form2 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD16_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ plasma_cell+ plasmacytoid_dendritic_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form2, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r2.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()




##############remove CD16_positive_monocyte
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil+CD16_positive_monocyte, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form3 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ plasma_cell+ plasmacytoid_dendritic_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form3, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r3.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()



##############remove plasmacytoid_dendritic_cell
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil+CD16_positive_monocyte+ plasmacytoid_dendritic_cell, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form4 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ plasma_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form4, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r4.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()




##############remove macrophage_m2
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil+CD16_positive_monocyte+ plasmacytoid_dendritic_cell+ macrophage_m2, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form5 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ gamma_delta_T_cell+ memory_B_cell+ naive_B_cell+ plasma_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form5, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r5.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()




##############remove gamma_delta_T_cell
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil+CD16_positive_monocyte+ plasmacytoid_dendritic_cell+ macrophage_m2+ gamma_delta_T_cell, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form6 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell + memory_B_cell+ naive_B_cell+ plasma_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form6, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r6.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()






##############remove CD14_positive_monocyte
# extract residuals directly without storing intermediate results
residList <- fitVarPartModel(vobjGenes, ~ neutrophil+eosinophil+CD16_positive_monocyte+ plasmacytoid_dendritic_cell+ macrophage_m2+ gamma_delta_T_cell+ CD14_positive_monocyte, dataset,
  fxn = residuals
)

# convert list to matrix
residMatrix <- do.call(rbind, residList)

form7 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric + CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell + memory_B_cell+ naive_B_cell+ plasma_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell 

varPartResid <- fitExtractVarPartModel(residMatrix, form7, dataset)


vp <- sortCols(varPartResid)


pdf("violin plot_r7.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()
































