#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#The aim is to try to use DESeq2 method to normalize read count(an alternative method) and explore the association between immune cell fractions and clinical variables



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


#deseq2
library("DESeq2")
group$design <- ifelse(group$cluster==0|group$cluster==6,0,1)
group$design <- as.factor(group$design)

# create DESeq2 object from gene-level counts and metadata
dds <- DESeqDataSetFromMatrix(
  countData = exprSet,
  colData = group,
  design = ~1
)

# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds) > 1) >= 0.5 * ncol(dds)

# compute log2 Fragments Per Million
# Alternatively, fpkm(), vst() or rlog() could be used
quantLog <- log2(fpm(dds)[isexpr, ] + 1)



#conduct PCA
library(FactoMineR)
library(factoextra)

result <- PCA(t(quantLog),scale.unit = T,graph = F)

eig.val <- get_eigenvalue(result)
head(eig.val)

pdf("elbow plot_deseq2.pdf")

fviz_screeplot(result)

dev.off()

head(scale(result$ind$coord))
dimension <- scale(result$ind$coord)

dataset <- read.table("data_no_multicol.txt",header=T,sep = "\t")

all_dataset <- cbind(dataset,dimension)


pdf("gender plot_deseq2.pdf")
fviz_cluster(list(data=dimension[,1:2], cluster = as.factor(all_dataset$gender)),geom = c("point"))
dev.off()


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
#colnames(all_dataset)[which(adjusted_p_values<0.05)]

pvalue <- c()
for (i in 1:27)
{
reg <- c()
reg <- lm(all_dataset$Dim.2~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]



###############regression on category
pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$Dim.1~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]




pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$Dim.2~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]





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






library("variancePartition")
# Define formula
form <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric+ CD14_positive_monocyte+ CD16_positive_monocyte+ CD4_positive_alpha_beta_T_cell+ CD56bright_natural_killer_cell+ eosinophil+ gamma_delta_T_cell+ macrophage_m2+ memory_B_cell+ naive_B_cell+ neutrophil+ plasma_cell+ plasmacytoid_dendritic_cell+ CD8_positive_alpha_beta_T_cell+ basophil+ macrophage_m0+ macrophage_m1+ myeloid_dendritic_cell                                                                                                                  
# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
#form1 <- ~ (1 | Ppreg4) + (1 | Agegroup)+ (1 | gender)+ Bmi_Current+ Timebloodnumeric
#varPart1 <- fitExtractVarPartModel(vobjGenes, form1, dataset)
varPart <- fitExtractVarPartModel(quantLog, form, dataset)


vp <- sortCols(varPart)

pdf("bar plot2.pdf")

plotPercentBars(vp[1:10, ])

dev.off()

pdf("violin plot2.pdf")

plotVarPart(vp,label.angle = 90)

dev.off()







############explore the association between immune cell fractions and clinical variables


#CD14_positive_monocyte

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$CD14_positive_monocyte~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$CD14_positive_monocyte~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]



#CD16_positive_monocyte

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$CD16_positive_monocyte~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$CD16_positive_monocyte~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]




#CD4_positive_alpha_beta_T_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$CD4_positive_alpha_beta_T_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$CD4_positive_alpha_beta_T_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]




#CD56bright_natural_killer_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$CD56bright_natural_killer_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$CD56bright_natural_killer_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]







#CD8_positive_alpha_beta_T_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$CD8_positive_alpha_beta_T_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$CD8_positive_alpha_beta_T_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]





#MAST_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$MAST_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$MAST_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#basophil

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$basophil~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$basophil~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#eosinophil

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$eosinophil~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$eosinophil~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#gamma_delta_T_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$gamma_delta_T_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$gamma_delta_T_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]







#hematopoietic_progenitor

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$hematopoietic_progenitor~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$hematopoietic_progenitor~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]







#macrophage_m0

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m0~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m0~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]








#macrophage_m1

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m1~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m1~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]








#macrophage_m2

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m2~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$macrophage_m2~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#memory_B_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$memory_B_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$memory_B_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]








#naive_B_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$naive_B_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$naive_B_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#neutrophil

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$neutrophil~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$neutrophil~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]






#plasma_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$plasma_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$plasma_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]





#plasmacytoid_dendritic_cell

pvalue <- c()
for (i in 1:9)
{
reg <- c()
reg <- lm(all_dataset$plasmacytoid_dendritic_cell~all_dataset[[i]])
pvalue <- c(pvalue,summary(reg)$coefficients[2,4])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")
#which(adjusted_p_values<0.05)
#colnames(all_dataset)[which(adjusted_p_values<0.05)]


pvalue <- c()
for (i in 28:312)
{
reg <- c()
reg <- lm(all_dataset$plasmacytoid_dendritic_cell~factor(all_dataset[[i]]))
pvalue <- c(pvalue,anova(reg)[1,5])
}

adjusted_p_values <- p.adjust(pvalue, method = "BH")

#which(adjusted_p_values<0.05)+27
#colnames(all_dataset)[which(adjusted_p_values<0.05)+27]



















































