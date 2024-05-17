#This script explore the difference of immune cell fractions among subclusters and clinical indicators in HPC

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

group$condition <- ifelse(group$cluster=="0"|group$cluster=="6","control","case")
group$condition <- as.factor(group$condition)

#delete duplicate
exprSet <- exprSet[,colnames(exprSet)%in%name]
exprSet <- exprSet[,rownames(group)]
all(rownames(group) %in% colnames(exprSet))
all(rownames(group) == colnames(exprSet))

condition <- group$condition

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


########################dataset
dataset <- read.table("data_no_multicol.txt",header=T,sep = "\t")

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

dataset$condition <- condition
dataset$condition <- as.factor(dataset$condition)


all(rownames(dataset) %in% rownames(group))
all(rownames(dataset) == rownames(group))



############Dunn's test

library(FSA)
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$CD14_positive_monocyte)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$CD16_positive_monocyte)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$CD4_positive_alpha_beta_T_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$CD56bright_natural_killer_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$CD8_positive_alpha_beta_T_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$MAST_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$basophil)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$eosinophil)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$gamma_delta_T_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$hematopoietic_progenitor)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$macrophage_m0)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")	
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$macrophage_m1)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$macrophage_m2)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$memory_B_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$myeloid_dendritic_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$naive_B_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")
		 
test_data <- data.frame(group = as.factor(group$cluster),value =dataset$neutrophil)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$plasma_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")	

test_data <- data.frame(group = as.factor(group$cluster),value =dataset$plasmacytoid_dendritic_cell)
dunnTest(value ~ group,
         data=test_data,
         method="bonferroni")		 




######################plot
C4_list <- rownames(group[group$cluster==4,])
data_C4 <- dataset[rownames(dataset) %in% C4_list,]

C5_list <- rownames(group[group$cluster==5,])
data_C5 <- dataset[rownames(dataset) %in% C5_list,]

mean_C4 <- mean(data_C4$hematopoietic_progenitor)

mean_C5 <- mean(data_C5$hematopoietic_progenitor)

		 
table(dataset$COMP_1)
mean_A1 <- mean(dataset[dataset$COMP_1<0.5,]$macrophage_m0)  #0.05374246
mean_A2 <- mean(dataset[dataset$COMP_1==0.5,]$macrophage_m0)  #0.06305216
mean_A3 <- mean(dataset[dataset$COMP_1==1,]$macrophage_m0)      #0.1694528


table(dataset$COMPULSIONS)
mean_B1 <- mean(dataset[dataset$COMPULSIONS==0,]$macrophage_m0)  #0.05419509
mean_B2 <- mean(dataset[dataset$COMPULSIONS==1,]$macrophage_m0)  #0.2083797

		 
table(dataset$PHQ_4)
mean_C1 <- mean(dataset[dataset$PHQ_4==0,]$memory_B_cell)  #0.007474954
mean_C2 <- mean(dataset[dataset$PHQ_4==0.3333,]$memory_B_cell)  #0.007208405
mean_C3 <- mean(dataset[dataset$PHQ_4==0.6667,]$memory_B_cell) #0.03281304
mean_C4 <- mean(dataset[dataset$PHQ_4==1,]$memory_B_cell)  #0.006839979


table(dataset$Phq_Cat)
mean_C1 <- mean(dataset[dataset$Phq_Cat==0,]$plasma_cell)  #0.007650758
mean_C2 <- mean(dataset[dataset$Phq_Cat==0.25,]$plasma_cell)  #0.004686133
mean_C3 <- mean(dataset[dataset$Phq_Cat==0.5,]$plasma_cell) #0.005859679
mean_C4 <- mean(dataset[dataset$Phq_Cat==0.75,]$plasma_cell)  #0.002446383
mean_C5 <- mean(dataset[dataset$Phq_Cat==1,]$plasma_cell)  #0.08981529