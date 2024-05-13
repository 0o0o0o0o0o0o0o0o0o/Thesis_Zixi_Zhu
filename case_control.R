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

# The variable to be tested must be a fixed effect
form <- ~ condition + neutrophil+eosinophil+CD16_positive_monocyte+ plasmacytoid_dendritic_cell+ macrophage_m2+ gamma_delta_T_cell+ CD14_positive_monocyte

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(gExpr, form, dataset)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, dataset)
fitmm <- eBayes(fitmm)

head(fitmm$design, 3)

top30 <- topTable(fitmm, coef = "conditioncontrol", number = 30,sort.by = "P")
top60 <- topTable(fitmm, coef = "conditioncontrol", number = 60,sort.by = "P")
top100 <- topTable(fitmm, coef = "conditioncontrol", number = 100,sort.by = "P")
top150 <- topTable(fitmm, coef = "conditioncontrol", number = 150,sort.by = "P")
top300 <- topTable(fitmm, coef = "conditioncontrol", number = 300,sort.by = "P")
top500 <- topTable(fitmm, coef = "conditioncontrol", number = 500,sort.by = "P")





######################GSEA
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

library(org.Hs.eg.db)

#########top30
df <- top30
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.1,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")


#########top60
df <- top60
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.05,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")


#########top100
df <- top100
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.05,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
			 
			 
			 
			 

#########top150
df <- top150
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.05,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")



#########top300
df <- top300
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.05,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
			 
			 


#########top500
df <- top500
# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

#gsea
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
			 pvalueCutoff = 0.1,
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

























