#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#The aim is to conduct cibersortx to get immune cell fractions


rm(list=ls())
setwd("/panfs/compbio/users/zzhu244/volumes")
matrix1 <- read.table("matrix_counts1.txt",header = T)
matrix2 <- read.table("matrix_counts2.txt",header = T)
matrix3 <- read.table("matrix_counts3.txt",header = T)
matrix4 <- read.table("matrix_counts4.txt",header = T)

matrix <- merge(matrix1,matrix2,by="gene_id")
matrix <- merge(matrix,matrix3,by="gene_id")
matrix <- merge(matrix,matrix4,by="gene_id")



#Then apply voom in limma to normalization
suppressMessages(library(limma))
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("limma")

mRNAmatrix2<- matrix[,-1]
mRNAmatrix3<-as.matrix(mRNAmatrix2)
rownames(mRNAmatrix3)<-matrix$gene_id
exprSet <- mRNAmatrix3


library("readxl")
set2 <- read_excel("DepressionGenesNetworks_case_control_phenotype_distribution.xlsx", sheet = "Data", col_names = T)
spname <- colnames(matrix) 
spname <- spname[-1]
set2 <- set2[,c("Cell_id","case_control_status")]
spname <- as.data.frame(spname)
group <- merge(spname,set2,by.x="spname",by.y="Cell_id")
name <- group$spname
group <- group$case_control_status
group <- as.factor(group)

exprSet <- exprSet[,colnames(exprSet)%in%name]

#run script genelength.R to get gene length
load("genelength.RData")

geneid_efflen_gtf1 <- geneid_efflen_gtf1[geneid_efflen_gtf1$geneid%in%row.names(exprSet),]
#remove LD0053
exprSet <- exprSet[,colnames(exprSet)!="LD0053"]

#select gene that whole expression >0
expr_df <- exprSet[rowSums(exprSet)>0,]
geneid_efflen_gtf1 <- geneid_efflen_gtf1[geneid_efflen_gtf1$geneid%in%row.names(expr_df),]


kb <- geneid_efflen_gtf1$efflen/1000
save.image("beforetpm.RData")

#run tpm_cal.R to normalize to TPM

load("aftertpm.RData")
write.table(tpm,file = "matrix_counts_TPM.txt",quote = F,row.names = T,sep="\t")


#this code should be run in the docker version of cibersortx. This is not R code!
docker run -v D:/Thesis/data:/src/data -v D:/Thesis/data/result_tpm:/src/outdir cibersortx/fractions --username zixi.zhu@emory.edu --token 1c649ba22f955aeeed670ec1fa7539c3 --mixture matrix_counts_TPM.txt --sigmatrix ISM.txt --rmbatchBmode TRUE --QN TRUE --perm 100









