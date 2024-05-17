#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#The aim is to get count matrix in each gap folder




rm(list=ls())
###################GAP1
#bulid date.frame about sample name and file name
sample_name<-read.table("//home/zzhu244/thesis/small_sample_test/GAP1.txt")
colnames(sample_name)[1]<-"sample_name"
sample_name$file_name <- paste(sample_name$sample_name, "ReadsPerGene.out.tab", sep = "")

#read file
setwd("/panfs/compbio/users/zzhu244/volumes/GAP1/qc_output/outcome")
require(data.table)
#expr_df <- fread(sample_name$file_name[30], skip = 4, select = 1:2)
#names(expr_df) <- c("gene_id",sample_name$sample_name[30])
expr_df <- fread(sample_name$file_name[1], skip = 4, select = 1:2)
names(expr_df) <- c("gene_id",sample_name$sample_name[1])

library(dplyr)
for (i in 2:nrow(sample_name)){
  
  dfnew <- fread(sample_name$file_name[i], skip = 4, select = 1:2)
  names(dfnew) <- c("gene_id",sample_name$sample_name[i])
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}

setwd("/panfs/compbio/users/zzhu244/volumes")

write.table(expr_df,file = "matrix_counts1.txt",quote = F,row.names = F,sep="\t")



###################GAP2
rm(list=ls())
#bulid date.frame about sample name and file name
sample_name<-read.table("//home/zzhu244/thesis/small_sample_test/GAP2.txt")
colnames(sample_name)[1]<-"sample_name"
sample_name$file_name <- paste(sample_name$sample_name, "ReadsPerGene.out.tab", sep = "")

#read file
setwd("/panfs/compbio/users/zzhu244/volumes/GAP2/qc_output/outcome")
require(data.table)
#expr_df <- fread(sample_name$file_name[30], skip = 4, select = 1:2)
#names(expr_df) <- c("gene_id",sample_name$sample_name[30])
expr_df <- fread(sample_name$file_name[1], skip = 4, select = 1:2)
names(expr_df) <- c("gene_id",sample_name$sample_name[1])

library(dplyr)
for (i in 2:nrow(sample_name)){
  
  dfnew <- fread(sample_name$file_name[i], skip = 4, select = 1:2)
  names(dfnew) <- c("gene_id",sample_name$sample_name[i])
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}

setwd("/panfs/compbio/users/zzhu244/volumes")

write.table(expr_df,file = "matrix_counts2.txt",quote = F,row.names = F,sep="\t")





###################GAP3
rm(list=ls())
#bulid date.frame about sample name and file name
sample_name<-read.table("//home/zzhu244/thesis/small_sample_test/GAP3.txt")
colnames(sample_name)[1]<-"sample_name"
sample_name$file_name <- paste(sample_name$sample_name, "ReadsPerGene.out.tab", sep = "")

#read file
setwd("/panfs/compbio/users/zzhu244/volumes/GAP3/qc_output/outcome")
require(data.table)
#expr_df <- fread(sample_name$file_name[30], skip = 4, select = 1:2)
#names(expr_df) <- c("gene_id",sample_name$sample_name[30])
expr_df <- fread(sample_name$file_name[1], skip = 4, select = 1:2)
names(expr_df) <- c("gene_id",sample_name$sample_name[1])

library(dplyr)
for (i in 2:nrow(sample_name)){
  
  dfnew <- fread(sample_name$file_name[i], skip = 4, select = 1:2)
  names(dfnew) <- c("gene_id",sample_name$sample_name[i])
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}

setwd("/panfs/compbio/users/zzhu244/volumes")

write.table(expr_df,file = "matrix_counts3.txt",quote = F,row.names = F,sep="\t")




###################GAP4
rm(list=ls())
#bulid date.frame about sample name and file name
sample_name<-read.table("//home/zzhu244/thesis/small_sample_test/GAP4.txt")
colnames(sample_name)[1]<-"sample_name"
sample_name$file_name <- paste(sample_name$sample_name, "ReadsPerGene.out.tab", sep = "")

#read file
setwd("/panfs/compbio/users/zzhu244/volumes/GAP4/qc_output/outcome")
require(data.table)
#expr_df <- fread(sample_name$file_name[30], skip = 4, select = 1:2)
#names(expr_df) <- c("gene_id",sample_name$sample_name[30])
expr_df <- fread(sample_name$file_name[1], skip = 4, select = 1:2)
names(expr_df) <- c("gene_id",sample_name$sample_name[1])

library(dplyr)
for (i in 2:nrow(sample_name)){
  
  dfnew <- fread(sample_name$file_name[i], skip = 4, select = 1:2)
  names(dfnew) <- c("gene_id",sample_name$sample_name[i])
  expr_df <- inner_join(expr_df,dfnew,by="gene_id")
}

setwd("/panfs/compbio/users/zzhu244/volumes")

write.table(expr_df,file = "matrix_counts4.txt",quote = F,row.names = F,sep="\t")











