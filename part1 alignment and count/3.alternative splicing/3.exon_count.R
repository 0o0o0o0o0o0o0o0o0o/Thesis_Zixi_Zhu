#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#aim is to merge exon count matrix together

rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")
load("annotation.RData")

load("gap1_1.RData")
load("gap1_2.RData")
count <- cbind(se1_1,se1_2)
load("gap1_3.RData")
count <- cbind(count,se1_3)
load("gap1_4.RData")
count <- cbind(count,se1_4)

load("gap2_1.RData")
load("gap2_2.RData")
load("gap2_3.RData")
load("gap2_4.RData")
load("gap2_5.RData")
load("gap2_6.RData")
count <- cbind(count,se2_1)
count <- cbind(count,se2_2)
count <- cbind(count,se2_3)
count <- cbind(count,se2_4)
count <- cbind(count,se2_5)
count <- cbind(count,se2_6)

load("gap3_1.RData")
load("gap3_2.RData")
load("gap3_3.RData")
load("gap3_4.RData")
load("gap3_5.RData")
count <- cbind(count,se3_1)
count <- cbind(count,se3_2)
count <- cbind(count,se3_3)
count <- cbind(count,se3_4)
count <- cbind(count,se3_5)

load("gap4_1.RData")
load("gap4_2.RData")
load("gap4_3.RData")
load("gap4_4.RData")
load("gap4_5.RData")
load("gap4_6.RData")
count <- cbind(count,se4_1)
count <- cbind(count,se4_2)
count <- cbind(count,se4_3)
count <- cbind(count,se4_4)
count <- cbind(count,se4_5)
count <- cbind(count,se4_6)
#936
colnames(count) = substr(colnames(count),1,6)


setwd("/panfs/compbio/users/zzhu244/volumes")

library("readxl")
set2 <- read_excel("DepressionGenesNetworks_case_control_phenotype_distribution.xlsx", sheet = "Data", col_names = T)

name <- set2$Cell_id
name <- name[name!="LD0053"]

count <- count[,colnames(count)%in%name]
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")
save(count,file="count.RData")





