#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/
#The aim is to run function summarizeOverlaps by part


rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")

load("gap1_pre.RData")

se1_1 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[1:50]), singleEnd=TRUE, ignore.strand=TRUE)
save(se1_1,file="gap1_1.RData")

rm(list=ls())
gc()
load("gap1_pre.RData")	

se1_2 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[51:100]), singleEnd=TRUE, ignore.strand=TRUE )
save(se1_2,file="gap1_2.RData")

rm(list=ls())
gc()
load("gap1_pre.RData")		

se1_3 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[101:150]), singleEnd=TRUE, ignore.strand=TRUE )
save(se1_3,file="gap1_3.RData")

rm(list=ls())
gc()
load("gap1_pre.RData")	
	
se1_4 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[151:183]), singleEnd=TRUE, ignore.strand=TRUE )
save(se1_4,file="gap1_4.RData")




