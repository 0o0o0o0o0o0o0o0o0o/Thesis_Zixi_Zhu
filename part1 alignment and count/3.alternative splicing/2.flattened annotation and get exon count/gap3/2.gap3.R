rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")

load("gap3_pre.RData")

se3_1 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[1:50]), singleEnd=TRUE, ignore.strand=TRUE)
save(se3_1,file="gap3_1.RData")

rm(list=ls())
gc()
load("gap3_pre.RData")	

se3_2 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[51:100]), singleEnd=TRUE, ignore.strand=TRUE )
save(se3_2,file="gap3_2.RData")

rm(list=ls())
gc()
load("gap3_pre.RData")	

se3_3 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[101:150]), singleEnd=TRUE, ignore.strand=TRUE )
save(se3_3,file="gap3_3.RData")

rm(list=ls())
gc()
load("gap3_pre.RData")	

se3_4 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[151:200]), singleEnd=TRUE, ignore.strand=TRUE )
save(se3_4,file="gap3_4.RData")

rm(list=ls())
gc()
load("gap3_pre.RData")	
	
se3_5 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[201:241]), singleEnd=TRUE, ignore.strand=TRUE )
save(se3_5,file="gap3_5.RData")
	
