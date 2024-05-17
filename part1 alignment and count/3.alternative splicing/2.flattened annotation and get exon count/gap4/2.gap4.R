rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")

load("gap4_pre.RData")

se4_1 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[1:50]), singleEnd=TRUE, ignore.strand=TRUE)
save(se4_1,file="gap4_1.RData")

rm(list=ls())
gc()
load("gap4_pre.RData")	

se4_2 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[51:100]), singleEnd=TRUE, ignore.strand=TRUE )
save(se4_2,file="gap4_2.RData")

rm(list=ls())
gc()
load("gap4_pre.RData")	

se4_3 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[101:150]), singleEnd=TRUE, ignore.strand=TRUE )
save(se4_3,file="gap4_3.RData")

rm(list=ls())
gc()
load("gap4_pre.RData")

se4_4 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[151:200]), singleEnd=TRUE, ignore.strand=TRUE )
save(se4_4,file="gap4_4.RData")

rm(list=ls())
gc()
load("gap4_pre.RData")

se4_5 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[201:250]), singleEnd=TRUE, ignore.strand=TRUE )
save(se4_5,file="gap4_5.RData")

rm(list=ls())
gc()
load("gap4_pre.RData")

se4_6 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[251:253]), singleEnd=TRUE, ignore.strand=TRUE )
save(se4_6,file="gap4_6.RData")








