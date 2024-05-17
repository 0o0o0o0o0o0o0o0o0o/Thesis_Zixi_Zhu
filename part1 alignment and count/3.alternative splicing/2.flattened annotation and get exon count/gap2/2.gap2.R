rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")

load("gap2_pre.RData")

se2_1 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[1:50]), singleEnd=TRUE, ignore.strand=TRUE)
save(se2_1,file="gap2_1.RData")

rm(list=ls())
gc()
load("gap2_pre.RData")	
	
se2_2 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[51:100]), singleEnd=TRUE, ignore.strand=TRUE )
save(se2_2,file="gap2_2.RData")

rm(list=ls())
gc()
load("gap2_pre.RData")	

se2_3 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[101:150]), singleEnd=TRUE, ignore.strand=TRUE )
save(se2_3,file="gap2_3.RData")

rm(list=ls())
gc()
load("gap2_pre.RData")

se2_4 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[151:200]), singleEnd=TRUE, ignore.strand=TRUE )
save(se2_4,file="gap2_4.RData")

rm(list=ls())
gc()
load("gap2_pre.RData")

se2_5 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[201:250]), singleEnd=TRUE, ignore.strand=TRUE )
save(se2_5,file="gap2_5.RData")

rm(list=ls())
gc()
load("gap2_pre.RData")

se2_6 = summarizeOverlaps(
    flattenedAnnotation, BamFileList(bamFiles[251:259]), singleEnd=TRUE, ignore.strand=TRUE )
save(se2_6,file="gap2_6.RData")


