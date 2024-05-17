#This R script should be run in HPC
#path: /panfs/compbio/users/zzhu244/

rm(list=ls())
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")
load("annotation.RData")

folders = c("/panfs/compbio/users/zzhu244/volumes/GAP3/qc_output/outcome/")
bamFiles = sapply(folders, function(folder) {
    list.files(folder, pattern = "\\.bam$", full.names = TRUE)
}, USE.NAMES = FALSE)
bamFiles = unlist(bamFiles)
bamFiles = BamFileList(bamFiles)
seqlevelsStyle(flattenedAnnotation) = "UCSC"
save.image("gap3_pre.RData")
