#This R script should be run in HPC


rm(list=ls())
library(GenomicFeatures)
txdb = makeTxDbFromGFF("/panfs/compbio/users/sshar28/refData/ucsc/hg38/hg38.refGene.gtf")

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) =
    sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
setwd("/panfs/compbio/users/zzhu244/volumes/exon/")	
save.image("annotation.RData")





