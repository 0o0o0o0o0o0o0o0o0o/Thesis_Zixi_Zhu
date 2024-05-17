library(parallel) 
library(GenomicFeatures)
setwd("/panfs/compbio/users/zzhu244/volumes/")

txdb <- makeTxDbFromGFF("/panfs/compbio/users/sshar28/refData/ucsc/hg38/hg38.refGene.gtf",format="gtf") 
exons_gene <- exonsBy(txdb, by = "gene") ###提取基因外显子


##计算总外显子长度：用reduce去除掉重叠冗余的部分，,width统计长度，最后计算总长度
exons_gene_lens <- parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))}) 
#exons_gene_lens[1:10]
    
##转换为dataframe
geneid_efflen <- data.frame(geneid=names(exons_gene_lens),efflen=as.numeric(exons_gene_lens))
geneid_efflen_gtf1 <- geneid_efflen

save.image("genelength.RData")