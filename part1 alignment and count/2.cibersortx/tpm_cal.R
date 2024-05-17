rm(list=ls())

load("beforetpm.RData")

rpk <- expr_df/kb
tpm <- t(t(rpk)/colSums(rpk)*1000000)

save.image("aftertpm.RData")