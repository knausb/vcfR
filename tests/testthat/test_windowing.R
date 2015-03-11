# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("windowing functions")

data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)


gq <- extract.gt(pinf_vcf, element="GQ", as.numeric=TRUE)

gqw <- windowize_NM(gq, pos=pinf_mt@var.info$POS, starts=pinf_mt@win.info$start, ends=pinf_mt@win.info$end)

head(gqw)
nrow(gqw)


