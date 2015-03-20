# detach(package:vcfR, unload=T)
library(vcfR)
context("gt_to_popsum functions")

data(vcfR_example)


pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)


pinf_mt <- gt_to_popsum(pinf_mt)

head(pinf_mt@var.info)


