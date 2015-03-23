# detach(package:vcfR, unload=T)
library(vcfR)
context("gt_to_popsum functions")

data(vcfR_example)


pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)


pinf_mt <- gt_to_popsum(pinf_mt)

head(pinf_mt@var.info)
tail(pinf_mt@var.info)


gt <- extract.gt(pinf_mt, element="GT")
is.na(gt[1,1:4]) <- TRUE
gt[371,29] <- "2/2"
gt[371,28] <- "3/2"

pinf_mt@var.info <- .Call('vcfR_gt_to_popsum', PACKAGE = 'vcfR', var_info=pinf_mt@var.info, gt=gt)
head(pinf_mt@var.info)
tail(pinf_mt@var.info)



