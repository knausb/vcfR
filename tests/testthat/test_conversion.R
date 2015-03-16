# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("windowing functions")

data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)

gt <- extract.gt(pinf_mt, element="PL", as.numeric=TRUE)

test_that("gt is numeric",{
  expect_is(gt, "matrix")
  expect_equal(is.numeric(gt), TRUE)
})


#head(pinf_mt@vcf.gt)

#head(gt)


