# extractgt devel

#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("extract_gt functions")

data(vcfR_example)


pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)


gt <- extract.gt(pinf_vcf, element="GT", as.numeric=FALSE)
pl <- extract.gt(pinf_vcf, element="PL", as.numeric=FALSE)
gq <- extract.gt(pinf_vcf, element="GQ", as.numeric=TRUE)


test_that("gt, pl ad gq are matrices",{
  expect_is(gt, "matrix")
  expect_is(pl, "matrix")
  expect_is(gq, "matrix")
})

test_that("gq is numeric",{
  expect_is(gq, "matrix")
  expect_equal(is.numeric(gq), TRUE)
})




#head(gt)
#head(pl)
#head(gq)

#head(pinf_mt@vcf.gt)


#head(pinf_vcf)


#ncol(pinf_vcf@gt)

#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GT")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="PL")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="DP")


#outm <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', pinf_vcf@gt, element="DP")

#outm <- extract.gt2(pinf_vcf, element="DP", as.numeric=F)
#outm <- extract.gt2(pinf_vcf, element="GQ", as.numeric=T)

#head(pinf_vcf@gt)
#head(outm)



#outm <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")

#head(outm)


#outm2 <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outm)
#head(outm2)

