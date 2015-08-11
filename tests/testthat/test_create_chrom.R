# create_chrom tests.

# 
detach(package:vcfR, unload=T)
library(vcfR)
context("create_chrom functions")

#data(vcfR_example)

test_that("we can create a Chrom",{
  x <- new(Class="Chrom")
  expect_is(x, "Chrom")
  
  

})




pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)

test_that("we can create a Chrom",{
  expect_that(pinf_mt, is_a("Chrom"))  
})


pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, verbose=F)

test_that("we can create a Chrom, no annotation",{
  expect_that(pinf_mt, is_a("Chrom"))  
})


pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, ann=pinf_gff, verbose=F)

test_that("we can create a Chrom, no sequence",{
  expect_that(pinf_mt, is_a("Chrom"))
})


pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, verbose=F)

test_that("we can create a Chrom, no sequence or annotation",{
  expect_that(pinf_mt, is_a("Chrom"))
})


