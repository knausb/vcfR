
context("vcfR methods")

#library(testthat)
library(vcfR)
data("vcfR_example")



##### ##### ##### ##### #####
# rbind



test_that("rbind works",{

  vcf2 <- rbind( vcf )
  expect_equal( 1 * nrow(vcf@fix), nrow(vcf2@fix) )
  
#  vcf2 <- rbind( vcf, vcf )
  vcf2 <- rbind( vcf, vcf, deparse.level = 0 )
  expect_equal( 2 * nrow(vcf@fix), nrow(vcf2@fix) )
  
#  vcf2 <- rbind( vcf, vcf, vcf, vcf )
#  expect_equal( 4 * nrow(vcf@fix), nrow(vcf2@fix) )
  
})


test_that("nrow works",{
  expect_equal( nrow(vcf@fix), nrow(vcf) )
})


##### ##### ##### ##### #####
# EOF.