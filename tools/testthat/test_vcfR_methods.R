

#
library(testthat)
context("vcfR methods")
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
# []

test_that("[samples = numeric]",{
  data("vcfR_test")
  vcf <- vcfR_test[samples=c(1,3)]
  expect_equal(colnames(vcf@gt), colnames(vcfR_test@gt)[c(1,2,4)])
})

test_that("[samples = integer]",{
  data("vcfR_test")
  vcf <- vcfR_test[samples=as.integer(c(1,3))]
  expect_equal(colnames(vcf@gt), colnames(vcfR_test@gt)[c(1,2,4)])
})

test_that("[samples = character]",{
  data("vcfR_test")
  vcf <- vcfR_test[samples=c('NA00001', 'NA00003')]
  expect_equal(colnames(vcf@gt), colnames(vcfR_test@gt)[c(1,2,4)])
})

test_that("[samples = logical]",{
  data("vcfR_test")
  vcf <- vcfR_test[samples=c(TRUE, FALSE, TRUE)]
  expect_equal(colnames(vcf@gt), colnames(vcfR_test@gt)[c(1,2,4)])
})


##### ##### ##### ##### #####
# EOF.
