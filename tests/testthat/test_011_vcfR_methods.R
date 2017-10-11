
#library(testthat)
library(vcfR)

context("vcfR methods")


##### ##### ##### ##### #####
# rbind

test_that("rbind works",{
  data("vcfR_example")
  vcf2 <- rbind( vcf )
  expect_equal( 1 * nrow(vcf@fix), nrow(vcf2@fix) )
  
  vcf2 <- rbind( vcf, vcf, deparse.level = 0 )
  expect_equal( 2 * nrow(vcf@fix), nrow(vcf2@fix) )
})


##### ##### ##### ##### #####
# nrow

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


test_that("show no gt",{
  data("vcfR_test")
  vcfR_test@gt <- matrix(nrow=0, ncol=0)
  myTest <-  capture.output(head(vcfR_test), type = c("output", "message"))
  expect_equal(grep('No gt slot present', myTest), 23)
})



##### ##### ##### ##### #####
# EOF.