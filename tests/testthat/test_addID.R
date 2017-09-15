

library(vcfR)
#library(testthat)
context("addID")



test_that("addID works",{
  data("vcfR_example")
  x <- addID(vcf)
  expect_equal( sum( is.na(x@fix[,'ID']) ), 0 )
})


