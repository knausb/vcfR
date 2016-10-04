

library(vcfR)
library(testthat)
context("addID")

data("vcfR_example")


test_that("addID works",{
  x <- addID(vcf)
  expect_equal( sum( is.na(x@fix[,'ID']) ), 0 )
})


