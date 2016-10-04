


library(vcfR)
library(testthat)
context("maf")

data("vcfR_example")


test_that("maf works on vcfR",{
  my.maf <- maf(vcf)

  expect_equal( ncol(my.maf), 4 )
})




