


library(vcfR)
#library(testthat)
context("maf")


test_that("maf works on vcfR",{
  data("vcfR_example")
  my.maf <- maf(vcf)

  expect_equal( ncol(my.maf), 4 )
})




