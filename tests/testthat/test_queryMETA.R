
#
library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("queryMETA")

##### ##### ##### ##### #####

test_that("queryMETA works",{
  data(vcfR_test)
  myMeta <- queryMETA(vcfR_test)

  expect_is(myMeta, "character")
})


test_that("queryMETA element works",{
  data(vcfR_test)
  myMeta <- queryMETA(vcfR_test, element = 'DP')

  expect_is(myMeta, "list")
})


