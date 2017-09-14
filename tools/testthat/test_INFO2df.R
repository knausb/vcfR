

#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)


#
context("INFO2df")


test_that("INFO2df example works",{
  data("vcfR_test")
  myDf <- INFO2df(vcfR_test)
  expect_true(class(myDf) == 'data.frame')
})


test_that("metaINFO2df example works",{
  data("vcfR_test")
  myDf <- metaINFO2df(vcfR_test)
  expect_true(class(myDf) == 'data.frame')
})