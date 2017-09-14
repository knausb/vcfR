
# detach(package:vcfR, unload=T)
#
library(testthat)
library(vcfR)
context("query_gt")


##### ##### ##### ##### #####
# is.polymorphic 


test_that("is.polymorphic works",{
  data("vcfR_test")
  myVars <- is.polymorphic(vcfR_test)
  
  expect_equal(sum(myVars), nrow(vcfR_test))  
})

test_that("is.polymorphic, not vcfR warning",{
  myMat <- matrix(1, ncol=3, nrow=4)
  expect_error(is.polymorphic(myMat), "Expected an object of class vcfR")
})


test_that("is.polymorphic, na.omit == TRUE",{
  data("vcfR_test")
  myVars <- is.polymorphic(vcfR_test, na.omit = TRUE)
  
  expect_equal(sum(myVars), nrow(vcfR_test))
})


##### ##### ##### ##### #####
# is.biallelic


test_that("is.biallelic works",{
  data("vcfR_test")
  myVars <- is.biallelic(vcfR_test)
  
  expect_equal(sum(myVars), 3)  
})


