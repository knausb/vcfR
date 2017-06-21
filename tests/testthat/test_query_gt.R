
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

##### ##### ##### ##### #####
# is.biallelic


test_that("is.biallelic works",{
  data("vcfR_test")
  myVars <- is.biallelic(vcfR_test)
  
  expect_equal(sum(myVars), 3)  
})


