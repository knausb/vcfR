

#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)

#
context("vcfR2hapmap")

##### ##### ##### ##### #####
# vcfR2hapmap

test_that("vcfR2hapmap works",{
  data("vcfR_test")
  
  myHapMap <- vcfR2hapmap(vcfR_test)
  
  expect_true(inherits(myHapMap, "data.frame"))
  expect_true(inherits(myHapMap, "hapMap"))
  expect_true( ncol(myHapMap) > 11 )
  expect_true( nrow(myHapMap) >= 1 )
  
})

