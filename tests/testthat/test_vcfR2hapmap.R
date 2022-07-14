

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


test_that("vcfR2hapmap works, missing data",{
  data("vcfR_test")
  
  vcfR_test@gt[1, 3] <- "./.:48:4:51,51"
  vcfR_test@gt[2, 3] <- ".|.:48:4:51,51"
  
  myHapMap <- vcfR2hapmap(vcfR_test)
  expect_true(myHapMap[2, 13] == "NN")
  expect_true(myHapMap[3, 13] == "NN")
  
})


