

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("is.het")


test_that("is.hets, na_is_false=TRUE returns all logicals",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test)
  hets <- is.het(gt, na_is_false = TRUE)
  
  expect_is(hets[,1], 'logical')
  expect_is(hets[,2], 'logical')
  expect_is(hets[,3], 'logical')
})


test_that("is.hets, na_is_false=FALSE returns NAs",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test)
  gt[1,1] <- NA
  gt[1,2] <- "1|."
  gt[1,3] <- "./1"
  
  hets <- is.het(gt, na_is_false = FALSE)
  
  expect_equal( sum( is.na(hets[1,]) ), 3)
})


##### ##### ##### ##### #####
#
# Compiled version
#
##### ##### ##### ##### #####


test_that("is_hets, na_is_false=TRUE returns all logicals",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test)
  
  hets <- is_het(gt, na_is_false = TRUE)
  expect_is(hets[,1], 'logical')
  expect_is(hets[,2], 'logical')
  expect_is(hets[,3], 'logical')
})


test_that("is_hets, na_is_false=TRUE does not return NAs",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test)
  
  gt[1,1] <- "."
  is.na(gt[2,1]) <- TRUE
  gt[1,2] <- "1|."
  gt[1,3] <- "./1"
  
  hets <- is_het(gt, na_is_false = TRUE)
  
  expect_equal( sum( is.na(hets[1,]) ), 0 )
  
})


test_that("is_hets, na_is_false=FALSE returns NAs",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test)
  
  gt[1,1] <- "."
  is.na(gt[2,1]) <- TRUE
  gt[1,2] <- "1|."
  gt[1,3] <- "./1"
  
  hets <- is_het(gt, na_is_false = FALSE)
  
  expect_equal( sum( is.na(hets[1,]) ), 3)

})


