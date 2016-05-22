

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("is.het")

data(vcfR_test)
gt <- extract.gt(vcfR_test)

test_that("is.hets, na_is_false=TRUE returns all logicals",{
  hets <- is.het(gt, na_is_false = TRUE)
  expect_is(hets[,1], 'logical')
  expect_is(hets[,2], 'logical')
  expect_is(hets[,3], 'logical')
})


test_that("is.hets, na_is_false=FALSE returns NAs",{
  gt[1,1] <- NA
  gt[1,2] <- "1|."
  gt[1,3] <- "./1"
  
  hets <- is.het(gt, na_is_false = FALSE)
  
  expect_equal( sum( is.na(hets[1,]) ), 3)

})


