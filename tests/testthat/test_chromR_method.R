

# detach(package:vcfR, unload=T)
#library(testthat)
library(vcfR)
context("chromR methods")


test_that("chromR show",{
  data("chromR_example")
  capture.output( tmp <- show(chrom) )
  expect_null(tmp)
})


##### ##### ##### ##### #####

test_that("chromR length",{
  data("chromR_example")
  expect_equal(length(chrom), chrom@len)
})


##### ##### ##### ##### #####

test_that("chromR print",{
#  tmp <- print(chrom)
#  expect_true(is.null(tmp))
})

test_that("chromR head",{
  data("chromR_example")
#  tmp <- head(chrom)
#  expect_null(tmp)
#  expect_true(is.null(tmp))
})


test_that("chromR names<-",{
  data("chromR_example")
  names(chrom) <- "bob"
  expect_identical(names(chrom), "bob")
})


# EOF.

