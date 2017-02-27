

# detach(package:vcfR, unload=T)
#
library(testthat)
library(vcfR)
context("chromR methods")

data("vcfR_example")

chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, verbose=FALSE)

test_that("chromR show",{
  capture.output( tmp <- show(chrom) )
  expect_null(tmp)
})


##### ##### ##### ##### #####

test_that("chromR length",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, verbose=FALSE)
  expect_equal(length(chrom), chrom@len)
  
#  tmp <- print(chrom)
#  expect_true(is.null(tmp))
})


##### ##### ##### ##### #####

test_that("chromR print",{
#  tmp <- print(chrom)
#  expect_true(is.null(tmp))
})

test_that("chromR head",{
#  tmp <- head(chrom)
#  expect_null(tmp)
#  expect_true(is.null(tmp))
})


test_that("chromR names<-",{
  names(chrom) <- "bob"
  expect_identical(names(chrom), "bob")
#  expect_true(is.null(tmp))
})


# EOF.

