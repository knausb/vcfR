library('vcfR')
library('testthat')

context("getFIX accessors")



#
data("vcfR_test")
#
chrom <- create.chromR(vcfR_test, verbose = FALSE)

test_that("getFIX returns a matrix", {
  expect_is(getFIX(vcfR_test), "matrix")
  expect_is(getFIX(chrom), "matrix")
})

test_that("getCHROM returns a character", {
  expect_is(getCHROM(vcfR_test), "character")
  expect_is(getCHROM(chrom), "character")
})

test_that("getPOS returns an integer", {
  expect_is(getPOS(vcfR_test), "integer")
  expect_is(getPOS(chrom), "integer")
})

test_that("getQUAL returns a numeric", {
  expect_is(getQUAL(vcfR_test), "numeric")
  expect_is(getQUAL(chrom), "numeric")
})

test_that("getALT returns a character", {
  expect_is(getALT(vcfR_test), "character")
  expect_is(getALT(chrom), "character")
})

test_that("getREF returns a character", {
  expect_is(getREF(vcfR_test), "character")
  expect_is(getREF(chrom), "character")
})

test_that("getID returns a character", {
  expect_is(getID(vcfR_test), "character")
  expect_is(getID(chrom), "character")
})

test_that("getFILTER returns a character", {
  expect_is(getFILTER(vcfR_test), "character")
  expect_is(getFILTER(chrom), "character")
})

