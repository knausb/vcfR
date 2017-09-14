
# detach(package:vcfR, unload=T)
#library(testthat)
library(vcfR)
context("create.chromR functions")

#library(testthat)
#data(vcfR_example)

test_that("Create a null chromR",{
  data("vcfR_example")
  chrom <- methods::new(Class="chromR")
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix), 0)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})



