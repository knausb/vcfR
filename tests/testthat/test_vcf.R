
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("vcf functions")

ex_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")


test_that("compiled input functions work",{
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  expect_equal(length(stats), 4)
  expect_is(stats, "numeric")
  
  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
  expect_equal(length(meta), as.numeric(stats["meta"]))
  expect_is(meta, "character")
  
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
  expect_is(body, "matrix")
  expect_equal(nrow(body), as.numeric(stats["variants"]))
  expect_equal(ncol(body), as.numeric(stats["columns"]))
})


test_that("read.vcf works",{
  vcf <- read.vcf(ex_file)
  expect_is(vcf, "vcfR")
  expect_is(vcf@meta, "character")
  expect_is(vcf@fix, "matrix")
  expect_is(vcf@gt, "matrix")
})


vcf <- read.vcf(ex_file)

test_that("write.vcf works",{
  
  
})
  
  
#debug(read.vcf)


