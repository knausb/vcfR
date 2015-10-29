
#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("vcf functions")

#ex_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")

data("vcfR_example")

##### ##### ##### ##### #####
# Manage directories.

#original_dir <- getwd()
test_dir <- tempdir()

ex_file <- paste(test_dir, "/test.vcf.gz", sep="")

write.vcf(vcf, file=ex_file)


##### ##### ##### ##### #####

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
  vcf <- read.vcf(ex_file, verbose=FALSE)
  expect_is(vcf, "vcfR")
  expect_is(vcf@meta, "character")
  expect_is(vcf@fix, "matrix")
  expect_is(vcf@gt, "matrix")
})


vcf <- read.vcf(ex_file, verbose=FALSE)

test_that("write.vcf works",{
  orig.dir <- getwd()
  temp.dir <- tempdir()
  setwd(temp.dir)
  write.vcf(vcf, file="temp.vcf.gz")  
  unlink("temp.vcf.gz")
  setwd(orig.dir)
})


test_that("vcfR subsetters works",{
  # Rows
  vcf2 <- vcf[1:10,]
  expect_equal(nrow(vcf2@fix), 10)
  expect_equal(nrow(vcf2@gt), 10)
  
  # Columns
  vcf2 <- vcf[,1:4]
  expect_equal(ncol(vcf2@fix), 8)
  expect_equal(ncol(vcf2@gt), 4)
})

  
#debug(read.vcf)

unlink(ex_file)

##### ##### ##### ##### #####
# EOF.