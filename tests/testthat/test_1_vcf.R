
#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("vcf functions")

#ex_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")

data("vcfR_example")
tot_var <- nrow(vcf@gt)

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
  
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)
  expect_is(body, "matrix")
  expect_equal(nrow(body), as.numeric(stats["variants"]))
  expect_equal(ncol(body), as.numeric(stats["columns"]))
})


test_that("compiled vcfR_read_body works when file contains no variants",{
  vcf2 <- vcf
  vcf2@fix <- vcf2@fix[0,]
  vcf2@gt <- vcf2@gt[0,]
  
  write.vcf(vcf2, ex_file)
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
#  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = 0, skip = 0, cols=1:stats['columns'], 0)

  unlink(ex_file)

})


##### ##### ##### ##### #####


data("vcfR_example")
write.vcf(vcf, file=ex_file)


test_that("read.vcfR works",{
  vcf <- read.vcfR(ex_file, verbose=FALSE)
  expect_is(vcf, "vcfR")
  expect_is(vcf@meta, "character")
  expect_is(vcf@fix, "matrix")
  expect_is(vcf@gt, "matrix")
})


test_that("read.vcfR nrows works",{
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, nrows=count)
  expect_equal(nrow(vcf), count)
})


test_that("read.vcfR skip works",{
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, skip=count)
  expect_equal(nrow(vcf), tot_var - count)
})


test_that("read.vcfR column selection works",{
  vcf <- read.vcfR(ex_file, verbose=FALSE, cols=11:12)
  expect_is(vcf@gt, "matrix")
  expect_equal(ncol(vcf@fix), 8)
  expect_equal(ncol(vcf@gt), 3)
})


test_that("read.vcfR works for vcf files which contain no variants",{  
  vcf2 <- vcf
  vcf2@fix <- vcf2@fix[0,]
  vcf2@gt <- vcf2@gt[0,]
  
  write.vcf(vcf2, ex_file)
  test <- read.vcfR(ex_file, verbose=FALSE)
  unlink(ex_file)


  expect_equal(ncol(test@fix), ncol(vcf2@fix))
  expect_equal(ncol(test@gt), ncol(vcf2@gt))
  expect_equal(nrow(test@fix), nrow(vcf2@fix))
  expect_equal(nrow(test@gt), nrow(vcf2@gt))
})


##### ##### ##### ##### #####

#vcf <- read.vcfR(ex_file, verbose=FALSE)

test_that("write.vcf works",{
  write.vcf(vcf, file=ex_file)  
  unlink(ex_file)
})


##### ##### ##### ##### #####


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

  
#debug(read.vcfR)

#unlink(ex_file)

##### ##### ##### ##### #####
# EOF.