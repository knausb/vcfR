

#detach(package:vcfR, unload=TRUE)
library(vcfR)
#library(testthat)
context("write.vcf R functions")


##### ##### ##### ##### #####
# write.vcf


test_that("read/write.vcf works for vcfR objects",{
  data("vcfR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  write.vcf(vcf, "test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
#  setwd(original_dir)
  
  expect_is(test, "vcfR")

  expect_identical(colnames(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), nrow(vcf@gt))
  expect_equal(ncol(test@gt), ncol(vcf@gt))
})



test_that("write.vcf APPEND=TRUE does not include header",{
  data("vcfR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  write.vcf(vcf, "test.vcf.gz", APPEND=TRUE)
  test <- read.vcfR("test.vcf.gz", checkFile = FALSE, verbose = FALSE)
  unlink("test.vcf.gz")
#  setwd(original_dir)
})


test_that("write.vcf.gz works for Chrom objects",{
  data("chromR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  write.vcf(chrom, "test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
#  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_identical(colnames(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), nrow(vcf@gt))
  expect_equal(ncol(test@gt), ncol(vcf@gt))
})


test_that("write.vcf.gz works for Chrom objects with mask",{
  data("chromR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  chrom@var.info$mask <- FALSE
  chrom@var.info$mask[1:50] <- TRUE
  
  write.vcf(chrom, "test.vcf.gz", mask=TRUE)
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
#  setwd(original_dir)
  chrom@var.info$mask <- TRUE

  expect_is(test, "vcfR")
  expect_identical(colnames(test@fix)[1], "CHROM")
  expect_equal(ncol(test@gt), ncol(vcf@gt))
  expect_equal( nrow(test@fix), 50 )
})


test_that("write.var.info works for chromR objects",{
  data("chromR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  write.var.info(chrom, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
#  setwd(original_dir)
  
  expect_is(test, "data.frame")
  expect_equal(nrow(test), nrow(vcf@fix))
  expect_equal(length(grep("CHROM", colnames(test))), 1)
  expect_equal(length(grep("POS", colnames(test))), 1)
  expect_equal(length(grep("mask", colnames(test))), 1)  
})


test_that("write.win.info works for Chrom objects",{
  data("chromR_example")
  chrom <- proc.chromR(chrom, verbose = FALSE)
  test_dir <- tempdir()
  setwd(test_dir)
  write.win.info(chrom, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
#  setwd(original_dir)

  expect_is(test, "data.frame")
  expect_equal(nrow(test), nrow(chrom@win.info))
#  expect_equal(ncol(test), 12)
  expect_equal(grep("CHROM", names(test), value=TRUE), "CHROM")
  expect_equal(grep("window", names(test), value=TRUE), "window")
  expect_equal(grep("start", names(test), value=TRUE), "start")
  expect_equal(grep("end", names(test), value=TRUE), "end")
#  expect_equal(length(grep("window", names(test))), 1)
#  expect_equal(length(grep("start", names(test))), 1)
#  expect_equal(length(grep("end", names(test))), 1)
})



##### ##### ##### ##### #####
# EOF.