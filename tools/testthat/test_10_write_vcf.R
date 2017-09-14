

#
library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
#
context("write vcf functions")

#data("vcfR_example")
#tot_var <- nrow(vcf@gt)


##### ##### ##### ##### #####
# Setup

myDir <- getwd()
myTmp <- tempdir()

setwd(myTmp)


##### ##### ##### ##### #####
#
# Compiled write_vcf_body works.
#
##### ##### ##### ##### #####




##### ##### ##### ##### #####
#
# Write funcitons work.
#
##### ##### ##### ##### #####

#vcf <- read.vcfR(ex_file, verbose=FALSE)

#test_that("write.vcf works on objects of class vcfR",{
#  write.vcf(vcf, file="test.vcf.gz")
#  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
#  unlink("test.vcf.gz")

#  expect_equal(nrow(vcf), nrow(test))
#  expect_equal(ncol(vcf@fix), ncol(test@fix))
#  expect_equal(ncol(vcf@gt), ncol(test@gt))
#})


#test_that("write.vcf works on objects of class chromR",{
#  suppressWarnings(chrom <- create.chromR(vcf=vcf, seq = dna, ann = gff))

#  write.vcf(chrom, file="test.vcf.gz")
#  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
#  unlink("test.vcf.gz")
  
#  expect_equal(nrow(vcf), nrow(test))
#  expect_equal(ncol(vcf@fix), ncol(test@fix))
#  expect_equal(ncol(vcf@gt), ncol(test@gt))
#})

#test_that("write.vcf works on objects of class chromR when mask=TRUE",{
#  suppressWarnings(chrom <- create.chromR(vcf=vcf, seq = dna, ann = gff))
#  chrom <- masker(chrom, min_DP = 250, max_DP = 750, min_MQ = 59.5, max_MQ = 60.5)
  
#  write.vcf(chrom, file="test.vcf.gz", mask = TRUE)
#  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
#  unlink("test.vcf.gz")

#  expect_equal(sum(chrom@var.info$mask), nrow(test))
#})




##### ##### ##### ##### #####

#setwd( myDir )
#unlink( myTmp, recursive = TRUE )

##### ##### ##### ##### #####
# EOF.
