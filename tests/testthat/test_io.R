#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("io functions")

data(vcfR_example)

# Manage directories.
original_dir <- getwd()
test_dir <- tempdir()

# Load data
data(vcfR_example)

test_that("vcf file io works",{
  setwd(test_dir)

  write.vcf(pinf_vcf, "test.vcf")
  test <- read.vcf("test.vcf")
  
  expect_is(test, "vcfR")
  expect_identical(names(test@fix)[1], "CHROM")
  
  unlink("test.vcf")
  
  setwd(original_dir)
})


