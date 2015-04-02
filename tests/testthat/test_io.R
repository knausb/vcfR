#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("io functions")

# Load data
data(vcfR_example)
pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf_mt <- proc_chrom(pinf_mt, win.size=1000, verbose=FALSE)


# Manage directories.
original_dir <- getwd()
test_dir <- tempdir()



test_that("read/write.vcf works for vcfR objects",{

  setwd(test_dir)
  write.vcf(pinf_vcf, "test.vcf")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
#  expect_is(test@fix$QUAL, "integer") # Old, pure R version.
  expect_is(test@fix$QUAL, "numeric")
  expect_identical(names(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), 371)
  expect_equal(ncol(test@gt), 30)

})


test_that("read/write.vcf works for Chrom objects",{
  
  setwd(test_dir)
  write.vcf(pinf_mt, "test.vcf")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
  #  expect_is(test@fix$QUAL, "integer") # Old, pure R version.
  expect_is(test@fix$QUAL, "numeric")
  expect_identical(names(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), 371)
  expect_equal(ncol(test@gt), 30)
  
})


test_that("write_var_info works for vcfR objects",{
  
  setwd(test_dir)
  write_var_info(pinf_mt, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
  setwd(original_dir)
  
  expect_is(test, "data.frame")
  expect_equal(nrow(test), 371)
#  expect_equal(ncol(test), 371)
  
})





