
#library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)

context("read vcfR R functions")

##### ##### ##### ##### #####


test_that("read.vcfR works",{
  data("vcfR_example")
  ex_file <- "myFile.vcf.gz"
  write.vcf(vcf, file=ex_file)
  vcf <- read.vcfR(ex_file, verbose=FALSE)
  unlink(ex_file)
  
  expect_is(vcf, "vcfR")
  expect_is(vcf@meta, "character")
  expect_is(vcf@fix, "matrix")
  expect_is(vcf@gt, "matrix")
})


test_that("read.vcfR nrows works",{
  data("vcfR_example")
  ex_file <- "myFile.vcf.gz"
  write.vcf(vcf, file=ex_file)
  
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, nrows=count)
  expect_equal(nrow(vcf), count)
  unlink(ex_file)
})


test_that("read.vcfR skip works",{
  data("vcfR_example")
  ex_file <- "myFile.vcf.gz"
  write.vcf(vcf, file=ex_file)
  tot_var <- nrow(vcf)
  
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, skip=count)
  expect_equal(nrow(vcf), tot_var - count)
  unlink(ex_file)
})


test_that("read.vcfR column selection works",{
  data("vcfR_example")
  ex_file <- "myFile.vcf.gz"
  write.vcf(vcf, file=ex_file)
  vcf <- read.vcfR(ex_file, verbose=FALSE, cols=c(9,11:12))
  
  expect_is(vcf@gt, "matrix")
  expect_equal(ncol(vcf@fix), 8)
  expect_equal(ncol(vcf@gt), 3)
  unlink(ex_file)
})


test_that("read.vcfR works for vcf files which contain no variants",{  
  data("vcfR_example")
  ex_file <- "myFile.vcf.gz"
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


test_that("read.vcfR works when file contains one variant",{
  data(vcfR_test)
  ex_file <- "myFile.vcf.gz"
  vcfR_test <- vcfR_test[1,]

  write.vcf(vcfR_test, ex_file)
  vcfR_test2 <- read.vcfR(ex_file, verbose=FALSE)
  
  expect_equal( nrow(vcfR_test2@fix), 1)
  unlink(ex_file)
})


test_that("read.vcfR works when file contains one variant, no meta",{
  data(vcfR_test)
  ex_file <- "myFile.vcf.gz"
  vcfR_test <- vcfR_test[1,]
  vcfR_test@meta <- vector(mode='character', length=0)

  write.vcf(vcfR_test, ex_file)
  vcfR_test2 <- read.vcfR(ex_file, checkFile = FALSE, verbose=FALSE)
  
  expect_equal( nrow(vcfR_test2@fix), 1)
  
  unlink(ex_file)
})


test_that("read.vcfR verbose works",{
  data(vcfR_test)
  test_file <- "test.vcf.gz"
  
  write.vcf(x = vcfR_test, file = test_file)
  testMessage <- capture.output(read.vcfR(test_file, verbose=TRUE))
  unlink(test_file)

  expect_equal( grep("File attributes:", testMessage), 2)
  expect_equal( grep("  meta lines:", testMessage), 3)
  expect_equal( grep("  header_line:", testMessage), 4)
  expect_equal( grep("  variant count:", testMessage), 5)
  expect_equal( grep("  column count:", testMessage), 6)
})


test_that("VCF with no GT",{
  data(vcfR_test)
  vcf@gt <- matrix("a", nrow=0, ncol=0)
  ex_file <- "test.vcf.gz"
  
  test <- write.vcf(vcf, file = ex_file)
  expect_equal(test, NULL)
  
#  debug(read.vcfR)
  test <- read.vcfR(ex_file, verbose = FALSE)
  unlink("test.vcf.gz")
    
  expect_is(test, "vcfR")  
  expect_is(test@fix, "matrix")
  
  expect_equal(nrow(test@fix), nrow(vcf@fix))
  expect_equal(ncol(test@fix), 8)
  
  expect_is(test@gt, "matrix")
  expect_equal( ncol(test@gt), 0 )
  expect_equal( nrow(test@gt), 0 )
})


test_that("read.vcfR works for files in other directories",{
  data("vcfR_example")
  test_dir <- tempdir()
  setwd(test_dir)
  
  
  if( !dir.exists('subdir') ){
    dir.create('subdir')
  }
  
  setwd('subdir')
  write.vcf(vcf, "test.vcf.gz")
  setwd(test_dir)
  vcf1 <- read.vcfR("./subdir/test.vcf.gz", verbose = FALSE)
  unlink("./subdir/test.vcf.gz")
  
  expect_equal(nrow(vcf@fix), nrow(vcf1@fix))
})


##### ##### ##### ##### #####



##### ##### ##### ##### #####
# EOF.