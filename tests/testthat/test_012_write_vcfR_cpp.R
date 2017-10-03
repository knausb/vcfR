
#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
#
context("write vcf C++ functions")


##### ##### ##### ##### #####
#
# Compiled write_vcf_body works.
#
##### ##### ##### ##### #####

myDir <- getwd()
myTmp <- tempdir()
setwd(myTmp)

test_that("write_vcf_body works on objects of class vcfR",{
  data(vcfR_test)
  test <- .write_vcf_body(fix = vcfR_test@fix,
                          gt = vcfR_test@gt, 
                          filename = 'myVars.vcf.gz')

  test2 <- scan('myVars.vcf.gz', what='character', sep = '\n', quiet = TRUE)
  unlink('myVars.vcf.gz')

  expect_true( is.null(test) )
  expect_equal( nrow(vcfR_test), length(test2) )
})


test_that("write_vcf_body works on objects of class vcfR with no variants",{
  data(vcfR_test)
  vcfR_test@fix <- vcfR_test@fix[0,]
  vcfR_test@gt <- vcfR_test@gt[0,]

  test <- .write_vcf_body(fix = vcfR_test@fix,
                          gt = vcfR_test@gt, 
                          filename = 'myVars.vcf.gz')
    
  test2 <- scan('myVars.vcf.gz', what='character', sep = '\n', quiet = TRUE)
  unlink('myVars.vcf.gz')
  
  expect_equal( length(test2), 0 )
})

  
test_that("write_vcf_body works on objects of class vcfR with no GT region",{
  data(vcfR_test)
  vcfR_test@gt <- vcfR_test@gt[0,0]
  
  test <- .write_vcf_body(fix = vcfR_test@fix,
                          gt = vcfR_test@gt, 
                          filename = 'myVars.vcf.gz')
  
  test2 <- read.table('myVars.vcf.gz', sep = '\t')
  unlink('myVars.vcf.gz')

  expect_equal( nrow(vcfR_test), nrow(test2) )
  expect_equal( as.numeric( ncol( vcfR_test ) ), 8 )
})



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
# EOF.