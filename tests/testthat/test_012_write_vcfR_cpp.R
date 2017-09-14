
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
#  test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', 
#                fix = vcfR_test@fix, gt = vcfR_test@gt, 
#                filename = 'myVars.vcf.gz', mask = 0)
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

#  test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR',
#                fix = vcfR_test@fix, gt = vcfR_test@gt, 
#                filename = 'myVars.vcf.gz', mask = 0)
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
  
#  test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR',
#                fix = vcfR_test@fix, gt = vcfR_test@gt, 
#                filename = 'myVars.vcf.gz', mask = 0)
  test <- .write_vcf_body(fix = vcfR_test@fix,
                          gt = vcfR_test@gt, 
                          filename = 'myVars.vcf.gz')
  
  test2 <- read.table('myVars.vcf.gz', sep = '\t')
  unlink('myVars.vcf.gz')

  expect_equal( nrow(vcfR_test), nrow(test2) )
  expect_equal( as.numeric( ncol( vcfR_test ) ), 8 )
})


