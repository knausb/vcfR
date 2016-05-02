
#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("genotype matrix functions")


##### ##### ##### ##### #####
#
# alleles2consensus
#
##### ##### ##### ##### #####

test_that("alleles2consensus works",{

gt <- structure(c("T/T", "./.", NA, "./G", "T/.", "A/A", "C/C", 
"A/A", "T/T", "T/T", "C/C", "T/T", "C/C", "T/T", "G/G", "T/T", 
"A/A", "C/C", "A/A", "T/T", "T/T", "C/C", "T/T", "C/C", "T/T", 
"G/G", "T/T", "T/T", "C/C", "A/A", "T/T", "T/T", "C/C", "T/T", 
"T/T", "T/T", "A/A", "T/G", "T/T", "C/C", "C/C", "T/T", "T/T", 
"C/C"), .Dim = c(11L, 4L), .Dimnames = list(NULL, c("1.Anne", 
"1.Nadine", "2.Heidi", "3.Agatha")))

  gt <- alleles2consensus(gt, NA_to_n = TRUE)
  
  expect_true(gt[2,1] == "n")
  expect_true(gt[3,1] == "n")
  expect_true(gt[4,1] == "n")
  expect_true(gt[5,1] == "n")
})


