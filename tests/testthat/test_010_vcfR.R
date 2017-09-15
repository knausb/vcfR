
#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)

<<<<<<< HEAD

#
context("vcf functions")
=======
#
context("vcfR functions")

##### ##### ##### ##### #####
#
# Tests for functions that work on vcfR objects.
# Tests for reading and writing vcfR objects
# occur elsewhere.
#
##### ##### ##### ##### #####
>>>>>>> d14ca14a2707d349a9faad2450863e09bbea4e70


##### ##### ##### ##### #####
#
# vcfR class.
#
##### ##### ##### ##### #####


test_that("We can create an empty vcfR object",{
  myvcf <- methods::new("vcfR")
  
  expect_is(myvcf@meta, "character")
  expect_equal(length(myvcf@meta), 0)
  
  expect_is(myvcf@fix, "matrix")
  expect_equal(nrow(myvcf@fix), 0)
  expect_equal(ncol(myvcf@fix), 8)
  
  expect_is(myvcf@gt, "matrix")
  expect_equal(nrow(myvcf@gt), 0)
  expect_equal(ncol(myvcf@gt), 0)
  
})


test_that("vcfR subsetters works",{
<<<<<<< HEAD
  data("vcfR_example")
  
=======
  data(vcfR_example)
>>>>>>> d14ca14a2707d349a9faad2450863e09bbea4e70
  # Rows
  vcf2 <- vcf[1:10,]
  expect_equal(nrow(vcf2@fix), 10)
  expect_equal(nrow(vcf2@gt), 10)
  
  # Columns
  vcf2 <- vcf[,1:4]
  expect_equal(ncol(vcf2@fix), 8)
  expect_equal(ncol(vcf2@gt), 4)
})


<<<<<<< HEAD
=======
##### ##### ##### ##### #####
# EOF.
>>>>>>> d14ca14a2707d349a9faad2450863e09bbea4e70
