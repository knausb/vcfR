

#
library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("conversion functions")

##### ##### ##### ##### #####

test_that("vcfR2genind works",{
  data(vcfR_test)
  suppressMessages(my_genind <- vcfR2genind(vcfR_test))

  expect_is(my_genind, "genind")  
})


##### ##### ##### ##### #####

test_that("vcfR2genlight works",{
  suppressMessages(library(adegenet))
  library(parallel)

  data(vcfR_test)
  vcfR_test <- vcfR_test[is.biallelic(vcfR_test),]
  my_genlight <- vcfR2genlight(vcfR_test)
  expect_is(my_genlight, "genlight")
})


##### ##### ##### ##### #####

test_that("vcfR2loci works",{
  data(vcfR_test)
  myLoci <- vcfR2loci(vcfR_test)
  expect_is(myLoci, "loci")
})

#debug(vcfR2DNAbin)

##### ##### ##### ##### #####
# EOF.
