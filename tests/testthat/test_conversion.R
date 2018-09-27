

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

test_that("vcfR2genind works, return.alleles = TRUE",{
  data(vcfR_test)
  suppressMessages(my_genind <- vcfR2genind(vcfR_test, return.alleles = TRUE))

  expect_is(my_genind, "genind")
  my_alleles <- unlist(lapply(strsplit(colnames(my_genind@tab), ".", fixed = TRUE), function(x){x[2]}))
  expect_equal(sum(grepl("[A|C|G|T]", my_alleles)), length(my_genind@loc.fac))
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

test_that("vcfR2loci works, return.alleles = TRUE",{
  data(vcfR_test)
  myLoci <- vcfR2loci(vcfR_test, return.alleles = TRUE)
  expect_is(myLoci, "loci")
  my_alleles <- as.character(unlist(myLoci))
  my_alleles <- unlist(strsplit(my_alleles, "[|/]"))
  expect_equal(sum(grepl("[A|C|G|T]", my_alleles)), length(my_alleles))
})

#debug(vcfR2DNAbin)

##### ##### ##### ##### #####
# EOF.
