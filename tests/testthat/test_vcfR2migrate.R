
#
library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("vcfR2migrate")

##### ##### ##### ##### #####

test_that("vcfR2migrate method N works",{
  data(vcfR_test)
  my_pop <- as.factor(paste("pop_", rep(c("A", "B", "C"), each = 6), sep = ""))
  
  setwd(tempdir())
  vcfR2migrate(vcf = vcf , pop = my_pop , in_pop = c("pop_A","pop_C"), out_file = "my2pop.txt", method = 'N')
  
  expect_true(file.exists("my2pop.txt"))
  unlink("my2pop.txt")
})


test_that("vcfR2migrate method H works",{
  data(vcfR_test)
  my_pop <- as.factor(paste("pop_", rep(c("A", "B", "C"), each = 6), sep = ""))
  
  setwd(tempdir())
  vcfR2migrate(vcf = vcf , pop = my_pop , in_pop = c("pop_A","pop_C"), out_file = "my2pop.txt", method = 'H')
  
  expect_true(file.exists("my2pop.txt"))
  unlink("my2pop.txt")
})


##### ##### ##### ##### #####
