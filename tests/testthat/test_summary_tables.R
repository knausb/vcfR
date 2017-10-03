


#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("summary_tables functions")

##### ##### ##### ##### #####

test_that("write.var.info works",{
#  data(vcfR_test)
  #myChrom <- create.chromR(vcfR_test, verbose = FALSE)
  data("chromR_example")
  
  setwd(tempdir())
  
  write.var.info(chrom, file = "test_var_info.csv")
  expect_true(file.exists("test_var_info.csv"))
  unlink("test_var_info.csv")
})


test_that("write.var.info works, mask == TRUE",{
  #data(vcfR_test)
  #myChrom <- create.chromR(vcfR_test, verbose = FALSE)
  data("chromR_example")
    
  setwd(tempdir())
  
  write.var.info(chrom, file = "test_var_info.csv", mask = TRUE)
  expect_true(file.exists("test_var_info.csv"))
  unlink("test_var_info.csv")
})


##### ##### ##### ##### #####

test_that("write.win.info works",{
#  data(vcfR_test)
#  myChrom <- create.chromR(vcfR_test, verbose = FALSE)
  data("chromR_example")
  myChrom <- proc.chromR(chrom, verbose = FALSE)
  
  setwd(tempdir())
  write.win.info(myChrom, file = "test_win_info.csv")
  expect_true(file.exists("test_win_info.csv"))
  unlink("test_win_info.csv")
})


