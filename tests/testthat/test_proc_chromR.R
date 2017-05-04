


# detach(package:vcfR, unload=T)
#
library(testthat)
library(vcfR)
context("proc_chromR")



##### ##### ##### ##### #####

test_that("vcfR_gt_to_popsum",{
  data("vcfR_test")
  gt <- extract.gt(vcfR_test)
  var_info <- as.data.frame(vcfR_test@fix[,1:2])
  var_info$mask <- TRUE
  tmp <- .Call('vcfR_gt_to_popsum', PACKAGE = 'vcfR', var_info = var_info, gt = gt)
  
  expect_equal(class(tmp), "data.frame")
  expect_equal(nrow(tmp), nrow(gt))
  expect_equal(ncol(tmp), 7)
  expect_equal(tmp$Allele_counts, structure(c(3L, 4L, 1L, 5L, 2L), .Label = c("0,2,4", "2,3,1", "3,3", "5,1", "6"), class = "factor"))
})


test_that("vcfR_gt_to_popsum, one variant",{
  data("vcfR_test")
  vcfR_test <- vcfR_test[1,]
  gt <- extract.gt(vcfR_test)
  var_info <- as.data.frame(vcfR_test@fix[,1:2, drop = FALSE])
  var_info$mask <- TRUE
  tmp <- .Call('vcfR_gt_to_popsum', PACKAGE = 'vcfR', var_info = var_info, gt = gt)
  
  expect_equal(class(tmp), "data.frame")
  expect_equal(nrow(tmp), 1)
})


test_that("vcfR_gt_to_popsum, counts missing alleles",{
  data("vcfR_test")
  gt <- extract.gt(vcfR_test)
  gt[1,] <- "2/2"
  gt[2,] <- "0/2"
  var_info <- as.data.frame(vcfR_test@fix[,1:2, drop = FALSE])
  var_info$mask <- TRUE
  tmp <- .Call('vcfR_gt_to_popsum', PACKAGE = 'vcfR', var_info = var_info, gt = gt)
  
  expect_identical(as.character(tmp$Allele_counts[1]), "0,0,6")
  expect_equal(as.character(tmp$Allele_counts[2]), "3,0,3")
})

##### ##### ##### ##### #####


