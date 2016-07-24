

#library(testthat)
library(vcfR)
context("vcfRtidy functions")

data("vcfR_example")
#data("vcfR_test")

##### ##### ##### ##### #####
# extract_gt_tidy

test_that("extract_gt_tidy works for GT element",{
#  Z <- extract_gt_tidy(vcf)
  suppressMessages( 
#
    Z <- extract_gt_tidy(vcf, format_fields = c('GT'))
#    Z <- extract_gt_tidy(vcf, format_fields = c('GT'), format_types = TRUE )
#    Z <- extract_gt_tidy( vcf, format_fields = c('GT'), format_types = character(0) )
  )
  expect_is(Z, 'tbl_df')
  expect_equal(names(Z)[1], 'Key')
  expect_equal(names(Z)[2], 'Indiv')
  expect_equal(names(Z)[3], 'gt_GT')
  expect_equal(names(Z)[4], 'gt_GT_alleles')
})


test_that("extract_gt_tidy works for all elements",{
  suppressMessages( Z <- extract_gt_tidy(vcf) )
  expect_is(Z, 'tbl_df')
  
})


##### ##### ##### ##### #####
# vcfR2tidy

test_that("vcfR2tidy works",{
  Z <- vcfR2tidy(vcf, info_only = FALSE)
  
  expect_is(Z, 'list')
  expect_equal( length(Z), 3 )
  
  expect_is(Z[['fix']], "tbl_df")
  expect_is(Z[['fix']], "tbl")
  expect_is(Z[['fix']], "data.frame")
  
})


##### ##### ##### ##### #####
# extract_info_tidy


test_that("extract_info_tidy works",{
  Z <- extract_info_tidy(vcf, info_fields = c("AC", "AN", "MQ"), info_types = c(AN = "i", MQ = "n"))

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
  
#  expect_equal( length(Z), 3 )

})


##### ##### ##### ##### #####
# vcf_field_names


test_that("vcf_field_names works",{
  Z <- vcf_field_names(vcf, tag = "INFO")
  
  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
  
  Z <- vcf_field_names(vcf, tag = "FORMAT")

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
})


