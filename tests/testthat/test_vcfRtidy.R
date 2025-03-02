

#
library(testthat)
library(vcfR)
context("vcfRtidy functions")

data("vcfR_example")
#data("vcfR_test")


##### ##### ##### ##### #####
# vcf_field_names


test_that("vcf_field_names works",{
  data("vcfR_example")
  Z <- vcf_field_names(vcf, tag = "INFO")
  
  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
  
  Z <- vcf_field_names(vcf, tag = "FORMAT")

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
})


test_that("vcf_field_names works, comma in quotes not parsed",{
   data("vcfR_test")
   myMeta <- vcfR_test@meta
   vcfR_test@meta <- c(myMeta[1:12], '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">', myMeta[13:18])

   Z <- vcf_field_names(vcfR_test, tag = "INFO")
   expect_is(Z, "tbl_df")
   expect_is(Z, "tbl")
   expect_is(Z, "data.frame")
  
   vcfR_test@meta[13] <- '##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">'
   Z <- vcf_field_names(vcfR_test)
   
   expect_is(Z, "tbl_df")
   expect_is(Z, "tbl")
   expect_is(Z, "data.frame")
})


test_that("vcf_field_names works, zero INFO records in meta",{
  data("vcfR_test")
  ## cause error ##
  # remove all INFO lines in @meta object
  INFO.meta.lines <- grepl("^##INFO", vcfR_test@meta);
  vcfR_test@meta <- vcfR_test@meta[!INFO.meta.lines];
  
  #debug(vcfR2tidy)
  #debug(vcf_field_names)
  tidyVCF <- vcfR2tidy(vcfR_test)
  
  expect_is(tidyVCF, "list")
  expect_identical(names(tidyVCF), c("fix", "gt", "meta"))
  expect_identical(
    names(tidyVCF$fix),
    c("ChromKey", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "VariantKey")
    #c("ChromKey", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  )
})


test_that("vcf_field_names works, equal sign (=) in key.",{
  library(vcfR)
  data("vcfR_test")
  vcfR_test@meta[ grep("##FORMAT=<ID=DP", vcfR_test@meta) ] <- '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'  
  my_fnames <- vcf_field_names(vcfR_test, tag = "FORMAT")
  my_fnames
  expect_identical(
    my_fnames$Description[3],
    "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"
  )
})


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
  #expect_equal(names(Z)[1], 'Key')
  expect_equal(names(Z)[1], 'VariantKey')
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
  data("vcfR_test")
  Z <- vcfR2tidy(vcfR_test, info_only = FALSE)
  
  expect_is(Z, 'list')
  expect_equal( length(Z), 3 )
  
  expect_is(Z[['fix']], "tbl_df")
  expect_is(Z[['fix']], "tbl")
  expect_is(Z[['fix']], "data.frame")
  
})

# test_that("vcfR2tidy works, ID=REF",{
#    data("vcfR_test")
#    myMeta <- vcfR_test@meta
#    vcfR_test@meta <- c(myMeta[1:12], '##INFO=<ID=REF,Number=0,Type=Flag,Description="Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8">', myMeta[13:18])  
#    Z <- vcfR2tidy(vcf, info_only = FALSE, verbose = FALSE)
#    Z$meta$ID
# })

##### ##### ##### ##### #####
# extract_info_tidy


test_that("extract_info_tidy works",{
  Z <- extract_info_tidy(vcf, info_fields = c("AC", "AN", "MQ"), info_types = c(AN = "i", MQ = "n"))

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
  
#  expect_equal( length(Z), 3 )

})

test_that("extract_info_tidy works with Flags",{
  data(vcfR_test)
  Z <- extract_info_tidy(vcfR_test, info_types = TRUE)
  expect_is(Z$DB, "logical")
  expect_is(Z$H2, "logical")
})




