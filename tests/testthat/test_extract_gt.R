# extractgt devel

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("extract.gt functions")

##### ##### ##### ##### #####
#
# extract.gt tests
#
##### ##### ##### ##### #####


test_that("gt is a matrix",{
  data(chromR_example)
  gt <- extract.gt(chrom, element="GT", as.numeric=FALSE)
  
  expect_is(gt, "matrix")
})

test_that("gq is numeric",{
  data(chromR_example)
  gq <- extract.gt(chrom, element="GQ", as.numeric=TRUE)
  
  expect_is(gq, "matrix")
  expect_equal(is.numeric(gq), TRUE)
})


test_that("extract.gt mask=TRUE works", {
  data(chromR_example)
  gt <- extract.gt(chrom, element="GT", as.numeric=FALSE, mask = TRUE)
  
  expect_equal(nrow(gt), sum(chrom@var.info$mask))
})


test_that("extract.gt extract parameter works",{
  data(chromR_example)
  gt <- extract.gt(chrom, element="GT", extract=FALSE)
  expect_is(gt, "matrix")
  expect_equal(class(gt[1,1]), "character")
  expect_true(nchar(gt[1,1]) > 4)
})


##### ##### ##### ##### #####
#
# extract.gt missing data.
#
##### ##### ##### ##### #####


test_that("extract.gt converts missing GT to NA",{
  data("vcfR_test")
  
  vcfR_test@gt[1,2] <- "./.:48:1:51,51"
  gt <- extract.gt(vcfR_test)
  expect_true( is.na(gt[1,1]) )
  
  vcfR_test@gt[1,2] <- ".|.:48:1:51,51"
  gt <- extract.gt(vcfR_test)
  expect_true( is.na(gt[1,1]) )
  
  vcfR_test@gt[1,2] <- ".|0:48:1:51,51"
  gt <- extract.gt(vcfR_test)
  expect_false( is.na(gt[1,1]) )
})


test_that("extract.gt convertNA = FALSE works",{
  data("vcfR_test")
  
  vcfR_test@gt[1,2] <- "./.:48:1:51,51"
  gt <- extract.gt(vcfR_test, convertNA = FALSE)
  expect_false( is.na( gt[1,1] ) )
  
  vcfR_test@gt[1,2] <- ".|.:48:1:51,51"
  gt <- extract.gt(vcfR_test, convertNA = FALSE)
  expect_false( is.na(gt[1,1]) )

})

##### ##### ##### ##### #####
#
# extract.gt return.alleles tests
#
##### ##### ##### ##### #####


test_that("extract.gt return.alleles works #1",{
  data(vcfR_test)
  gt <- extract.gt(vcfR_test, element="GT")
  gt <- extract.gt(vcfR_test, element="GT", return.alleles=TRUE)
  expect_is(gt, "matrix")
})


test_that("extract.gt return.alleles works #2",{
  data(vcfR_example)
#  gt <- extract.gt(vcf[438,1:2], element="GT", return.alleles=TRUE)
  gt <- extract.gt(vcf[439,1:2], element="GT", return.alleles=TRUE)

  expect_is(gt, "matrix")
})

# This is my stack overflow.
test_that("extract.gt return.alleles works",{
  data(chromR_example)
  gt <- extract.gt(chrom, 
                   element="GT", 
                   return.alleles = TRUE
#                   allele.sep="|"
                   )
  expect_is(gt, "matrix")
})


test_that("extract.gt return.alleles works for multiallelic variants",{
  data(vcfR_example)
#  vcf2 <- vcf[nchar(vcf@fix[,'ALT']) > 1,]
  vcf2 <- vcf[grep(",", vcf@fix[,'ALT']),]
  gt <- extract.gt(vcf2, 
                   element="GT", 
                   return.alleles = TRUE
#                   allele.sep="|"
                   )
  
  # Locus 1: 0,1,2 = T,A,G
  # 0|0
  expect_equal( as.character(gt[1,'BL2009P4_us23']), "T|T")
  # 0|1
  expect_equal( as.character(gt[1,'DDR7602']), "T|A")
  # 0|2
  expect_equal( as.character(gt[1,'blue13']), "T|G")
  
  # Locus 6: 0,1,2 = G,GTCTAATAGAGGCTCGAACTC,GTCTAATAGAGGCTCGAGCTC
  # 0|2
  expect_equal( as.character(gt[6,'BL2009P4_us23']), "G|GTCTAATAGAGGCTCGAGCTC")
  # 1|2
  expect_equal( as.character(gt[6,'DDR7602']), "GTCTAATAGAGGCTCGAACTC|GTCTAATAGAGGCTCGAGCTC" )

})



##### ##### ##### ##### #####
#
# extract.haps tests
#
##### ##### ##### ##### #####


test_that("extract_haps compiled code works",{
  data(vcfR_example)
  gt <- extract.gt(vcf, element="GT", extract=TRUE)

  is.na(gt[1:5,1]) <- TRUE
  gt[1:5,2] <- ".|."
  gt[4,2] <- "0|."
  gt[5,2] <- ".|1"
    
  haps <- .extract_haps(vcf@fix[,'REF'], vcf@fix[,'ALT'], gt, 0, 0)
  
  expect_is(haps, "matrix")
  expect_true( is.na(haps[1,1]) )
  expect_true( is.na(haps[1,2]) )

  expect_true( is.na(haps[1,3]) )
  expect_true( is.na(haps[1,4]) )
  expect_true( !is.na(haps[4,3]) )
  expect_true( is.na(haps[4,4]) )
  expect_true( is.na(haps[5,3]) )
  expect_true( !is.na(haps[5,4]) )
  
  expect_equal(ncol(haps), 2 * ncol(gt))
  expect_equal(nrow(haps), nrow(gt))
})


test_that("extract_haps R code works on vcfR objects",{
  data(vcfR_example)
  
  vcf@gt[1,3] <- ".|.:12,0:12:39:0,39,585"
  vcf@gt[2,3] <- "0|.:12,0:12:39:0,39,585"
  vcf@gt[3,3] <- ".|1:12,0:12:39:0,39,585"
  
  haps <- extract.haps(vcf, verbose = FALSE)
  haps[1:6,1:8]
  
  expect_true( is.na(haps[1,3]) )
  expect_true( is.na(haps[1,4]) )
  
  expect_true( !is.na(haps[2,3]) )
  expect_true( is.na(haps[2,4]) )
  
  expect_true( is.na(haps[3,3]) )
  expect_true( !is.na(haps[3,4]) )
})


test_that("extract_haps R code works on chromR objects",{
  data(chromR_example)
  haps <- extract.haps(chrom, verbose = FALSE)
  expect_is(haps, "matrix")
  expect_equal(ncol(haps), 2 * (ncol(chrom@vcf@gt) - 1) )
  expect_equal(nrow(haps), nrow(chrom@vcf@gt))
})


test_that("extract_haps unphased_as_NA works",{
  data(vcfR_test)
  haps <- extract.haps(vcfR_test, unphased_as_NA = FALSE, verbose = FALSE)
  expect_equal(sum(is.na(haps)), 0)
  haps <- extract.haps(vcfR_test, unphased_as_NA = TRUE, verbose = FALSE)
  expect_equal(sum(is.na(haps)), 14)
})


##### ##### ##### ##### #####
#
# extract.indels tests
#
##### ##### ##### ##### #####


test_that("extract.indels works",{
  data(vcfR_example)
  indels <- extract.indels(vcf, return.indels=TRUE)
  expect_equal(nrow(indels@fix), 328)
})

test_that("extract.indels works, <NON_REF>",{
  data(vcfR_test)
  indels <- extract.indels(vcfR_test, return.indels=FALSE)
  
  data(vcfR_test)
  vcfR_test@fix[1,'ALT'] <- "<NON_REF>"
  indels2 <- extract.indels(vcfR_test, return.indels=FALSE)
  expect_equal(nrow(indels), nrow(indels2))
  expect_equal(sum(is.na(getPOS(indels2))), 0)
  
  data(vcfR_test)
  vcfR_test@fix[1,'ALT'] <- "A,<NON_REF>"
  indels2 <- extract.indels(vcfR_test, return.indels=FALSE)
  expect_equal(nrow(indels), nrow(indels2))
})


##### ##### ##### ##### #####
#
# gt2alleles tests
#
##### ##### ##### ##### #####


#data(vcfR_example)

#vcf <- vcf[1:4, 1:3]
#vcf@gt[2,3] <- ".|0:12,0:12:39:0,39,585"
#vcf@gt[3,3] <- "0|.:12,0:12:39:0,39,585"
#vcf@gt[4,3] <- ".|.:12,0:12:39:0,39,585"

test_that("gt2alleles works",{
  data("vcfR_test")
  vcfR_test@gt
  gt <- extract.gt(vcfR_test, return.alleles=TRUE)
#  gt
  
  expect_equal( as.character(gt[1,1]), "G|G")
  expect_equal( as.character(gt[1,3]), "A/A")
  expect_equal( as.character(gt[5,1]), "GTC/G")
  
#  cbind(vcf@fix[,4:5], vcf@gt)
#  gt <- extract.gt(vcf, return.alleles=TRUE, allele.sep = "|")
#  expect_equal( as.character(gt[1,2]), "T|T")
#  expect_equal( as.character(gt[2,2]), ".|C")
#  expect_equal( as.character(gt[3,2]), "A|.")
#  expect_equal( as.character(gt[4,2]), ".|.")
})



##### ##### ##### ##### #####
#
#
#
##### ##### ##### ##### #####



test_that("extract_gt_to_CM compiled code works",{
  data(vcfR_example)
#  gt <- .extract_GT_to_CM(fix = vcf@fix, gt = vcf@gt,
#                          element = 'GT',
#                          alleles = 0, extract = 1, convertNA = 1 )

#  gt <- .extract_GT_to_CM(fix = vcf@fix, gt = vcf@gt, 
#                          element = 'GT',
#                          alleles = 1, extract = 1, convertNA = 1 )

  gt <- .extract_GT_to_CM(fix = vcf@fix, gt = vcf@gt,
                          element = 'DP',
                          alleles = 0, extract = 1, convertNA = 1 )

  expect_equal(length( grep("\\.", gt) ), 0)
})


##### ##### ##### ##### #####


##### ##### ##### ##### #####
#
# extract.info tests
#
##### ##### ##### ##### #####


test_that("extract.info handles missing elements",{
  data(vcfR_test)
  
  # Element is present in some, but not all, samples.
  info <- extract.info(vcfR_test, element = "AF")
  expect_equal(length(info), nrow(vcfR_test))
  expect_equal( length(grep("AF", vcfR_test@fix[,'INFO'])), sum( !is.na(info) ) )
  
  # Element does not exist.
  info <- extract.info(vcfR_test, element = "XX")
  expect_equal( length(grep("XX", vcfR_test@fix[,'INFO'])), sum( !is.na(info) ) )
})


##### ##### ##### ##### #####
# EOF.
