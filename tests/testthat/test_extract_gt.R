# extractgt devel

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("extract.gt functions")

#
#
# extract.gt tests ----
#
#


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


#
#
# extract.gt missing data. ----
#
#


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

#
#
# extract.gt return.alleles tests ----
#
#


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


# determine_ploidy tests ----

test_that(".determine_ploidy works on mixed ploid data",{
  # Fabricate data set.
  data(vcfR_test)
  gt4 <- matrix(nrow=5, ncol=2)
  colnames(gt4) <- c("NA00004", "NA00005")
  gt4[1,] <- c("1|0|0|1:48:8:51,51,51,51", "1|0|0|1:48:8:51,51,51,51")
  gt4[2,] <- c("1/0/0/1:48:8:51,51,51,51", "1/0/0/1:48:8:51,51,51,51")
  gt4[3,] <- c("1/2/2/1:48:8:51,51,51,51", "1/2/2/1:48:8:51,51,51,51")
  gt4[4,] <- c("0/0/1/2:48:8:51,51,51,51", "1/2/0/1:48:8:51,51,51,51")
  gt4[5,] <- c("1/0/0/1:48:8", "1|0|0|1:48:8")
  vcfR_test@gt <- cbind(vcfR_test@gt, gt4)
  is.na(vcfR_test@gt[1,3]) <- TRUE
  gt <- extract.gt(vcfR_test)
  
  my_ploidies <- .determine_ploidy(gt)
  expect_equal(my_ploidies, c(2, 2, 2, 4, 4))
})


#
#
# extract.haps tests ----
#
#


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

test_that("extract_haps works on haploid data",{
  data(vcfR_test)
  
  # Haploidize test data
  my_non_gt <- extract.gt(vcfR_test, extract = FALSE)
  my_gt <- extract.gt(vcfR_test)
  my_gt <- unlist(lapply(strsplit(my_gt, split = "[/|]"), function(x){x[1]}))
  # https://stackoverflow.com/a/35589023
  my_non_gt <- paste(my_gt, my_non_gt, sep = ":")
  dim(my_non_gt) <- dim(vcfR_test@gt[,-1])
  vcfR_test@gt[,-1] <- my_non_gt
  
  my_alleles <- extract.haps(vcfR_test)
  expect_equal(grep("[ACGT]", my_alleles, invert = TRUE), integer(0))
})


test_that(".extract_haps2 works on mixed ploid data",{
  data(vcfR_test)
  gt4 <- matrix(nrow=5, ncol=2)
  colnames(gt4) <- c("NA00004", "NA00005")
  gt4[1,] <- c("1|0|0|1:48:8:51,51,51,51", "1|0|0|1:48:8:51,51,51,51")
  gt4[2,] <- c("1/0/0/1:48:8:51,51,51,51", "1/0/0/1:48:8:51,51,51,51")
  gt4[3,] <- c("1/2/2/1:48:8:51,51,51,51", "1/2/2/1:48:8:51,51,51,51")
  gt4[4,] <- c("0/0/1/2:48:8:51,51,51,51", "1/2/0/1:48:8:51,51,51,51")
  gt4[5,] <- c("1/0/0/1:48:8", "1|0|0|1:48:8")
  vcfR_test@gt <- cbind(vcfR_test@gt, gt4)
  is.na(vcfR_test@gt[1,3]) <- TRUE
  vcfR_test@gt[2,3] <- "./."
  #  vcfR_test@fix[4,'ALT'] <- "A"
  gt <- extract.gt(vcfR_test)
  #    
  # #   myHaps <- .extract_haps(getREF(vcfR_test), getALT(vcfR_test), gt, 0, 1)
  # #   myHaps <- 
  #
  
  myStart <- 1
  myEnd <- 3
#  
  my_haps <- .extract_haps2(getREF(vcfR_test)[myStart:myEnd], 
                            getALT(vcfR_test)[myStart:myEnd], 
                            gt[myStart:myEnd, , drop = FALSE], 0, 0)
  
  expect_true(is.na(my_haps[1,3]))
  expect_true(is.na(my_haps[1,4]))
  
#  first_variant <- c(NA00001_0 = "G", NA00001_1 = "G", NA00002_0 = NA, NA00002_1 = NA, 
#                     NA00003_0 = "A", NA00003_1 = "A", NA00004_0 = "A", NA00004_1 = "G", 
#                     NA00004_2 = "G", NA00004_3 = "A", NA00005_0 = "A", NA00005_1 = "G", 
#                     NA00005_2 = "G", NA00005_3 = "A")

  
#  expect_equal(my_haps[1,], first_variant)
  # Verbose output.
  # 
  myStart <- 1
  myEnd <- 1
#  .extract_haps2(getREF(vcfR_test)[myStart:myEnd], 
#                 getALT(vcfR_test)[myStart:myEnd], 
#                 gt[myStart:myEnd, , drop = FALSE], 0, 1)
})


#
#
# extract.indels tests ----
#
#


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


#
#
# gt2alleles tests ----
#
#


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



#
#
# extract_gt_to_CM ----
#
#



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



#
#
# extract.info tests ----
#
#


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


####
# EOF.
