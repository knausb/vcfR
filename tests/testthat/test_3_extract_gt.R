# extractgt devel

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("extract.gt functions")

#
data(vcfR_example)

#vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
#seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
#gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

#vcf <- read.vcfR(vcf_file, verbose = FALSE)
#dna <- ape::read.dna(seq_file, format = "fasta")
#gff <- read.table(gff_file, sep="\t")

chrom <- create.chromR(name="Chrom", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_DP = 1e3, max_DP = 2e3)

##### ##### ##### ##### #####


gt <- extract.gt(chrom, element="GT", as.numeric=FALSE)
gt2 <- extract.gt(chrom, element="GT", as.numeric=FALSE, mask = TRUE)
gt3 <- extract.gt(chrom, element="GT", mask = c(TRUE, FALSE)) # Recycled vector


pl <- extract.gt(chrom, element="PL", as.numeric=FALSE)
gq <- extract.gt(chrom, element="GQ", as.numeric=TRUE)


##### ##### ##### ##### #####
#
# extract.gt tests
#
##### ##### ##### ##### #####


test_that("gt, pl ad gq are matrices",{
  expect_is(gt, "matrix")
  expect_is(gt2, "matrix")
  expect_is(pl, "matrix")
  expect_is(gq, "matrix")
})

test_that("gq is numeric",{
  expect_is(gq, "matrix")
  expect_equal(is.numeric(gq), TRUE)
})


test_that("extract.gt mask=TRUE works", {
  expect_equal(nrow(gt), nrow(vcf@gt))
  expect_equal(nrow(gt2), sum(chrom@var.info$mask))
  expect_equal(nrow(gt3), nrow(gt[c(TRUE, FALSE),]))
})


test_that("extract.gt extract parameter works",{
  gt <- extract.gt(chrom, element="GT", extract=TRUE)
  expect_is(gq, "matrix")
  expect_equal(is.numeric(gq), TRUE)
})


##### ##### ##### ##### #####
#
# extract.gt return.alleles tests
#
##### ##### ##### ##### #####

test_that("extract.gt return.alleles works",{
  gt <- extract.gt(chrom, element="GT", return.alleles = TRUE, allele.sep="|")
  expect_is(gt, "matrix")
  
})


test_that("extract.gt return.alleles works for multiallelic variants",{
#  vcf2 <- vcf[nchar(vcf@fix[,'ALT']) > 1,]
  vcf2 <- vcf[grep(",", vcf@fix[,'ALT']),]
  gt <- extract.gt(vcf2, element="GT", return.alleles = TRUE, allele.sep="|")
  
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

data(vcfR_example)
gt <- extract.gt(vcf, element="GT", extract=TRUE)

test_that("extract_haps compiled code works",{
  is.na(gt[1:5,1]) <- TRUE
  gt[1:5,2] <- ".|."
  gt[4,2] <- "0|."
  gt[5,2] <- ".|1"
    
  haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', vcf@fix[,'REF'], vcf@fix[,'ALT'], gt, '|', 0)
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
  vcf@gt[1,3] <- ".|.:12,0:12:39:0,39,585"
  vcf@gt[2,3] <- "0|.:12,0:12:39:0,39,585"
  vcf@gt[3,3] <- ".|1:12,0:12:39:0,39,585"
  
  haps <- extract.haps(vcf, gt.split="|", verbose = FALSE)
  haps[1:6,1:8]
  
  expect_true( is.na(haps[1,3]) )
  expect_true( is.na(haps[1,4]) )
  
  expect_true( !is.na(haps[2,3]) )
  expect_true( is.na(haps[2,4]) )
  
  expect_true( is.na(haps[3,3]) )
  expect_true( !is.na(haps[3,4]) )
})


chrom <- create.chromR(name="Chrom", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_DP = 1e3, max_DP = 2e3)

test_that("extract_haps R code works on chromR objects",{
  haps <- extract.haps(chrom, gt.split="|", verbose = FALSE)
  expect_is(haps, "matrix")
  expect_equal(ncol(haps), 2 * ncol(gt))
  expect_equal(nrow(haps), nrow(gt))
})




##### ##### ##### ##### #####
#
# extract.indels tests
#
##### ##### ##### ##### #####


test_that("extract.indels works",{
  indels <- extract.indels(vcf, return.indels=TRUE)
  expect_equal(nrow(indels@fix), 328)
})


##### ##### ##### ##### #####
#
# gt2alleles tests
#
##### ##### ##### ##### #####


data(vcfR_example)

vcf <- vcf[1:4, 1:3]
vcf@gt[2,3] <- ".|0:12,0:12:39:0,39,585"
vcf@gt[3,3] <- "0|.:12,0:12:39:0,39,585"
vcf@gt[4,3] <- ".|.:12,0:12:39:0,39,585"

test_that("gt2alleles works",{
#  cbind(vcf@fix[,4:5], vcf@gt)
  gt <- extract.gt(vcf, return.alleles=TRUE, allele.sep = "|")
  expect_equal( as.character(gt[1,2]), "T|T")
  expect_equal( as.character(gt[2,2]), ".|C")
  expect_equal( as.character(gt[3,2]), "A|.")
  expect_equal( as.character(gt[4,2]), ".|.")
})



##### ##### ##### ##### #####
#
#
#
##### ##### ##### ##### #####



test_that("extract_gt_to_CM2 compiled code works",{
  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'GT', '|', 0, 1 )
#  head(gt)
  # Return alleles
  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'GT', '|', 1, 1 )
#  head(gt)
  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'DP', '|', 0, 1 )
#  head(gt)
  
  # Manage NA_STRING
#  grep("\\.", gt[,1], value = TRUE)
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