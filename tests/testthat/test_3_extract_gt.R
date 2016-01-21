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
# extract.haps tests
#
##### ##### ##### ##### #####


test_that("extract_haps compiled code works",{
  is.na(gt[1:5,1]) <- TRUE
  haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', vcf@fix[,'REF'], vcf@fix[,'ALT'], gt, '|', 0)
  expect_is(haps, "matrix")
  expect_true( is.na(haps[1,1]) )
  expect_true( is.na(haps[1,2]) )
  expect_equal(ncol(haps), 2 * ncol(gt))
  expect_equal(nrow(haps), nrow(gt))
})


test_that("extract_haps R code works",{
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
  

#  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'GT', '/', 0, 0 )
#  
#  head(gt)
#  head(vcf@gt)
#  expect_equal(nchar(gt[1,1]), 17)
  
#  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'AD', '/', 0, 0 )
#  expect_equal(nchar(gt[1,1]), 17)
  
#  gt <- .Call( 'vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', vcf@fix, vcf@gt, 'PL', '/', 0, 0 )
#  expect_equal(nchar(gt[1,1]), 12)
  
#  expect_is(haps, "matrix")
#  expect_equal(ncol(haps), 2 * ncol(gt))
#  expect_equal(nrow(haps), nrow(gt))
})


##### ##### ##### ##### #####


##### ##### ##### ##### #####
#
# extract.info tests
#
##### ##### ##### ##### #####


test_that("extract.info works",{
  

})


##### ##### ##### ##### #####
# EOF.