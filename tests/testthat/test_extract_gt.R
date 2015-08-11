# extractgt devel

#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("extract_gt functions")

#data(vcfR_example)

vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

vcf <- read.vcf(vcf_file, verbose = FALSE)
dna <- ape::read.dna(seq_file, format = "fasta")
gff <- read.table(gff_file, sep="\t")

chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)

#####


chrom <- create_chrom('Supercontig_1.100', seq=dna, vcf=vcf, ann=gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)
pinf_mt <- masker(pinf_mt, min_QUAL = 990, min_DP = 6000, max_DP = 8000, min_MQ = 40, max_MQ = 100)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)

gt <- extract.gt(pinf_mt, element="GT", as.numeric=FALSE)
gt2 <- extract.gt(pinf_mt, element="GT", as.numeric=FALSE, mask = TRUE)
gt3 <- extract.gt(pinf_vcf, element="GT", mask = c(TRUE, FALSE))


pl <- extract.gt(pinf_mt, element="PL", as.numeric=FALSE)
gq <- extract.gt(pinf_mt, element="GQ", as.numeric=TRUE)


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


test_that("extract_gt mask=TRUE works", {
  expect_equal(nrow(gt), 371)
  expect_equal(nrow(gt2), 212)
  expect_equal(nrow(gt3), 186)
})


test_that("extract_indels works",{
  indels <- extract_indels(pinf_vcf, return_indels=TRUE)
  expect_equal(min(nchar(indels@fix$REF)), 2)
  expect_equal(min(nchar(indels@fix$ALT)), 2)
  
  indels <- extract_indels(pinf_vcf, return_indels=FALSE)
  expect_equal(max(nchar(indels@fix$REF)), 1)
  expect_equal(max(nchar(indels@fix$ALT)), 1)
})


test_that("extract_haps works",{
#  .Call('vcfR_extract_haps', PACKAGE = 'vcfR', ref, alt, gt, vebosity)
#  haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', pinf_vcf@fix$REF, pinf_vcf@fix$ALT, gt, 1)
  haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', pinf_vcf@fix$REF, pinf_vcf@fix$ALT, gt, '/', 1)
#  head(haps)
})


#head(gt)
#head(pl)
#head(gq)

#head(pinf_mt@vcf.gt)


#head(pinf_vcf)


#ncol(pinf_vcf@gt)

#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GT")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="PL")
#outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="DP")


#outm <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', pinf_vcf@gt, element="DP")

#outm <- extract.gt2(pinf_vcf, element="DP", as.numeric=F)
#outm <- extract.gt2(pinf_vcf, element="GQ", as.numeric=T)

#head(pinf_vcf@gt)
#head(outm)



#outm <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")

#head(outm)


#outm2 <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outm)
#head(outm2)

