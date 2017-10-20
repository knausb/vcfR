

#library(testthat)
library(vcfR)
context("rePOS functions")


##### ##### ##### ##### #####


test_that("rePOS works",{
  # Create some VCF data.
  data(vcfR_example)
  vcf1 <-vcf[1:500,]
  vcf2 <-vcf[500:1500,]
  vcf3 <- vcf[1500:2533]
  vcf1@fix[,'CHROM'] <- 'chrom1'
  vcf2@fix[,'CHROM'] <- 'chrom2'
  vcf3@fix[,'CHROM'] <- 'chrom3'
  vcf2@fix[,'POS'] <- as.character(getPOS(vcf2) - 21900)
  vcf3@fix[,'POS'] <- as.character(getPOS(vcf3) - 67900)
  vcf <- rbind2(vcf1, vcf2)
  vcf <- rbind2(vcf, vcf3)
  rm(vcf1, vcf2, vcf3)

  # Create lens
  lens <- data.frame(matrix(nrow=3, ncol=2))
  lens[1,1] <- 'chrom1'
  lens[2,1] <- 'chrom2'
  lens[3,1] <- 'chrom3'
  lens[1,2] <- 22000
  lens[2,2] <- 47000
  lens[3,2] <- 32089
  
  newPOS <- rePOS(vcf, lens)
  
  expect_equal(length(newPOS), nrow(vcf))
  expect_true(newPOS[length(newPOS)] > getPOS(vcf)[nrow(vcf)])
})



