



# detach(package:vcfR, unload=T)
#
library(testthat)
library(vcfR)

context("freq_peak")


#data("vcfR_example")

test_that("freq_peak works",{
  data(vcfR_example)
  gt <- extract.gt(vcf)
  hets <- is_het(gt)
  # Censor non-heterozygous positions.
  is.na(vcf@gt[,-1][!hets]) <- TRUE
  # Extract allele depths.
  ad <- extract.gt(vcf, element = "AD")
  ad1 <- masplit(ad, record = 1)
  ad2 <- masplit(ad, record = 2)
  freq1 <- ad1/(ad1+ad2)
  freq2 <- ad2/(ad1+ad2)
  
  
  myPeaks1 <- freq_peak(freq1, getPOS(vcf))
  
  expect_is(myPeaks1, 'list')
  expect_equal(length(myPeaks1), 2)
  expect_equal( sum(myPeaks1$peaks >= 0), length(myPeaks1$peaks))
  expect_equal( sum(myPeaks1$peaks <= 1), length(myPeaks1$peaks))
  
})
