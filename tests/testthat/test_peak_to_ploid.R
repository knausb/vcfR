

#library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)

#
context("peak_to_ploid")


test_that("compiled vcfR_vcf_stats_gz works",{
  data(vcfR_example)
  gt <- extract.gt(vcf)
  # Censor non-heterozygous positions.
  hets <- is_het(gt)
  is.na(vcf@gt[,-1][!hets]) <- TRUE
  # Extract allele depths.
  ad <- extract.gt(vcf, element = "AD")
  ad1 <- masplit(ad, record = 1)
  ad2 <- masplit(ad, record = 2)
  freq1 <- ad1/(ad1+ad2)
  freq2 <- ad2/(ad1+ad2)
  myPeaks1 <- freq_peak(freq1, getPOS(vcf))
  # Censor windows with fewer than 20 heterozygous positions
  is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
  # Convert peaks to ploidy call
  myCalls <- peak_to_ploid(myPeaks1)
  
  expect_equal(class(myCalls), "list")
  expect_equal(length(myCalls), 2)
  expect_true(max(myCalls$dfe, na.rm = TRUE) < 1)
  expect_true(min(myCalls$dfe, na.rm = TRUE) > -1)
  expect_equal(ncol(myCalls$calls), ncol(myPeaks1$peaks))
  expect_equal(nrow(myCalls$calls), nrow(myPeaks1$peaks))
})


