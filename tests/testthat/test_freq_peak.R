



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
  expect_equal(length(myPeaks1), 3)
  expect_equal( sum(myPeaks1$peaks >= 0), length(myPeaks1$peaks))
  expect_equal( sum(myPeaks1$peaks <= 1), length(myPeaks1$peaks))
  
})


test_that("freq_peak works, zero variants",{
  data("vcfR_test")
  vcfR_test <- vcfR_test[0,]
  dp <- extract.gt(vcfR_test, element = "DP", as.numeric = TRUE)
  
  myPeaks1 <- freq_peak(dp, pos = getPOS(vcfR_test))
#  myPeaks1 <- freq_peak(dp[0,], pos = getPOS(vcfR_test))
  
  expect_is(myPeaks1, 'list')
  expect_equal(nrow(myPeaks1$wins), 0)
  expect_equal(nrow(myPeaks1$peaks), 0)
})


test_that("freq_peak works, one variant",{
  data("vcfR_test")
  vcfR_test <- vcfR_test[1,]
  dp <- extract.gt(vcfR_test, element = "DP", as.numeric = TRUE)
  
  myPeaks1 <- freq_peak(dp, pos = getPOS(vcfR_test))
  #  myPeaks1 <- freq_peak(dp[0,], pos = getPOS(vcfR_test))
  
  expect_is(myPeaks1, 'list')
})

