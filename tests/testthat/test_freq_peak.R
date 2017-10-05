



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
  expect_equal(length(myPeaks1), 5)
  expect_equal( sum(myPeaks1$peaks >= 0), length(myPeaks1$peaks))
  expect_equal( sum(myPeaks1$peaks <= 1), length(myPeaks1$peaks))
  
  expect_equal( colSums(myPeaks1$counts), apply(freq1, MARGIN = 2, function(x){sum(!is.na(x))}) )
  
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
  
  myPeaks1 <- freq_peak(dp, pos = getPOS(vcfR_test), lhs = FALSE)
  #  myPeaks1 <- freq_peak(dp[0,], pos = getPOS(vcfR_test))
  
  expect_is(myPeaks1, 'list')
})


test_that("freq_peak works, one variant",{
  data("vcfR_test")
  dp <- extract.gt(vcfR_test, element = "DP", as.numeric = TRUE)
  
  set.seed(1)
  dp[1:nrow(dp), 1:ncol(dp)] <- as.matrix( jitter( rep( 0.5, times = length(dp) ) ) )
  dp[1:2,2] <- 1/4
  
  myPeaks1 <- freq_peak(dp, pos = getPOS(vcfR_test), lhs = FALSE)

#  myPeaks1$wins[c(1,112,124),]
#  myPeaks1$peaks[c(1,112,124),]
#  myPeaks1$counts[c(1,112,124),]
  
  expect_is(myPeaks1, 'list')
  
})




# pos[0]: 848853, pos[ pos.size() - 1 ]: 4.50364e+06, min_pos: 800001, max_pos: 4600000

test_that("freq_peak works, first window doesn't begin at 1",{
  data("vcfR_test")
  myPos <- getPOS(vcfR_test)
  myPos <- myPos + 848853 - myPos[1]
  myPos[ length(myPos) ] <- 4503640
  vcfR_test@fix[,'POS'] <- myPos
  
  dp <- extract.gt(vcfR_test, element = "DP", as.numeric = TRUE)
  myPeaks1 <- freq_peak(dp, pos = getPOS(vcfR_test), winsize = 1e5, lhs = FALSE)
  
  expect_equal( sum(myPeaks1$wins[-1,"START"] > myPeaks1$wins[-nrow(myPeaks1$wins),"END"]), nrow(myPeaks1$wins) - 1 )
})



##### ##### ##### ##### #####
# EOF