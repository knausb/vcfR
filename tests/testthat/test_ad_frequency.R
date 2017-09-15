# detach(package:vcfR, unload=T)
#library(testthat)
library(vcfR)

context("AD_frequency")

#data("vcfR_example")

test_that("AD_frequency works",{
  set.seed(999)
  x1 <- round(rnorm(n=9, mean=10, sd=2))
  x2 <- round(rnorm(n=9, mean=20, sd=2))
  ad <- matrix(paste(x1, x2, sep=","), nrow=3, ncol=3)
  colnames(ad) <- paste('Sample', 1:3, sep="_")
  rownames(ad) <- paste('Variant', 1:3, sep="_")
  ad[1,1] <- "9,23,12"
  my_ad <- AD_frequency(ad=ad)
  
  test <- as.numeric(unlist(strsplit(ad[1,1], split=",")))
  test <- sort(test, decreasing = TRUE)
  
  expect_equal( round(my_ad[1,1], digits=6), round(test[1]/sum(test[1:2]), digits=6) )
  
  my_ad <- AD_frequency(ad = ad, allele = 2)
  expect_equal( round(my_ad[1,1], digits=6), round(test[2]/sum(test[1:2]), digits=6) )
  
  my_ad <- AD_frequency(ad = ad, allele = 3)
  expect_equal( round(my_ad[1,1], digits=6), round(test[3]/sum(test[1:2]), digits=6) )
  
  my_ad <- AD_frequency(ad = ad, allele = 1, sum_type=1)
  expect_equal( round(my_ad[1,1], digits=6), round(test[1]/sum(test), digits=6) )

})


