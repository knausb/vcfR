


#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)


#
context("genetic_diff")

##### ##### ##### ##### #####
# Jost's D

test_that("Jost's example works",{
  data("vcfR_test")
  
  # Create VCF data.
  jost <- vcfR_test[1,]
  jost@gt <- matrix(nrow=1, ncol=221)
  jost@gt[1,1] <- "GT"
  jost@gt[,2:11] <- "0/1"
  jost@gt[,12:21] <- "2/3"
  jost@gt[,22:121] <- "2/3"
  jost@gt[,122:221] <- "4/5"
  colnames(jost@gt) <- c("FORMAT", paste("sample", 1:220, sep="_"))
  
  # Pop factor
  myPops <- rep("b", times=220)
  myPops[1:20] <- "a"
  myPops <- as.factor(myPops)
  
  # genetic_diff
  tmp <- genetic_diff(jost, myPops, method = "jost")
  
  expect_equal(trunc(1e2*tmp$a), 25)
  expect_equal(trunc(1e7*tmp$b), 4788895)
  expect_equal(trunc(1e7*tmp$Dest_Chao), 4779589)
})


##### ##### ##### ##### #####
# Nei's Gst


test_that("Nei's method works",{
  #
  devtools::load_all(".")
  #
  debug("calc_nei")
  data("vcfR_test")
  vcfR_test@gt <- cbind(vcfR_test@gt, vcfR_test@gt[,2:4])
  myPops <- as.factor(rep(c('a','b'), each=3))
  
  tmp <- genetic_diff(vcfR_test, myPops, method = "nei")
  

  
})



