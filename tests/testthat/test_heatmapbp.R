# detach(package:vcfR, unload=T)
#
library(testthat)
library(vcfR)

context("heatmap.bp function")

#data("vcfR_example")

test_that("heatmap.bp works",{
  tmp <- matrix( rnorm(100,1,10), ncol=10, nrow=10)
  test <- heatmap.bp(tmp)
  expect_equal(test, NULL)
})



test_that("heatmap.bp works, labels with small number of colors",{
  tmp <- matrix( rnorm(100,1,10), ncol=10, nrow=10)
  test <- heatmap.bp(tmp, col.ramp = viridisLite::viridis(n = 8, alpha = 1))
  expect_equal(test, NULL)
})


