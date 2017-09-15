#
# rank_variants tests.
# detach(package:vcfR, unload=T)
library(vcfR)
#library(testthat)
context("rank_variants functions")

#data(vcfR_example)


##### ##### ##### ##### #####


test_that("rank.variants.chromR works",{
  data(vcfR_test)
  chrom <- create.chromR(name="Supercontig", vcf=vcfR_test, verbose=FALSE)
  chrom <- proc.chromR(chrom, verbose = FALSE)
  set.seed(9)
  scores <- runif(n=nrow(vcfR_test))
  chrom <- rank.variants.chromR(chrom, scores)
  expect_equal( length(chrom@var.info$rank), nrow(vcfR_test) )
})


test_that("rank.variants.chromR throws error when nvars != nscores",{
  data(vcfR_test)
  chrom <- create.chromR(name="Supercontig", vcf=vcfR_test, verbose=FALSE)
  chrom <- proc.chromR(chrom, verbose = FALSE)
  set.seed(9)
  scores <- runif(n=nrow(vcfR_test) + 1)
  
  msg <- "The number of variants and scores do not match."
#  msg <- paste(msg, " nrow(x@vcf): ", nrow(chrom@vcf), sep = "")
#  msg <- paste(msg, ", length(scores): ", length(scores), sep = "")
  
  expect_error(chrom <- rank.variants.chromR(chrom, scores), msg )
})



