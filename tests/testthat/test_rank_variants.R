# rank_variants tests.

# detach(package:vcfR, unload=T)


library(vcfR)
#library(testthat)
context("rank_variants functions")

data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE, win.size=1e3)

set.seed(1)

testv <- runif(nrow(pinf_mt@var.info), 0, 40)
vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', pinf_mt@var.info, pinf_mt@win.info$end, testv)


#.Call('vcfR_rank_variants', PACKAGE = 'vcfR', variants, ends, score)

#vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', pinf_mt@var.info, pinf_mt@win.info$end, pinf_mt@var.info$MQ)




#cbind(vars[,c(1:4, 24:25)], testv)[1:28,]



#cbind(vars[,c(1:4,24:25)], testv)[1:18,]


#vars <- .Call('vcfR_pair_sort', PACKAGE = 'vcfR')



# cbind(vars$window_number, pinf_mt@var.info$MQ, order(vars$window_number, pinf_mt@var.info$MQ, decreasing=T))
# cbind(pinf_mt@var.info$MQ, vars$window_number, order(vars$window_number, pinf_mt@var.info$MQ, decreasing=T))


test_that("Rank variants binary is working",{
#  expect_equal(names(vars)[ncol(vars)], "window_number")
  expect_equal(names(vars)[ncol(vars)], "rank")
  expect_equal(length(unique(vars$window_number)), 39)
  expect_equal(length(vars$window_number[vars$window_number == 19]), 5)
})


test_that("rank_variants_chrom is working",{
  expect_is(pinf_mt, "Chrom")

  set.seed(1)
  scores <- runif(nrow(pinf_mt@var.info), min=1, max=10)
  
  tmp <- rank_variants_chrom(pinf_mt, scores)
#  head(tmp@var.info)
  expect_equal(length(grep("window_number", names(tmp@var.info))), 1)
  expect_equal(length(grep("rank", names(tmp@var.info))), 1)
})


#head(vars)






