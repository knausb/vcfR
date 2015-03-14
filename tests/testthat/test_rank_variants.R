# rank_variants tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("rank_variants functions")


data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)


#.Call('vcfR_rank_variants', PACKAGE = 'vcfR', variants, ends, score)

vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', pinf_mt@var.info, pinf_mt@win.info$end, pinf_mt@var.info$MQ)

test_that("Rank variants binary is working",{
  expect_equal(names(vars)[ncol(vars)], "window_number")
  expect_equal(length(vars$window_number), 247)
  expect_equal(length(vars$window_number[vars$window_number == 19]), 1)
})

#head(vars)






