library(vcfR)
context("var_window functions")

data(vcfR_example)

win1 <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(pinf_dna))

test_that("vcfR_window_init creates ", {
#  win1 <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(pinf_dna))
  expect_equal(ncol(win1), 4)
  expect_equal(nrow(win1), 40)
})


test_that("vcfR_window_init creates ", {
#  win1 <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(pinf_dna))
  win2 <- .Call('vcfR_windowize_fasta', PACKAGE = 'vcfR', wins=win1, seq=as.character(pinf_dna)[1,])
  expect_equal(ncol(win2), 10)
  expect_equal(nrow(win2), 40)
  expect_equal(rowSums(win2[,5:10]), c(rep(1000, times=39), 870))
})



win3 <- .Call('vcfR_windowize_variants', PACKAGE = 'vcfR', wins=win2, pos=pinf_vcf@fix$POS)

