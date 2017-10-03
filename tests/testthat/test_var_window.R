
library(vcfR)
context("var_window functions")

data(vcfR_example)


win1 <- .window_init(window_size=1e3, max_bp=length(dna))
win2 <- .windowize_fasta(wins=win1, seq=as.character(dna)[1,])
win3 <- .windowize_variants(windows=win2, variants=chrom@var.info[c('POS','mask')])
win4 <- .windowize_annotations(wins=win3,
              ann_starts=as.numeric(as.character(gff[,4])), 
              ann_ends=as.numeric(as.character(gff[,5])),
              chrom_length=length(dna))


test_that("vcfR_window_init works", {
  expect_equal(ncol(win1), 4)
  expect_equal(nrow(win1), ceiling(length(dna)/1e3))
})


test_that("vcfR_windowize_fasta works ", {
  expect_equal(ncol(win2), 10)
  expect_equal(nrow(win2), ceiling(length(dna)/1e3))
})


test_that("vcfR_windowize_variants works ", {
  expect_equal(ncol(win3), 11)
  expect_equal(nrow(win2), ceiling(length(dna)/1e3))
})


test_that("vcfR_windowize_annotations works ", {
  expect_equal(ncol(win4), 12)
  expect_equal(nrow(win2), ceiling(length(dna)/1e3))
})


