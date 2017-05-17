library(vcfR)
context("var_window functions")

data(vcfR_example)

# Global code.
#pinf_mt <- create.chromR('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf <- create.chromR('pinf', seq=dna, vcf=vcf, ann=gff, verbose=FALSE)


win1 <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(dna))
win2 <- .Call('vcfR_windowize_fasta', PACKAGE = 'vcfR', wins=win1, seq=as.character(dna)[1,])
win3 <- .Call('vcfR_windowize_variants', PACKAGE = 'vcfR', windows=win2, variants=pinf@var.info[c('POS','mask')])
win4 <- .Call('vcfR_windowize_annotations', PACKAGE = 'vcfR', wins=win3,
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


