
# detach(package:vcfR, unload=T)
#library(testthat)
library(vcfR)
context("create.chromR functions")

#library(testthat)
#data(vcfR_example)

test_that("Create a null chromR",{
  data("vcfR_example")
  chrom <- methods::new(Class="chromR")
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix), 0)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})


test_that("We can create a Chrom, no sequence or annotation",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, verbose=FALSE)
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")
  expect_is(chrom@var.info, "data.frame")
  expect_is(chrom@var.info[,'POS'], "integer")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
  expect_equal(ncol(chrom@var.info), 5)
  expect_equal(nrow(chrom@var.info)>0, TRUE)
})


test_that("We can create a chromR, no annotation",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, seq=dna, verbose=FALSE)
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "DNAbin")
  expect_is(chrom@ann, "data.frame")
  expect_is(chrom@var.info, "data.frame")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq)>0, TRUE)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})


test_that("We can create a chromR, no sequence",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, ann=gff, verbose=FALSE)
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")
  expect_is(chrom@var.info, "data.frame")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann)>0, TRUE)
})


test_that("We can create a chromR, no sequence, annotation greater than vcf POS",{
  data("vcfR_example")
  gff2 <- gff
  gff2[23,5] <- 100000
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, ann=gff2, verbose=FALSE)
  expect_equal( chrom@len, 100000)
})


test_that("We can create a chromR",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  expect_is(chrom, "chromR")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "DNAbin")
  expect_is(chrom@ann, "data.frame")
  expect_is(chrom@var.info, "data.frame")
  
  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq)>0, TRUE)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann)>0, TRUE)
})


##### ##### ##### ##### #####
# masker


test_that("We implemented the mask",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  chrom <- masker(chrom, min_DP = 300, max_DP = 700)
  expect_true( sum(chrom@var.info[,'mask']) < nrow(chrom@var.info) )
})


test_that("preserve = TRUE preserves the original mask ", {
  data("chromR_example")
  chrom_preserve_mask <- masker(chrom, min_QUAL = 40, preserve = TRUE)
  chrom_new_mask <- masker(chrom, min_QUAL = 40, preserve = FALSE)
  ori_mask <- chrom@var.info[, 'mask']
  add_mask <- chrom_preserve_mask@var.info[, 'mask']
  # adding restrictions reduces variants
  expect_true(sum(add_mask) < sum(ori_mask))
  # on comparison, the original masked variants are kept
  expect_equal(sum(ori_mask | add_mask), sum(ori_mask))
})



##### ##### ##### ##### #####
# proc.chromR


test_that("proc.chromR works",{
  data("chromR_example")  
  chrom <- proc.chromR(chrom, verbose = FALSE)
  
  expect_true( ncol(chrom@var.info) >= 3 )
  chrom <- proc.chromR(chrom, verbose = FALSE)
  expect_true( ncol(chrom@var.info) >= 3 )
})


##### ##### ##### ##### #####
# seq2rects


test_that("seq2rects works for test data",{
  data("chromR_example")
  rects1 <- seq2rects(chrom)
  
  expect_is( rects1, "matrix" )
  expect_true( nrow(rects1) > 0 )
})


test_that("seq2rects works with ns",{
  data("chromR_example")
  rects2 <- seq2rects(chrom, chars="n")
  
  expect_is( rects2, "matrix" )
  expect_true( nrow(rects2) > 0 )
})


test_that("seq2rects works when seq has no Ns",{
  data("chromR_example")
  # Replace n with a.
  seq2 <- as.character( chrom@seq )
  seq2[ seq2 == 'n' ] <- 'a'
  chrom@seq <- ape::as.DNAbin(seq2)
  
  rects2 <- seq2rects(chrom, chars="n")
  
  expect_is( rects2, "matrix" )
  expect_true( nrow(rects2) == 0 )
})


##### ##### ##### ##### #####
# chromR2vcfR

test_that("chromR2vcfR works",{
  data('vcfR_example')
  chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  chrom <- masker(chrom, min_DP = 300, max_DP = 700)
  
  test <- chromR2vcfR(chrom, use.mask = TRUE)

  expect_equal(sum(chrom@var.info$mask), nrow(test))
})

##### ##### ##### ##### #####
# EOF.