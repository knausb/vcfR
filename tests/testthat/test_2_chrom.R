# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("create_chrom functions")

#data(vcfR_example)

test_that("we can create null a Chrom",{
  chrom <- new(Class="Chrom")
  expect_is(chrom, "Chrom")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")

  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix), 0)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})



vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

vcf <- read.vcf(vcf_file, verbose = FALSE)
dna <- ape::read.dna(seq_file, format = "fasta")
gff <- read.table(gff_file, sep="\t")


test_that("We can create a Chrom, no sequence or annotation",{
  chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, verbose=FALSE)
  expect_is(chrom, "Chrom")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")

  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})


test_that("We can create a Chrom, no annotation",{
  chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, seq=dna, verbose=FALSE)
  expect_is(chrom, "Chrom")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "DNAbin")
  expect_is(chrom@ann, "data.frame")

  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq)>0, TRUE)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann), 0)
})


test_that("We can create a Chrom, no sequence",{
  chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, ann=gff, verbose=FALSE)
  expect_is(chrom, "Chrom")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "NULL")
  expect_is(chrom@ann, "data.frame")

  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq), 0)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann)>0, TRUE)
})


test_that("We can create a Chrom",{
  chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  expect_is(chrom, "Chrom")
  expect_is(chrom@vcf, "vcfR")
  expect_is(chrom@seq, "DNAbin")
  expect_is(chrom@ann, "data.frame")

  expect_equal(ncol(chrom@vcf@fix), 8)
  expect_equal(nrow(chrom@vcf@fix) > 0, TRUE)
  expect_equal(length(chrom@seq)>0, TRUE)
  expect_equal(ncol(chrom@ann), 9)
  expect_equal(nrow(chrom@ann)>0, TRUE)
})


