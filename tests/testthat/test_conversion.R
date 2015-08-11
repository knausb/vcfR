# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("windowing functions")


vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

vcf <- read.vcf(vcf_file, verbose = FALSE)
dna <- ape::read.dna(seq_file, format = "fasta")
gff <- read.table(gff_file, sep="\t")

chrom <- create_chrom(name="Supercontig_1.100", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)



gt <- extract.gt(chrom, element="DP", as.numeric=TRUE)

test_that("gt is numeric",{
  expect_is(gt, "matrix")
  expect_equal(is.numeric(gt), TRUE)
})


#head(pinf_mt@vcf.gt)

#head(gt)


