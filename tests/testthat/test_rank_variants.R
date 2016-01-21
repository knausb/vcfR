# rank_variants tests.
# detach(package:vcfR, unload=T)
library(vcfR)
#library(testthat)
context("rank_variants functions")

#data(vcfR_example)

vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

vcf <- read.vcfR(vcf_file, verbose = FALSE)
dna <- ape::read.dna(seq_file, format = "fasta")
gff <- read.table(gff_file, sep="\t")


chrom <- create_chrom('pinf', seq=dna, vcf=vcf, ann=gff, verbose=FALSE)
pinf_mt <- proc_chrom(chrom, verbose=FALSE, win.size=1e3)

set.seed(1)
testv <- runif(nrow(chrom@var.info), 0, 40)
#vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', chrom@var.info, chrom@win.info$end, testv)


#.Call('vcfR_rank_variants', PACKAGE = 'vcfR', variants, ends, score)

#vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', pinf_mt@var.info, pinf_mt@win.info$end, pinf_mt@var.info$MQ)




#cbind(vars[,c(1:4, 24:25)], testv)[1:28,]



#cbind(vars[,c(1:4,24:25)], testv)[1:18,]


#vars <- .Call('vcfR_pair_sort', PACKAGE = 'vcfR')



# cbind(vars$window_number, pinf_mt@var.info$MQ, order(vars$window_number, pinf_mt@var.info$MQ, decreasing=T))
# cbind(pinf_mt@var.info$MQ, vars$window_number, order(vars$window_number, pinf_mt@var.info$MQ, decreasing=T))


test_that("Rank variants binary is working",{
#  expect_equal(names(vars)[ncol(vars)], "window_number")
#  expect_equal(names(vars)[ncol(vars)], "rank")
#  expect_equal(length(unique(vars$window_number)), 39)
#  expect_equal(length(vars$window_number[vars$window_number == 19]), 5)
})


test_that("rank_variants_chrom is working",{
  expect_is(chrom, "Chrom")

  set.seed(1)
  scores <- runif(nrow(chrom@var.info), min=1, max=10)
  
#  tmp <- rank_variants_chrom(chrom, scores)
#  head(tmp@var.info)
#  expect_equal(length(grep("window_number", names(tmp@var.info))), 1)
#  expect_equal(length(grep("rank", names(tmp@var.info))), 1)
})


#head(vars)






