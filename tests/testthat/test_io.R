#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("io functions")

# Load data

data(vcfR_example)

#vcf_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
#seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")
#gff_file <- system.file("extdata", "pinf_sc100.gff", package = "vcfR")

#vcf <- read.vcfR(vcf_file, verbose = FALSE)
#dna <- ape::read.dna(seq_file, format = "fasta")
#gff <- read.table(gff_file, sep="\t")

chrom <- create.chromR(name="Supercontig_1.50", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)


# Manage directories.
original_dir <- getwd()
test_dir <- tempdir()

##### ##### ##### ##### #####

test_that("vcfR_vcf_stats_gz works",{
  setwd(test_dir)
  write.vcf(chrom, "test.vcf.gz")
  x <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "test.vcf.gz")
  #  test <- read.vcfR("test.vcf")
  unlink("test.vcf.gz")
  setwd(original_dir)
  
  expect_equal(as.numeric(x["meta"]), length(chrom@vcf@meta))
  expect_equal(as.numeric(x["header"]), c(length(chrom@vcf@meta) + 1) )
  expect_equal(as.numeric(x["variants"]), nrow(chrom@vcf@fix))
  expect_equal(as.numeric(x["columns"]), ncol(chrom@vcf@fix) + ncol(chrom@vcf@gt))
})


test_that("vcfR_vcf_meta_gz works",{
  setwd(test_dir)
  write.vcf(chrom, "test.vcf.gz")
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "test.vcf.gz")
  x <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', "test.vcf.gz", stats, 0)
  unlink("test.vcf.gz")
  setwd(original_dir)
  
  expect_equal(length(x), length(chrom@vcf@meta))
})




test_that("vcfR_read_body_gz works",{
  setwd(test_dir)
  write.vcf(chrom, "test.vcf.gz")
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "test.vcf.gz")
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', "test.vcf.gz", stats,
                nrows = -1, skip = 0, cols=1:stats['columns'], 0)
  unlink("test.vcf.gz")
  setwd(original_dir)

  expect_equal(colnames(body)[1], "CHROM")
  expect_equal(ncol(body), as.integer(stats['columns']))
  expect_equal(nrow(body), as.integer(stats['variants']))
  
})


test_that("read/write.vcf works for vcfR objects",{
  setwd(test_dir)
  write.vcf(vcf, "test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
  setwd(original_dir)
  
  expect_is(test, "vcfR")

  expect_identical(colnames(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), nrow(vcf@gt))
  expect_equal(ncol(test@gt), ncol(vcf@gt))
})



test_that("write.vcf APPEND=TRUE does not include header",{
  setwd(test_dir)
  write.vcf(vcf, "test.vcf.gz", APPEND=TRUE)
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
  setwd(original_dir)

})

#test_that("read/write.vcf works for Chrom objects",{
#  setwd(test_dir)
#  write.vcf(chrom, "test.vcf")
#  test <- read.vcfR("test.vcf", verbose = FALSE)
#  unlink("test.vcf")
#  setwd(original_dir)
  
#  expect_is(test, "vcfR")
#  expect_identical(colnames(test@fix)[1], "CHROM")
#  expect_equal(nrow(test@gt), nrow(vcf@gt))
#  expect_equal(ncol(test@gt), ncol(vcf@gt))  
#})


test_that("write.vcf.gz works for Chrom objects",{
  
  setwd(test_dir)
  write.vcf(chrom, "test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_identical(colnames(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), nrow(vcf@gt))
  expect_equal(ncol(test@gt), ncol(vcf@gt))
})


test_that("write.var.info works for Chrom objects",{
  setwd(test_dir)
  write.var.info(chrom, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
  setwd(original_dir)
  
  expect_is(test, "data.frame")
  expect_equal(nrow(test), nrow(vcf@fix))
  expect_equal(length(grep("CHROM", colnames(test))), 1)
  expect_equal(length(grep("POS", colnames(test))), 1)
  expect_equal(length(grep("mask", colnames(test))), 1)  
})


test_that("write.win.info works for Chrom objects",{
  setwd(test_dir)
  write.win.info(chrom, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
  setwd(original_dir)

  expect_is(test, "data.frame")
  expect_equal(nrow(test), nrow(chrom@win.info))
#  expect_equal(ncol(test), 12)
  expect_equal(grep("CHROM", names(test), value=TRUE), "CHROM")
  expect_equal(grep("window", names(test), value=TRUE), "window")
  expect_equal(grep("start", names(test), value=TRUE), "start")
  expect_equal(grep("end", names(test), value=TRUE), "end")
#  expect_equal(length(grep("window", names(test))), 1)
#  expect_equal(length(grep("start", names(test))), 1)
#  expect_equal(length(grep("end", names(test))), 1)
})





test_that("write_fasta works",{

#invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, seqname, filename, rowlength, verbose))
#invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', as.character(pinf_dna)[1,], "myseq", "test.fasta", 10, 1))

#  write_fasta(pinf_mt, file="pinf_mt.fasta", gt_split="/", rowlength=1141)
})



#data(vcfR_example)
#write.vcf.gz(pinf_vcf, "test.vcf.gz")

#.Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "test.vcf.gz")
#.Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "../vcf_data/gatk_hc/sc_1.100.vcf.gz")
#.Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', pinf_vcf@fix, pinf_vcf@gt, "test.vcf.gz", 0)



