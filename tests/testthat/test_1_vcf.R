
#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("vcf functions")

#ex_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")

data("vcfR_example")
tot_var <- nrow(vcf@gt)

##### ##### ##### ##### #####
# Manage directories.

original_dir <- getwd()
test_dir <- tempdir()

setwd( test_dir )

#ex_file <- paste(test_dir, "/test.vcf.gz", sep="")
ex_file <- "test.vcf.gz"

write.vcf(vcf, file=ex_file)


##### ##### ##### ##### #####
#
# vcfR class.
#
##### ##### ##### ##### #####


test_that("We can create an empty vcfR object",{
  myvcf <- new("vcfR")
  
  expect_is(myvcf@meta, "character")
  expect_equal(length(myvcf@meta), 0)
  
  expect_is(myvcf@fix, "matrix")
  expect_equal(nrow(myvcf@fix), 0)
  expect_equal(ncol(myvcf@fix), 8)
  
  expect_is(myvcf@gt, "matrix")
  expect_equal(nrow(myvcf@gt), 0)
  expect_equal(ncol(myvcf@gt), 0)
  
})


##### ##### ##### ##### #####
#
# Input functions
#
##### ##### ##### ##### #####

test_that("compiled input functions work",{
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  expect_equal(length(stats), 4)
  expect_is(stats, "numeric")
  
  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
  expect_equal(length(meta), as.numeric(stats["meta"]))
  expect_is(meta, "character")
  
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)
  expect_is(body, "matrix")
  expect_equal(nrow(body), as.numeric(stats["variants"]))
  expect_equal(ncol(body), as.numeric(stats["columns"]))
  
  # Check for cariage return.
#  body[1,ncol(body)] <- "0/0:32,0:32:90:0,90,999\r" # Test example.
  expect_equal(grep("\r$", body[1,ncol(body)]), integer(0))
  
})



test_that("compiled vcfR_read_body works when file contains no variants",{
  vcf2 <- vcf
  vcf2@fix <- vcf2@fix[0,]
  vcf2@gt <- vcf2@gt[0,]
  
  write.vcf(vcf2, ex_file)
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
#  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = 0, skip = 0, cols=1:stats['columns'], 0)
  
  expect_equal( ncol(body), ncol(vcf2@fix) + ncol(vcf2@gt) )
  expect_equal( nrow(body), nrow(vcf2@fix) )
  
  unlink(ex_file)

})


##### ##### ##### ##### #####


data("vcfR_example")
write.vcf(vcf, file=ex_file)


test_that("read.vcfR works",{
  vcf <- read.vcfR(ex_file, verbose=FALSE)
  expect_is(vcf, "vcfR")
  expect_is(vcf@meta, "character")
  expect_is(vcf@fix, "matrix")
  expect_is(vcf@gt, "matrix")
})


test_that("read.vcfR nrows works",{
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, nrows=count)
  expect_equal(nrow(vcf), count)
})


test_that("read.vcfR skip works",{
  count <- 100
  vcf <- read.vcfR(ex_file, verbose=FALSE, skip=count)
  expect_equal(nrow(vcf), tot_var - count)
})


test_that("read.vcfR column selection works",{
  vcf <- read.vcfR(ex_file, verbose=FALSE, cols=c(9,11:12))
  expect_is(vcf@gt, "matrix")
  expect_equal(ncol(vcf@fix), 8)
  expect_equal(ncol(vcf@gt), 3)
})


test_that("read.vcfR works for vcf files which contain no variants",{  
  vcf2 <- vcf
  vcf2@fix <- vcf2@fix[0,]
  vcf2@gt <- vcf2@gt[0,]
  
  write.vcf(vcf2, ex_file)
  test <- read.vcfR(ex_file, verbose=FALSE)
  unlink(ex_file)


  expect_equal(ncol(test@fix), ncol(vcf2@fix))
  expect_equal(ncol(test@gt), ncol(vcf2@gt))
  expect_equal(nrow(test@fix), nrow(vcf2@fix))
  expect_equal(nrow(test@gt), nrow(vcf2@gt))
})


##### ##### ##### ##### #####
#
# Write funcitons work.
#
##### ##### ##### ##### #####

#vcf <- read.vcfR(ex_file, verbose=FALSE)

test_that("write.vcf works on objects of class vcfR",{
  write.vcf(vcf, file="test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")

  expect_equal(nrow(vcf), nrow(test))
  expect_equal(ncol(vcf@fix), ncol(test@fix))
  expect_equal(ncol(vcf@gt), ncol(test@gt))
})


test_that("write.vcf works on objects of class chromR",{
  suppressWarnings(chrom <- create.chromR(vcf=vcf, seq = dna, ann = gff))

  write.vcf(chrom, file="test.vcf.gz")
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")
  
  expect_equal(nrow(vcf), nrow(test))
  expect_equal(ncol(vcf@fix), ncol(test@fix))
  expect_equal(ncol(vcf@gt), ncol(test@gt))
})

test_that("write.vcf works on objects of class chromR when mask=TRUE",{
  suppressWarnings(chrom <- create.chromR(vcf=vcf, seq = dna, ann = gff))
  chrom <- masker(chrom, min_DP = 250, max_DP = 750, min_MQ = 59.5, max_MQ = 60.5)
  
  write.vcf(chrom, file="test.vcf.gz", mask = TRUE)
  test <- read.vcfR("test.vcf.gz", verbose = FALSE)
  unlink("test.vcf.gz")

  expect_equal(sum(chrom@var.info$mask), nrow(test))
})


##### ##### ##### ##### #####


test_that("vcfR subsetters works",{
  # Rows
  vcf2 <- vcf[1:10,]
  expect_equal(nrow(vcf2@fix), 10)
  expect_equal(nrow(vcf2@gt), 10)
  
  # Columns
  vcf2 <- vcf[,1:4]
  expect_equal(ncol(vcf2@fix), 8)
  expect_equal(ncol(vcf2@gt), 4)
})

  
#debug(read.vcfR)
#unlink(ex_file)


##### ##### ##### ##### #####
#
# IO works for VCF data that 
# contains no gt region.
#
##### ##### ##### ##### #####


test_that("VCF with no GT, compiled functions",{
#  vcf@gt <- vcf@gt[-c(1:nrow(vcf@gt)),]
  vcf@gt <- matrix("a", nrow=0, ncol=0)
  
  ex_file <- "test.vcf.gz"
  
  test <- .Call("vcfR_write_vcf_body", PACKAGE = "vcfR", 
            fix = vcf@fix, gt = vcf@gt, filename = ex_file, mask = 0)
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)

  expect_equal(test, NULL)
  expect_is( body, "matrix" )
  expect_equal( nrow(body), nrow(vcf@fix))
  expect_equal( ncol(body), ncol(vcf@fix))
  
  unlink("test.vcf.gz")
  
})


test_that("VCF with no GT",{
  vcf@gt <- matrix("a", nrow=0, ncol=0)
  ex_file <- "test.vcf.gz"
  
  test <- write.vcf(vcf, file = ex_file)
  expect_equal(test, NULL)
  
#  debug(read.vcfR)
  test <- read.vcfR(ex_file, verbose = FALSE)
  unlink("test.vcf.gz")
    
  expect_is(test, "vcfR")  
  expect_is(test@fix, "matrix")
  
  expect_equal(nrow(test@fix), nrow(vcf@fix))
  expect_equal(ncol(test@fix), 8)
  
  expect_is(test@gt, "matrix")
  expect_equal( ncol(test@gt), 0 )
  expect_equal( nrow(test@gt), 0 )
})


##### ##### ##### ##### #####
# Manage directories.

setwd(original_dir)

##### ##### ##### ##### #####
# EOF.