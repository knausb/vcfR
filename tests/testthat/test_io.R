#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("io functions")

# Load data
data(vcfR_example)
pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf_mt <- proc_chrom(pinf_mt, win.size=1000, verbose=FALSE)


# Manage directories.
original_dir <- getwd()
test_dir <- tempdir()

##### ##### ##### ##### #####

test_that("read/write.vcf works for vcfR objects",{

  setwd(test_dir)
  write.vcf(pinf_vcf, "test.vcf")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
#  expect_is(test@fix$QUAL, "integer") # Old, pure R version.
  expect_is(test@fix$QUAL, "numeric")
  expect_identical(names(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), 371)
  expect_equal(ncol(test@gt), 30)

})


test_that("read/write.vcf works for Chrom objects",{
  
  setwd(test_dir)
  write.vcf(pinf_mt, "test.vcf")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
  #  expect_is(test@fix$QUAL, "integer") # Old, pure R version.
  expect_is(test@fix$QUAL, "numeric")
  expect_identical(names(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), 371)
  expect_equal(ncol(test@gt), 30)
  
})


test_that("write.vcf.gz works for Chrom objects",{
  
  setwd(test_dir)
  write.vcf.gz(pinf_mt, "test.vcf.gz")
  system("gunzip test.vcf.gz")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
  #  expect_is(test@fix$QUAL, "integer") # Old, pure R version.
  expect_is(test@fix$QUAL, "numeric")
  expect_identical(names(test@fix)[1], "CHROM")
  expect_equal(nrow(test@gt), 371)
  expect_equal(ncol(test@gt), 30)
  
})



test_that("write_var_info works for Chrom objects",{
  
  setwd(test_dir)
  write_var_info(pinf_mt, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
  setwd(original_dir)
  
  expect_is(test, "data.frame")
  expect_equal(nrow(test), 371)
  expect_equal(ncol(test), 9)
  expect_equal(length(grep("CHROM", names(test))), 1)
  expect_equal(length(grep("POS", names(test))), 1)
  expect_equal(length(grep("mask", names(test))), 1)
  
})


test_that("write_var_info works for Chrom objects",{
  
  setwd(test_dir)
  write_win_info(pinf_mt, "test.csv")
  test <- read.table("test.csv", header=TRUE, sep=",")
  unlink("test.csv")
  setwd(original_dir)

  expect_is(test, "data.frame")
  expect_equal(nrow(test), 40)
  expect_equal(ncol(test), 12)
  expect_equal(length(grep("window", names(test))), 1)
  expect_equal(length(grep("start", names(test))), 1)
  expect_equal(length(grep("end", names(test))), 1)

})


test_that("read.vcf works for vcf files which contain no variants",{
  
  
  pinf_vcf2 <- pinf_vcf
  pinf_vcf2@fix <- pinf_vcf2@fix[0,]
  pinf_vcf2@gt <- pinf_vcf2@gt[0,]
  
  setwd(test_dir)
  
  write.vcf.gz(pinf_vcf2, "test.vcf.gz")
  system("gunzip test.vcf.gz")
  test <- read.vcf("test.vcf")
  unlink("test.vcf")
  
  setwd(original_dir)

  expect_equal(ncol(test@fix), 8)
  expect_equal(ncol(test@gt), 30)
  expect_equal(nrow(test@fix), 0)
  expect_equal(nrow(test@gt), 0)
  
})





#.Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "test.vcf.gz")
#.Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', "../vcf_data/gatk_hc/sc_1.100.vcf.gz")
#.Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', pinf_vcf@fix, pinf_vcf@gt, "test.vcf.gz", 0)



