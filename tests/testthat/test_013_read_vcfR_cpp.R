
#library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)

#
context("read vcfR C++ functions")



##### ##### ##### ##### #####
#
# vcf_stats
#
##### ##### ##### ##### #####

myDir <- getwd()
temp_dir <- tempdir()

setwd( temp_dir )

test_that("compiled vcfR_vcf_stats_gz works",{
  data("vcfR_test")
  write.vcf(vcfR_test, file="myFile.vcf.gz")
  
  stats <- .vcf_stats_gz("myFile.vcf.gz", verbose = 0)
  unlink("myFile.vcf.gz")
  
  expect_equal( as.numeric(stats['meta']), length(vcfR_test@meta) )
  expect_equal( as.numeric(stats['variants']), nrow(vcfR_test@fix) )
  expect_equal( as.numeric(stats['columns']), ncol(vcfR_test@fix) + ncol(vcfR_test@gt))
})


test_that("compiled vcfR_vcf_stats_gz works, Windows carriage return",{
  data("vcfR_test")
  ex_file <- "myFile.vcf.gz"
  
  # Create e file with carriage returns.
  # Write it to a file.
  cat(vcfR_test@meta[1], file=ex_file, append = FALSE, sep = '\t')
  for(i in 2:length(vcfR_test@meta) ){
    cat('\r\n', file=ex_file, append = TRUE, sep = '\t')
#    cat('\n', file=ex_file, append = TRUE)
    cat(vcfR_test@meta[i], file=ex_file, append = TRUE, sep = '\t')
  }
  cat('\r\n', file=ex_file, append = TRUE, sep = '\t')
#  cat('\n', file=ex_file, append = TRUE)
  cat('#', file=ex_file, append = TRUE, sep = '\t')
  cat(colnames(vcfR_test@fix), file=ex_file, append = TRUE, sep = '\t')
  cat('\t', file=ex_file, append = TRUE)
  cat(colnames(vcfR_test@gt), file=ex_file, append = TRUE, sep = '\t')
  for(i in 1:nrow(vcfR_test@fix)){
    cat('\r\n', file=ex_file, append = TRUE, sep = '\t')
#    cat('\n', file=ex_file, append = TRUE)
    cat(vcfR_test@fix[i,], file=ex_file, append = TRUE, sep = '\t')
    cat(vcfR_test@gt[i,], file=ex_file, append = TRUE, sep = '\t')
  }
  cat('\r\n', file=ex_file, append = TRUE, sep = '\t')
#  cat('\n', file=ex_file, append = TRUE)
  
  # Query file.
  stats <- .vcf_stats_gz(ex_file, verbose = 0)

  # Test.
  expect_equal( as.numeric(stats['meta']), length(vcfR_test@meta) )
  expect_equal( as.numeric(stats['variants']), nrow(vcfR_test@fix) )
  expect_equal( as.numeric(stats['columns']), ncol(vcfR_test@fix) + ncol(vcfR_test@gt))

})


test_that("compiled vcf_stats_gz nrows works",{
#  data("vcfR_test")
  data("vcfR_example")
  write.vcf(vcf, file="myFile.vcf.gz")
  
  stats <- .vcf_stats_gz("myFile.vcf.gz", nrows = 10, verbose = 0)
  unlink("myFile.vcf.gz")
  
  expect_equal( as.numeric(stats['variants']), 10 )
})


test_that("compiled vcf_stats_gz skip works",{
#  data("vcfR_test")
  data("vcfR_example")
  write.vcf(vcf, file="myFile.vcf.gz")
  
  stats <- .vcf_stats_gz("myFile.vcf.gz", skip = 2000, verbose = 0)
  unlink("myFile.vcf.gz")
  
  expect_equal( as.numeric(stats['variants']), 2533 )
  # Skip doesn't do anything in this context.
  # It is added to nrows so it is needed.
})



##### ##### ##### ##### #####
#
# vcfR_read_body
#
##### ##### ##### ##### #####

test_that("compiled vcfR_read_body works when file contains no variants",{

  data("vcfR_test")
  vcf2 <- vcfR_test[0,]
  
  myFile <- "myFile.vcf.gz"
  write.vcf(vcf2, myFile)
  stats <- .vcf_stats_gz(myFile, verbose = 0)
  
  body <- .read_body_gz(myFile, stats,
                nrows = 0, skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
#  
  expect_equal( ncol(body), ncol(vcf2@fix) + ncol(vcf2@gt) )
#
  expect_equal( nrow(body), nrow(vcf2@fix) )
  
})


##### ##### ##### ##### #####
#
# vcfR_read
#
##### ##### ##### ##### #####

test_that("compiled input functions work",{
  data("vcfR_example")
  ex_file <- "test.vcf.gz"
  write.vcf(vcf, file=ex_file)

  stats <- .vcf_stats_gz(ex_file, verbose = 0)

  expect_equal(length(stats), 4)
  expect_is(stats, "numeric")
  
  meta <- .read_meta_gz(ex_file, stats, 0)
  expect_equal(length(meta), as.numeric(stats["meta"]))
  expect_is(meta, "character")
  
  body <- .read_body_gz(ex_file, stats,
                        nrows = -1, skip = 0, cols=1:stats['columns'],
                        convertNA = 1, verbose = 0)

  expect_is(body, "matrix")
  expect_equal(nrow(body), as.numeric(stats["variants"]))
  expect_equal(ncol(body), as.numeric(stats["columns"]))
  
  # Check for cariage return.
#  body[1,ncol(body)] <- "0/0:32,0:32:90:0,90,999\r" # Test example.
  expect_equal(grep("\r$", body[1,ncol(body)]), integer(0))
  
  unlink(ex_file)
})


test_that("compiled vcfR_read_body converts VCF missing to NA",{
  data("vcfR_test")
  vcfR_test@gt[1,3] <- "./."
  
  myFile <- "myFile.vcf.gz"
  write.vcf(vcfR_test, myFile)

  stats <- .vcf_stats_gz(myFile, verbose = 0)
  body <- .read_body_gz(myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
  expect_true( is.na(body[1,11]) )
  
  vcfR_test@gt[1,3] <- ".|."
  write.vcf(vcfR_test, myFile)

  stats <- .vcf_stats_gz(myFile, verbose = 0)
  body <- .read_body_gz(myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
  expect_true( is.na(body[1,11]) )
#  body
  
  vcfR_test@gt[1,3] <- "./0"
  write.vcf(vcfR_test, myFile)
  
  stats <- .vcf_stats_gz(myFile, verbose = 0)
  body <- .read_body_gz(myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)

  expect_false( is.na(body[1,11]) )
})


test_that("compiled vcfR_read_body convertNA = FALSE works",{
  data("vcfR_test")
#  vcfR_test@fix[1,6] <- "."
  vcfR_test@gt[1,3] <- "./."
  vcfR_test@gt[1,4] <- ".|."
  
  myFile <- "myFile.vcf.gz"
  write.vcf(vcfR_test, myFile)

  stats <- .vcf_stats_gz(myFile, verbose = 0)
  body <- .read_body_gz(myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 0, verbose = 0)
  unlink(myFile)

  expect_false( is.na(body[2,3]) )
  expect_false( is.na(body[1,11]) )
  expect_false( is.na(body[1,12]) )
})


test_that("VCF with no GT, compiled functions",{
  data("vcfR_test")
#  vcf@gt <- vcf@gt[-c(1:nrow(vcf@gt)),]
  vcf@gt <- matrix("a", nrow=0, ncol=0)
  
  ex_file <- "test.vcf.gz"
  
  test <- .write_vcf_body(fix = vcf@fix, gt = vcf@gt, filename = ex_file, mask = 0)
  stats <- .vcf_stats_gz(ex_file, verbose = 0)
  body <- .read_body_gz(ex_file, stats,
                nrows = -1, skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)

  expect_equal(test, NULL)
  expect_is( body, "matrix" )
  expect_equal( nrow(body), nrow(vcf@fix))
  expect_equal( ncol(body), ncol(vcf@fix))
  
  unlink(ex_file)
})


##### ##### ##### ##### #####
# EOF.