
#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)


#
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
# vcf_stats
#
##### ##### ##### ##### #####


test_that("compiled vcfR_vcf_stats_gz works",{
  data("vcfR_test")
  write.vcf(vcfR_test, file=ex_file)
  
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)

  expect_equal( as.numeric(stats['meta']), length(vcfR_test@meta) )
  expect_equal( as.numeric(stats['variants']), nrow(vcfR_test@fix) )
  expect_equal( as.numeric(stats['columns']), ncol(vcfR_test@fix) + ncol(vcfR_test@gt))
})


test_that("compiled vcfR_vcf_stats_gz works, Windows carriage return",{
  data("vcfR_test")
  ex_file <- "test.vcf"
  
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
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  
  # Test.
  expect_equal( as.numeric(stats['meta']), length(vcfR_test@meta) )
  expect_equal( as.numeric(stats['variants']), nrow(vcfR_test@fix) )
  expect_equal( as.numeric(stats['columns']), ncol(vcfR_test@fix) + ncol(vcfR_test@gt))

})


##### ##### ##### ##### #####
#
# Input functions
#
##### ##### ##### ##### #####


data("vcfR_example")
ex_file <- "test.vcf.gz"
write.vcf(vcf, file=ex_file)

test_that("compiled input functions work",{
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)
  expect_equal(length(stats), 4)
  expect_is(stats, "numeric")
  
  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
  expect_equal(length(meta), as.numeric(stats["meta"]))
  expect_is(meta, "character")
  
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats,
                nrows = -1, skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  expect_is(body, "matrix")
  expect_equal(nrow(body), as.numeric(stats["variants"]))
  expect_equal(ncol(body), as.numeric(stats["columns"]))
  
  # Check for cariage return.
#  body[1,ncol(body)] <- "0/0:32,0:32:90:0,90,999\r" # Test example.
  expect_equal(grep("\r$", body[1,ncol(body)]), integer(0))
  
})



test_that("compiled vcfR_read_body works when file contains no variants",{

  data("vcfR_test")
  vcf2 <- vcfR_test[0,]
  
  myFile <- paste(test_dir, "myFile.vcf.gz", sep = "/")
  write.vcf(vcf2, myFile)
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', myFile)
#  meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', ex_file, stats, 0)
##  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats, nrows = -1, skip = 0, cols=1:stats['columns'], 0)
#
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', myFile, stats,
                nrows = 0, skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
#  
  expect_equal( ncol(body), ncol(vcf2@fix) + ncol(vcf2@gt) )
#
  expect_equal( nrow(body), nrow(vcf2@fix) )
  
})


test_that("compiled vcfR_read_body converts VCF missing to NA",{
  data("vcfR_test")
  vcfR_test@gt[1,3] <- "./."
  
  myFile <- paste(test_dir, "myFile.vcf.gz", sep = "/")
  write.vcf(vcfR_test, myFile)

  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', myFile)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
  expect_true( is.na(body[1,11]) )
  
  vcfR_test@gt[1,3] <- ".|."
  write.vcf(vcfR_test, myFile)

  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', myFile)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)
  unlink(myFile)
  
  expect_true( is.na(body[1,11]) )
#  body
  
  vcfR_test@gt[1,3] <- "./0"
  write.vcf(vcfR_test, myFile)
  
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', myFile)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', myFile, stats,
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
  
  myFile <- paste(test_dir, "myFile.vcf.gz", sep = "/")
  write.vcf(vcfR_test, myFile)

  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', myFile)
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', myFile, stats,
                nrows = stats['variants'], skip = 0, cols=1:stats['columns'],
                convertNA = 0, verbose = 0)
  unlink(myFile)

  expect_false( is.na(body[2,3]) )
  expect_false( is.na(body[1,11]) )
  expect_false( is.na(body[1,12]) )
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
  
#  write.vcf(vcf2, ex_file)
#  test <- read.vcfR(ex_file, verbose=FALSE)
#  unlink(ex_file)


#  expect_equal(ncol(test@fix), ncol(vcf2@fix))
#  expect_equal(ncol(test@gt), ncol(vcf2@gt))
#  expect_equal(nrow(test@fix), nrow(vcf2@fix))
#  expect_equal(nrow(test@gt), nrow(vcf2@gt))
})


test_that("read.vcfR works when file contains one variant",{
  data(vcfR_test)
  vcfR_test <- vcfR_test[1,]
  
  setwd( test_dir )
  
  write.vcf(vcfR_test, ex_file)
  vcfR_test2 <- read.vcfR(ex_file, verbose=FALSE)
  
  expect_equal( nrow(vcfR_test2@fix), 1)
  
  unlink(ex_file)
  setwd( original_dir )
})


test_that("read.vcfR works when file contains one variant, no meta",{
  data(vcfR_test)
  vcfR_test <- vcfR_test[1,]
  vcfR_test@meta <- vector(mode='character', length=0)

  setwd( test_dir )
  
  write.vcf(vcfR_test, ex_file)
  vcfR_test2 <- read.vcfR(ex_file, verbose=FALSE)
  
  expect_equal( nrow(vcfR_test2@fix), 1)
  
  unlink(ex_file)
  setwd( original_dir )
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
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', ex_file, stats,
                nrows = -1, skip = 0, cols=1:stats['columns'],
                convertNA = 1, verbose = 0)

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

