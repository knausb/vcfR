#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("io functions")

#data(vcfR_example)

# Manage directories.
original_dir <- getwd()
test_dir <- tempdir()

# Load data
data(vcfR_example)

test_that("vcf file io works",{
  setwd(test_dir)

  write.vcf(pinf_vcf, "test.vcf")
  test <- read.vcf("test.vcf")

  unlink("test.vcf")
  setwd(original_dir)
  
  expect_is(test, "vcfR")
  expect_is(test@fix$POS, "integer")
  expect_is(test@fix$QUAL, "integer")
  expect_is(test@fix$QUAL, "numeric")
  
  expect_identical(names(test@fix)[1], "CHROM")
  
  
})


# Devel

setwd(test_dir)
data(vcfR_example)
write.vcf(pinf_vcf, "test.vcf")

x <- .Call('vcfR_read_to_line', PACKAGE = 'vcfR', "test.vcf")

setwd(original_dir)


# Devel2

original_dir <- getwd()
test_dir <- tempdir()
setwd(test_dir)
data(vcfR_example)
write.vcf(pinf_vcf, "test.vcf")

x1 <- .Call('vcfR_vcf_stats', PACKAGE = 'vcfR', "test.vcf")
x1
x2 <- .Call('vcfR_vcf_meta', PACKAGE = 'vcfR', "test.vcf", x1)
x2
x3 <- .Call('vcfR_vcf_body', PACKAGE = 'vcfR', "test.vcf", x1)
#x3
#x4 <- .Call('vcfR_vcf_body2', PACKAGE = 'vcfR', "test.vcf", x1)
x4 <- read.vcf.devel3("test.vcf")

unlink("test.vcf")
setwd(original_dir)


names(x3)[1]
class(x3$CHROM)
class(x3$POS)
class(x3$QUAL)

class(x4@fix$POS)
class(x4@fix$QUAL)







