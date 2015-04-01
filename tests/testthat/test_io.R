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
#  expect_is(test@fix$QUAL, "numeric")
  
  expect_identical(names(test@fix)[1], "CHROM")
  
  
})


# Devel

#setwd(test_dir)
#data(vcfR_example)
#write.vcf(pinf_vcf, "test.vcf")

#x <- .Call('vcfR_read_to_line', PACKAGE = 'vcfR', "test.vcf")

#setwd(original_dir)


# Devel2

#original_dir <- getwd()
#test_dir <- tempdir()
#setwd(test_dir)
#data(vcfR_example)
#write.vcf(pinf_vcf, "test.vcf")


#x1 <- .Call('vcfR_vcf_stats', PACKAGE = 'vcfR', "test.vcf")
#x1
#x2 <- .Call('vcfR_vcf_meta', PACKAGE = 'vcfR', "test.vcf", x1)
#x2
#x3 <- .Call('vcfR_vcf_body', PACKAGE = 'vcfR', "test.vcf", x1)
#x3
#x4 <- .Call('vcfR_vcf_body2', PACKAGE = 'vcfR', "test.vcf", x1)
#x4 <- read.vcf.devel3("test.vcf", limit=1e2)
#x4 <- read.vcf.devel3("test.vcf", limit=1e9)


#unlink("test.vcf")
#setwd(original_dir)

#print(object.size(x4), units='Mb')


#names(x3)[1]
#class(x3$CHROM)
#class(x3$POS)
#class(x3$QUAL)

#class(x4@fix$POS)
#class(x4@fix$QUAL)


#vcf <- "/home/likewise-open/USDA-ARS/knausb/gits/vcf_data/1kgenomes/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"

#x7 <- read.vcf.devel3(vcf, limit=1e9)


#x5 <- .Call('vcfR_vcf_stats', PACKAGE = 'vcfR', vcf)

#x6 <- .Call('vcfR_vcf_body', PACKAGE = 'vcfR', vcf, x5)

#x7 <- .Call('vcfR_ram_test', PACKAGE = 'vcfR')

#print(object.size(x7), units='b')
#print(object.size(x7), units='Mb')

#print(object.size(.Call('vcfR_ram_test', PACKAGE = 'vcfR', nrow=1e4, ncol=2000)), units="Gb")
#print(object.size(.Call('vcfR_ram_test', PACKAGE = 'vcfR', nrow=3e5, ncol=2000)), units="Gb")

#write_vcf_body( fix = pinf_vcf@fix, gt = pinf_vcf@gt, filename = "test.vcf", mask = 0L)



#.Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = pinf_vcf@fix, gt = pinf_vcf@gt, filename = "test.vcf", mask = 0)

#write.vcf.devel3(pinf_vcf, file = "test.vcf", mask = FALSE, APPEND = FALSE)



