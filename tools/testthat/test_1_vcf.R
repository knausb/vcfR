
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
# Input functions
#
##### ##### ##### ##### #####














##### ##### ##### ##### #####
#
# R version
#
##### ##### ##### ##### #####























##### ##### ##### ##### #####



  
#debug(read.vcfR)
#unlink(ex_file)


##### ##### ##### ##### #####
#
# IO works for VCF data that 
# contains no gt region.
#
##### ##### ##### ##### #####






##### ##### ##### ##### #####
# Manage directories.

setwd(original_dir)

##### ##### ##### ##### #####
# EOF.

