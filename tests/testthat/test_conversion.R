

#library(testthat)
#detach(package:vcfR, unload=TRUE)
library(vcfR)
context("conversion functions")

##### ##### ##### ##### #####
# Example data.

data(vcfR_example)

gene <- dna[,gff[1,4]:gff[1,5]]



#ex_file <- system.file("extdata", "pinf_sc1_100_sub.vcf.gz", package = "vcfR")
#seq_file <- system.file("extdata", "pinf_sc100.fasta", package = "vcfR")

#vcf <- read.vcfR(ex_file, verbose=FALSE)
#dna <- ape::read.dna(seq_file, format = "fasta")

#gene <- dna[,159730:160890]
#pos <- as.numeric(vcf@fix[,'POS'])

#vcf.gene <- vcf
#vcf.gene@gt <- vcf.gene@gt[pos >= 159730 & pos <= 160890,]
#vcf.gene@fix <- vcf.gene@fix[pos >= 159730 & pos <= 160890,]


##### ##### ##### ##### #####

test_that("vcfR2genlight works",{
  suppressWarnings( gl <- vcfR2genlight(vcf) )
  ma <- t(as.matrix(gl))
  
  vcf2 <- vcf[is.biallelic(vcf),]
  gt <- extract.gt(vcf2)
  gt[gt=='0|0'] <- 0
  gt[gt=='0|1'] <- 1
  gt[gt=='1|0'] <- 1
  gt[gt=='1|1'] <- 2

  expect_equal(sum(ma == gt, na.rm=TRUE), sum(!is.na(gt)))
})


##### ##### ##### ##### #####



#debug(vcfR2DNAbin)

##### ##### ##### ##### #####
# EOF.