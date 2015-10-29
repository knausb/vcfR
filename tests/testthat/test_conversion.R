

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

#vcf <- read.vcf(ex_file, verbose=FALSE)
#dna <- ape::read.dna(seq_file, format = "fasta")

#gene <- dna[,159730:160890]
#pos <- as.numeric(vcf@fix[,'POS'])

#vcf.gene <- vcf
#vcf.gene@gt <- vcf.gene@gt[pos >= 159730 & pos <= 160890,]
#vcf.gene@fix <- vcf.gene@fix[pos >= 159730 & pos <= 160890,]


##### ##### ##### ##### #####

test_that("vcfR2DNAbin with consensus works",{
  my_DNAbin <- vcfR2DNAbin( vcf, consensus = TRUE, gt.split = "|" )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], ncol(vcf@gt) - 1 )
  expect_equal( dim(my_DNAbin)[2], nrow(extract.indels(vcf)@gt) )
})


test_that("vcfR2DNAbin with extract.haps works",{
  my_DNAbin <- vcfR2DNAbin( vcf, consensus = FALSE, extract.haps = TRUE, gt.split = "|", verbose = FALSE )
  
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], 2*(ncol(vcf@gt) - 1) )
  expect_equal( dim(my_DNAbin)[2], nrow(extract.indels(vcf)@gt) )
})


test_that("vcfR2DANbin with extract.haps works",{
  my_DNAbin <- vcfR2DNAbin( vcf, consensus = FALSE, extract.haps = TRUE, 
                            gt.split = "|", verbose = FALSE,
                            ref.seq = gene, start.pos = gff[1,4] )
  
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], 2*(ncol(vcf@gt) - 1) )
  expect_equal( dim(my_DNAbin)[2], dim(gene)[2] )
})


#debug(vcfR2DNAbin)

##### ##### ##### ##### #####
# EOF.