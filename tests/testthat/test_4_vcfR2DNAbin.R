
#library(testthat)
library(vcfR)
context("vcfR2DNAbin functions")

data(vcfR_example)
gene <- dna[1,gff[1,4]:gff[1,5]]

##### ##### ##### ##### #####
#
# Missing data handling
#
##### ##### ##### ##### #####


test_that("Works with no variants",{
  vcf <- vcf[0,]
  my_DNAbin <- vcfR2DNAbin( vcf )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal(dim(my_DNAbin)[1], ncol(vcf@gt) - 1 )
})

test_that("Works with only indels, no SNPs",{
  vcf <- extract.indels( vcf, return.indels = TRUE)
  my_DNAbin <- vcfR2DNAbin( vcf )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal(dim(my_DNAbin)[1], ncol(vcf@gt) - 1 )
})

test_that("Works with no variants, ref.seq is not NULL",{
  vcf <- vcf[0,]
  my_DNAbin <- vcfR2DNAbin( vcf, ref.seq = gene, start.pos = gff[1,4], verbose = FALSE )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], ncol(vcf@gt) - 1 )
  expect_equal( dim(my_DNAbin)[2], length(gene) )
  expect_equal( length(ape::seg.sites(my_DNAbin)), 0 )
  expect_equal( sum(ape::base.freq(my_DNAbin) > 0), 4 )
})



##### ##### ##### ##### #####
#
# Fabricate data which varies in ploidy
#
##### ##### ##### ##### #####

data(vcfR_example)

# Fabricate a vcf of haploid data.
#haps <- extract.haps(vcf)
haps <- extract.gt(vcf, return.alleles=FALSE, allele.sep="|")
haps <- apply(haps, MARGIN=2, function(x){unlist(lapply(strsplit(x, split="|"), function(x){x[1]}))})
gt2 <- extract.gt(vcf, extract = FALSE)
gt2 <- matrix( paste( haps, gt2, sep=":"), nrow=nrow(haps), dimnames=dimnames(haps) )
is.na(gt2[is.na(haps)]) <- TRUE
vcf1 <- vcf
vcf1@gt <- cbind(vcf@gt[,'FORMAT'], gt2)
colnames(vcf1@gt)[1] <- 'FORMAT'
rm(haps)
rm(gt2)


# Fabricate a vcf of triploid data.
#haps <- extract.haps(vcf)
haps <- extract.gt(vcf, return.alleles=FALSE, allele.sep="|")
haps2 <- apply(haps, MARGIN=2, function(x){unlist(lapply(strsplit(x, split="|"), function(x){x[1]}))})
haps <- matrix( paste( haps, haps2, sep="|"), nrow=nrow(haps), dimnames=dimnames(haps) )
is.na(haps[is.na(haps2)]) <- TRUE
gt2 <- extract.gt(vcf, extract = FALSE)
gt2 <- matrix( paste( haps, gt2, sep=":"), nrow=nrow(haps), dimnames=dimnames(haps) )
is.na(gt2[is.na(haps)]) <- TRUE
vcf3 <- vcf
vcf3@gt <- cbind(vcf@gt[,'FORMAT'], gt2)
colnames(vcf3@gt)[1] <- 'FORMAT'
rm(haps)
rm(gt2)


##### ##### ##### ##### #####
#
# Test haploid data
#
##### ##### ##### ##### #####

test_that("vcfR2DNAbin works for haploid data, no ref.seq",{
  my_DNAbin <- vcfR2DNAbin( vcf1, gt.split = "|" )
  expect_equal( dim(my_DNAbin)[2], nrow( extract.indels(vcf)@gt ) )
})




##### ##### ##### ##### #####
#
# Test diploid data
#
##### ##### ##### ##### #####


test_that("vcfR2DNAbin works for diploid data, no ref.seq",{
  my_DNAbin <- vcfR2DNAbin( vcf, gt.split = "|", verbose = FALSE )
  expect_true( inherits(my_DNAbin, "DNAbin") )  
  expect_equal( dim(my_DNAbin)[2], nrow( extract.indels(vcf)@gt ) )
  expect_equal( dim(my_DNAbin)[1], 2 * (ncol(vcf@gt) - 1) )
})


test_that("vcfR2DNAbin works for diploid data, no ref.seq, wrong gt.split",{
  
  

})

test_that("vcfR2DNAbin with consensus works",{
  my_DNAbin <- vcfR2DNAbin( vcf, consensus = TRUE, extract.haps = FALSE, gt.split = "|" )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], ncol(vcf@gt) - 1 )
  expect_equal( dim(my_DNAbin)[2], nrow(extract.indels(vcf)@gt) )
})


test_that("vcfR2DNAbin works for diploid data, ref.seq is not NULL",{
  my_DNAbin <- vcfR2DNAbin( vcf, gt.split = "|", ref.seq = gene, start.pos = gff[1,4], verbose = FALSE )
  expect_true( inherits(my_DNAbin, "DNAbin") )
  expect_equal( dim(my_DNAbin)[1], 2 * (ncol(vcf@gt) - 1) )
  expect_equal( dim(my_DNAbin)[2], dim(gene)[2] )
  expect_true( length(ape::seg.sites(my_DNAbin)) > 0 )
#  expect_equal( dim(my_DNAbin)[2], nrow( extract.indels(vcf)@gt ) )
})



##### ##### ##### ##### #####
#
# Test triploid data
#
##### ##### ##### ##### #####


test_that("vcfR2DNAbin works for triploid data, no ref.seq",{
  my_DNAbin <- vcfR2DNAbin( vcf3, gt.split = "|", verbose = FALSE )
  expect_true( inherits(my_DNAbin, "DNAbin") )  
  expect_equal( dim(my_DNAbin)[2], nrow( extract.indels(vcf)@gt ) )
  expect_equal( dim(my_DNAbin)[1], 3 * (ncol(vcf@gt) - 1) )
})





##### ##### ##### ##### #####
# EOF.