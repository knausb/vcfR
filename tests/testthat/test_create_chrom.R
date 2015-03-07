# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("create_chrom functions")

data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
expect_that(pinf_mt, is_a("Chrom"))

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, verbose=F)
expect_that(pinf_mt, is_a("Chrom"))

pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, ann=pinf_gff, verbose=F)
expect_that(pinf_mt, is_a("Chrom"))

pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, verbose=F)
expect_that(pinf_mt, is_a("Chrom"))


