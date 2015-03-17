# detach(package:vcfR, unload=T)
library(vcfR)
context("proc_chrom functions")

data(vcfR_example)



pinf_mt <- create_chrom(name='pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=FALSE)
pinf_mt <- proc_chrom(pinf_mt, win.size=1000, verbose=FALSE)

test_that("Created Chrom", {
  expect_is(pinf_mt, "Chrom")
})

pinf_mt2 <- proc_chrom(pinf_mt, win.size=1e4, verbose=FALSE)

test_that("proc_chrom creates different window sizes", {
  expect_equal(length(pinf_mt@win.info$end), 40)
  expect_equal(length(pinf_mt2@win.info$end), 4)
})

# Binaries
#wins <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(pinf_dna))
#win2 <- .Call('vcfR_windowize_fasta', PACKAGE = 'vcfR', wins=wins, seq=as.character(pinf_dna)[1,])






#head(pinf_mt)
#pinf_mt
#names(pinf_mt)
#plot(pinf_mt)
#pinf_mt <- masker(pinf_mt)





# proc_chrom2
#pinf_mt2 <- proc_chrom2(pinf_mt, win.size=1000)

#chromoqc(pinf_mt2)



# Variants, sequence and annotations.
#expect_that(pinf_mt, is_a("Chrom"))



# Variants and sequence.

#pinf_mt <- proc_chrom(pinf_mt, win.size=1000, verbose=F)


#expect_that(pinf_mt, is_a("Chrom"))



# Variants and annotations.
#pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, ann=pinf_gff, verbose=F)
#pinf_mt <- masker(pinf_mt)
#pinf_mt <- proc_chrom(pinf_mt, win.size=1000)
#expect_that(pinf_mt, is_a("Chrom"))

#pinf_mt <- create_chrom('pinf_mt', vcf=pinf_vcf, verbose=F)
#expect_that(pinf_mt, is_a("Chrom"))


#print("proc_chrom functions finished")
