# create_chrom tests.

# detach(package:vcfR, unload=T)
library(vcfR)
context("windowing functions")

data(vcfR_example)

pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff, verbose=F)
pinf_mt <- proc_chrom(pinf_mt, verbose=FALSE)


gq <- extract.gt(pinf_vcf, element="GQ", as.numeric=TRUE)

gqw <- windowize_NM(gq, pos=pinf_mt@var.info$POS, 
                    starts=pinf_mt@win.info$start, 
                    ends=pinf_mt@win.info$end,
                    summary = "count")


test_that("gq is windowized", {
  expect_is(gqw, "matrix")
  expect_equal(sum(gqw[20,]), 5*ncol(gqw))  
})



#gqw[19:21,]

#gqw[1:5,]

#pinf_mt@var.info$POS[pinf_mt@var.info$POS >= 1 & pinf_mt@var.info$POS <= 1000]
#pinf_mt@var.info$POS[pinf_mt@var.info$POS >= 1001 & pinf_mt@var.info$POS <= 2000]
#pinf_mt@var.info$POS[pinf_mt@var.info$POS >= 2001 & pinf_mt@var.info$POS <= 3000]

#pinf_mt@var.info$POS[pinf_mt@var.info$POS >= 19001 & pinf_mt@var.info$POS <= 20000]


#nrow(gq[pinf_mt@var.info$POS >= 1 & pinf_mt@var.info$POS <= 1000,])
#nrow(gq[pinf_mt@var.info$POS >= 1001 & pinf_mt@var.info$POS <= 2000,])
#nrow(gq[pinf_mt@var.info$POS >= 2001 & pinf_mt@var.info$POS <= 3000,])
#nrow(gq[pinf_mt@var.info$POS >= 3001 & pinf_mt@var.info$POS <= 4000,])
#nrow(gq[pinf_mt@var.info$POS >= 4001 & pinf_mt@var.info$POS <= 5000,])
#nrow(gq[pinf_mt@var.info$POS >= 5001 & pinf_mt@var.info$POS <= 6000,])


#head(gqw[,1:10])



#gq[pinf_mt@var.info$POS >= 18001 & pinf_mt@var.info$POS < 19000,]
#gq[pinf_mt@var.info$POS >= 19001 & pinf_mt@var.info$POS < 20000,]
#gq[pinf_mt@var.info$POS >= 20001 & pinf_mt@var.info$POS < 21000,]


#apply(gq[pinf_mt@var.info$POS >= 20001 & pinf_mt@var.info$POS < 21000,], MARGIN=2, sum)



#cbind(pinf_mt@var.info$POS, gq)[118:130,1:4]



#head(gqw)
#nrow(gqw)


