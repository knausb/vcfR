

context("chromo function")

library(vcfR)
data("vcfR_example")


#chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
#chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
#chrom <- proc.chromR(chrom, verbose = FALSE)


##### ##### ##### ##### #####
# chromo, vcf only

chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)

chromo( chrom, boxp = FALSE )
chromo( chrom, boxp = TRUE )

##### ##### ##### ##### #####
# chromo, vcf and seq

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
chromo( chrom ) # Should error!
chrom <- proc.chromR(chrom, verbose = FALSE)

chromo( chrom, boxp = FALSE )
chromo( chrom, boxp = TRUE )


##### ##### ##### ##### #####
# chromo, vcf and annotation

chrom <- create.chromR(name="Supercontig", vcf=vcf, ann=gff, verbose=FALSE)
chromo( chrom )
chrom <- proc.chromR(chrom, verbose = FALSE)

chromo( chrom, boxp = FALSE )
chromo( chrom, boxp = TRUE )

##### ##### ##### ##### #####
# chromo, vcf, seq and annotation

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chromo( chrom )
chrom <- proc.chromR(chrom, verbose = FALSE)

chromo( chrom, boxp = FALSE )
chromo( chrom, boxp = TRUE )


##### ##### ##### ##### #####
# Create lists of dots and rectangles

rlst1 <- cbind( gff[ seq(1,23, by=2) , 4 ],
                0, 
                gff[ seq(1,23, by=2) , 5 ],
                500
              )
rlist <- list(rlst1)

myList1 <- list(title = "Track1",
                dmat  = chrom@var.info[,2:4],
                rlst = rlst1,
                rcol=4,
                bwcol=1:2
                )
# rlist, dcol, rcol, rbcol and bwcol.)
#names(myList1)


##### ##### ##### ##### #####
# chromo, vcf and drlist1


chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)

chromo( chrom, boxp = FALSE , chrom.e = chrom@len, drlist1 = myList1 )

chrom <- proc.chromR(chrom, verbose = FALSE)

chromo( chrom, boxp = FALSE , chrom.e = chrom@len, drlist1 = myList1 )
chromo( chrom, boxp = TRUE , chrom.e = chrom@len, drlist1 = myList1  )

chromo( chrom, boxp = TRUE , chrom.e = chrom@len, drlist1 = myList1, drlist2 = myList1  )
chromo( chrom, boxp = TRUE , chrom.e = chrom@len, drlist1 = myList1, drlist2 = myList1, drlist3 = myList1  )

##### ##### ##### ##### #####
# xlim

chromo( chrom, boxp = TRUE , chrom.e = chrom@len, drlist1 = myList1, drlist2 = myList1, xlim=c(2e4, 5e4) )



##### ##### ##### ##### #####
# EOF.