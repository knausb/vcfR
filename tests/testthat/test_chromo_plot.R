

library(vcfR)
data("vcfR_example")
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = FALSE)




##### ##### ##### ##### #####
# dot.plot

head(chrom@var.info)

par(oma=c(1,1,1,1))
dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom", hline=seq(0,2000,by=500), mwidth=6, layout=T )
par(oma=c(0,0,0,0))


##### ##### ##### ##### #####
# bar.plot

head(chrom@win.info)

wmat <- as.matrix( cbind(chrom@win.info[,2], 
                         chrom@win.info$A + chrom@win.info$T, 
                         chrom@win.info$C + chrom@win.info$G
                         ))
head(wmat)

bar.plot( wmat, title="MyBar" )
bar.plot( wmat, title="MyBar", scale=T )

##### ##### ##### ##### #####

mwidth <- 4
layout( matrix( 1:4, nrow=2, ncol=2, byrow = TRUE ), widths = c(mwidth,1) )
dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom", hline=seq(0,2000,by=500), mwidth=6, layout=F )
bar.plot( wmat, title="MyBar", layout = FALSE )

##### ##### ##### ##### #####
# EOF.