

library(vcfR)
data("vcfR_example")
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = FALSE)




##### ##### ##### ##### #####
# dot.plot

head(chrom@var.info)

par(oma=c(1,1,1,1))
dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom", hline=seq(0,2000,by=500), mwidth=6, layout=T )

dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom",
          hline=seq(0,2000,by=500), 
          col=c( rgb(34,139,34,10,maxColorValue=255), 
                 rgb(64,224,208,10,maxColorValue=255) ),
          mwidth=6, layout=T )

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
# rect.plot

mylist <- list(as.matrix(chrom@ann[,4:5]))
mylist[[2]] <- chrom@ann[seq(from=0, to=nrow(chrom@ann), by=2),4:5]

rect.plot( lst = mylist, xmax=1e5, heights = c(0.5,1))

rect.plot(chrom@seq.info, xmax=1e5, heights = c(1, 0.5), col=c('green', 'red'))



##### ##### ##### ##### #####
# Combo plot


mwidth <- 8
layout( matrix( 1:4, nrow=2, ncol=2, byrow = TRUE ), widths = c(mwidth,1) )
dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom", hline=seq(0,2000,by=500), mwidth=6, layout=F )
bar.plot( wmat, title="MyBar", layout = FALSE )


mwidth <- 5
layout( matrix( 1:6, nrow=3, ncol=2, byrow = TRUE ), widths = c(mwidth,1), heights = c(1,1,0.25) )
par(oma=c(2.1,0.1,0.1,0.1))
dot.plot( as.matrix( chrom@var.info[,c(2,3,4)] ), title="MyChrom", hline=seq(0,2000,by=500), mwidth=6, layout=F )
bar.plot( wmat, title="MyBar", layout = FALSE, scale = FALSE )
rect.plot(chrom@seq.info, xmax=1e5, title="Nucleotides", heights = c(1, 0.5), col=c('green', 'red'))
axis(side=1)
null.plot()
par(oma=c(0,0,0,0))


##### ##### ##### ##### #####
# chromoqc

chromoqc(chrom)




##### ##### ##### ##### #####
# EOF.