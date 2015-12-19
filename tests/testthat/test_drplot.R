

context("dr.plot functions")


library(vcfR)
data("vcfR_example")
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)

chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
chrom <- proc.chromR(chrom, verbose = FALSE)


##### ##### ##### ##### #####
# dr.plot

dr.plot( chrom.e = 100001 )


head(chrom@win.info)

mat1 <- cbind(chrom@win.info[,'start'],
              0, 
              chrom@win.info[,'end'], 
              chrom@win.info[,'variants']
              )

mat2 <- cbind(chrom@win.info[,'start'],
              0, 
              chrom@win.info[,'end'],
              rowSums(chrom@win.info[,c('A', 'T')])
)

mat3 <- cbind(chrom@win.info[,'start'],
              rowSums(chrom@win.info[,c('A', 'T')]),
              chrom@win.info[,'end'],
              rowSums(chrom@win.info[,c('A', 'C', 'G', 'T')])
)

dr.plot( rlst = mat1, chrom.e = 100001 )

dr.plot( rlst = list(mat1), chrom.e = 100001 )

dr.plot( rlst = list(mat2, mat3), chrom.e = 100001 )

dr.plot( rlst = list(mat2, mat3), chrom.e = 100001, xaxt = "n", yaxt = "n") #, frame.plot = TRUE )


dr.plot( dmat = as.matrix( chrom@var.info[,c(2,3)] ), chrom.e = 100001, 
         dcol=c(rgb(34,139,34, maxColorValue = 255), rgb(0,206,209, maxColorValue = 255)) )

dr.plot( dmat = as.matrix( chrom@var.info[chrom@var.info[,3] > 59,c(2,3)] ), chrom.e = 100001, 
         dcol=c(rgb(34,139,34, maxColorValue = 255), rgb(0,206,209, maxColorValue = 255)) )


dr.plot( dmat = as.matrix( chrom@var.info[,c(2,3,4)] ), chrom.e = 100001, 
         dcol=c(rgb(34,139,34, maxColorValue = 255), rgb(0,206,209, maxColorValue = 255)) )

dr.plot( dmat = as.matrix( chrom@var.info[,c(2,3,4)] ),
         rlst = list(mat2, mat3), chrom.e = 100001,
         dcol=c(rgb(34,139,34, maxColorValue = 255), rgb(0,206,209, maxColorValue = 255)) )
axis(side=1)

dr.plot( dmat = as.matrix( chrom@var.info[,c(2,3,4)] ),
         rlst = list(mat2, mat3), chrom.e = 100001,
         dcol=c(rgb(34,139,34,100, maxColorValue = 255), rgb(0,206,209,100, maxColorValue = 255)),
         rbcol = c("#000000", "#000000"),
         xlim=c(1e4,3e4), title="MyPlot",
         hline = seq(0, 5e3, by=5e2) )
axis(side=1)




