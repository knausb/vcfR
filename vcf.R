# R functions to work with vcf format files.

##### ##### Read in vcf format file ##### #####
# Gzipped files can be read in with:
# vcf <- read.vcf(gzfile(vcfsgz[i]))

read.vcf<-function(x){
  i <- -1 # Line counter.
  j <- 0 # Success?
  tmp <- scan(x, what="character", sep="\n", skip=0, nlines=1, quiet=T, comment.char="")
  if(length(grep('^##', tmp)) >0 ){
    while(j == 0){
      i <- i+1
      tmp <- scan(x, what="character", sep="\n", skip=i, nlines=1, quiet=T, comment.char="")
      if(length(grep('^##',tmp)) == 0){j <- 1}
    }
    read.table(x,header=T,sep='\t',skip=i,comment.char='')
  } else if (length(grep('^#',tmp)) >0 ){
    read.table(x,header=T,sep='\t',skip=0,comment.char='')
  } else {
    read.table(x,header=F,sep='\t',skip=0,comment.char='')
  }
}

##### ##### Extract genotypes from vcf data ##### #####

vcf2gt <- function(x, cell = 1) {
  get.gt <- function(y) {unlist(lapply(strsplit(as.character(y), split = ":"), function(z) {z[cell]}))}
  apply(x[,10:ncol(x)], 2, get.gt)
}

##### ##### Extract gene qualities from vcf data ##### #####

vcf2gq <- function(x, cell = 3) {
  get.gq <- function(y) {unlist(lapply(strsplit(as.character(y), split = ":"), function(z) {z[cell]}))}
  apply(x[,10:ncol(x)], 2, get.gq)
}

##### ##### Get Allele Frquency Spectrum ##### #####

get.af <- function (x) {
  af.sp1 <- vcf2gt(x)
  af.sp2 <- af.sp1[,26:34]
  af.sp1 <- af.sp1[,1:25]
  af.sp1 <- cbind(apply(af.sp1,MARGIN=1,sum.gt), apply(af.sp2,MARGIN=1,sum.gt))
  rownames(af.sp1) <- apply(x[,1:2], MARGIN=1, FUN=paste, collapse="_")
  af.sp1
}

##### ##### Remove the monomorphic SNPS ##### #####

rm.mono <- function (x){
  gt <- vcf2gt(x)
  gt <- x[apply(gt,FUN =function (x){length(unique(x))>1} ,MARGIN=1), ]
  if (nrow(gt) >= 1){
    return(gt)
  }
}



