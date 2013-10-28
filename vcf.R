# R functions to work with vcf format files.

# Read in a vcf format file.
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

vcf2gt <- function(x, cell = 1) {
  get.gt <- function(y) {unlist(lapply(strsplit(as.character(y), split = ":"), function(z) {z[cell]}))}
  apply(x[,10:ncol(x)], 2, get.gt)
}

vcf2gq <- function(x, cell = 3) {
  get.gq <- function(y) {unlist(lapply(strsplit(as.character(y), split = ":"), function(z) {z[cell]}))}
  apply(x[,10:ncol(x)], 2, get.gq)
}

