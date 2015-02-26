# chromR.

#### Class definition. ####

setOldClass("DNAbin")

#' @title Chrom class
#'
#' @name Chrom-class
#' @rdname Chrom-class
#'
#' @description
#' A class for representing chromosomes (or contigs).
#'
#' @details Defines a class for chromosomal or contig data.
#' 
#' This object has quite a few slots.
#' 
#' \itemize{
#'   \item name name of the object  
#'   \item seq object of class DNAbin (ape)
#'   \item len length of sequence (integer)
#'   \item vcf.meta vcf meta data
#'   \item vcf.fix vcf fixed data
#'   \item vcf.gt vcf genotype data
#'   \item ann annotation data in a gff-like data.frame
#'
#'   \item var.info a data.frame containing information on variants
#'   \item win.info a data.frame containing information on windows
#'   \item seq.info a list containing information on the sequence
#'      
#   \item pop1 vector indicating members of pop1
#   \item pop2 vector indicating members of pop2
#'   
#   \item acgt.w matrix indicating range of chromosome # of nucleotide compositions
#   \item n.w matrix indicating locations of blocks of Ns in chromosome
#'   
#   \item windows matrix of windows
#   \item nuccomp.w data.frame of nucleotide composition windows
#   \item snpden.w data.frame of snp density windows
#'   
#'   \item gt.m matrix of genotypes
#'   \item sfs matrix for the site frequency spectrum
#'   \item link matrix for linkages
#'   
#   \item mask a logical vector to indicate masked variants
#' }
#' 
#' More descriptions can be put here.
#' 
#' @export
#' @import methods
setClass(
  Class="Chrom",
  representation=representation(
    name = "character",
    seq = "DNAbin",
    len = "integer",
    vcf.meta = "character",
    vcf.fix = "data.frame",
    vcf.gt = "data.frame",
    ann = "data.frame",
    #
    var.info = "data.frame",
    win.info = "data.frame",
    seq.info = "list",
    #
    #    pop1 = "vector",
    #    pop2 = "vector",
    #
    #    acgt.w = "matrix",
    #    n.w = "matrix",
    #
    #    windows = "matrix",
    #    nuccomp.w = "data.frame",
    #    snpden.w = "data.frame",
    #
    gt.m = "matrix",
    vcf.stat = "data.frame",
    sfs = "matrix",
    link = "matrix",
    #
    mask = "logical"
  ),
  prototype=prototype(
    vcf.fix = data.frame(matrix(ncol=8, nrow=0,
                                dimnames=list(c(),
                                              c('chrom','pos','id','ref','alt','qual','filter','info'))
    ),
    stringsAsFactors=FALSE),
    vcf.stat = data.frame(matrix(ncol=11, nrow=0, 
                                 dimnames=list(c(),
                                               c('Allele_num','R_num','A_num','Ho','He','Ne','theta_pi','theta_w','theta_h','tajimas_d','fw_h'))
    ),
    stringsAsFactors=FALSE),
    ann = data.frame(matrix(ncol=9, nrow=0,
                            dimnames=list(c(),
                                          c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
    ),
    stringsAsFactors=FALSE)
  )
  #  )
)

##### ##### Generic methods. #####

setMethod(
  f="show",
  signature = "Chrom",
  definition=function(object){
    message("*** Class Chrom, method Show *** \n")
    message(paste("Name: ", object@name, "\n"))
    message(paste("Length: ", object@len, "\n"))
    message("Use head(object) for more details.\n")
    message("******* End Show (Chrom) ******* \n")
  }
)

setMethod(
  f="print",
  signature="Chrom",
  definition=function (x,y,...){
    message("***** Object of class 'Chrom' *****\n")
    message(paste("Name: ", x@name, "\n"))
    message(paste("Length: ", x@len, "\n"))
    message("\nVCF fixed data:\n")
    message("Last column (info) omitted.\n")
    message("\nVCF variable data:\n")
    message(paste("Columns: ", ncol(x@vcf.gt), "\n"))
    message(paste("Rows: ", nrow(x@vcf.gt), "\n"))
    message("(First column is format.)\n")
    message("\nAnnotation data:\n")
    if(length(x@ann)>0){
      print(head(x@ann[,1:8], n=4))
      message("Last column (attributes) omitted.\n")
    } else {
      message("Empty slot.\n")
    }
    message("***** End print (Chrom) ***** \n")
  }
)







##### ##### Set populations #####

#' @rdname Chrom-methods
#' @export
#' @aliases set.pop1 
#' 
#' @param pop1 a numeric vector indicating the samples in population 1
#' 
set.pop1 <- function(x, pop1){
  x@pop1 <- pop1
  return(x)  
}

#' @rdname Chrom-methods
#' @export
#' @aliases set.pop2
#' 
#' @param pop2 a numeric vector indicating the samples in population 2
#' 
set.pop2 <- function(x, pop2){
  x@pop2 <- pop2
  return(x)  
}


##### ##### seq.info functions #####

#acgt.win <- function(x, max.win=1000, regex="[acgtwsmkrybdhv]"){
regex.win <- function(x, max.win=1000, regex="[acgtwsmkrybdhv]"){
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  if(is.matrix(as.character(x@seq))){
    seq <- as.character(x@seq)[1:length(x@seq)]    
  }
  if(is.list(as.character(x@seq))){
    seq <- as.character(x@seq)[[1]]
  }
  # Subset to nucleotides of interest.
  seq <- grep(regex, seq, ignore.case=T, perl=TRUE)
  if(length(seq) == 0){
    return(matrix(NA, ncol=2))
    break
  }
  #
  bp.windows <- matrix(NA, ncol=2, nrow=max.win)
  bp.windows[1,1] <- seq[1]
  i <- 1
  # Scroll through the sequence looking for 
  # gaps (nucledotides not in the regex).
  # When you find them make a window.
  # Sequences with no gaps will have no
  # windows.
  for(j in 2:length(seq)){
    if(seq[j]-seq[j-1] > 1){
      bp.windows[i,2] <- seq[j-1]
      i <- i+1
      bp.windows[i,1] <- seq[j]
    }
  }
  bp.windows[i,2] <- seq[j]
  if(i == 1){
    # If there is one row we get an integer.
    # We need a matrix.
    bp.windows <- bp.windows[1:i,]
    bp.windows <- matrix(bp.windows, ncol=2)
  } else {
    bp.windows <- bp.windows[1:i,]
  }
  #  x@acgt.w <- bp.windows
  #  return(x)
  return(bp.windows)
}

##### ##### win.info functions #####

#' @rdname Chrom-methods
#' @export
#' @aliases windowize
#'
# @description
# Creates windows
#'
#' @param win.size window size, in base pairs
#' @param max.win maximum window size
#'
#' @details
#' Reads in a vcf file and stores it in a Chrom class.
#' 
#'
windowize <- function(x, win.size=1000, max.win=10000){
  #  acgt.w <- x@acgt.w
  acgt.w <- x@seq.info$nuc.win
  windows <- matrix(NA, ncol=2, nrow=max.win)
  i <- 1
  for(j in 1:nrow(acgt.w)){
    while(acgt.w[j,2]-acgt.w[j,1] > win.size){
      windows[i,1] <- acgt.w[j,1]
      windows[i,2] <- acgt.w[j,1] + win.size - 1
      acgt.w[j,1] <- acgt.w[j,1] + win.size + 0
      i <- i+1
      if(i > max.win){
        print(paste("max i equals", max.win))
        print(paste("i equals", i))
        print(paste("j equals", j))
        message("chrom.r error: max.win is too small.\n")
        break
      }
    }
    windows[i,1] <- acgt.w[j,1]
    windows[i,2] <- acgt.w[j,2]
    i <- i+1
  }
  x@windows <- windows[1:i-1,]
  return(x)
}

gc.win <- function(x){
  win <- matrix(ncol=7,
                nrow=nrow(x@windows),
                dimnames=list(c(),
                              c('index','start','stop','gc','at','gcf','atf'))
  )
  win[,1] <- 1:nrow(win)
  win[,2] <- x@windows[,1]
  win[,3] <- x@windows[,2]
  chrom <- as.character(x@seq)[[1]]
  #
  count.nucs <- function(x){
    chrom <- chrom[x[2]:x[3]]
    win[x[1],4] <<- length(grep("[GgCc]", chrom, perl=TRUE))
    win[x[1],5] <<- length(grep("[AaTt]", chrom, perl=TRUE))
  }
  #
  apply(win, MARGIN=1, count.nucs)
  win[,6] <- win[,4]/(win[,3]-win[,2])
  win[,7] <- win[,5]/(win[,3]-win[,2])
  x@nuccomp.w <- as.data.frame(win)
  return(x)
}

snp.win <- function(x){
  snp <- matrix(ncol=5,
                nrow=nrow(x@windows),
                dimnames=list(c(),
                              c('index','start','stop','count','density'))
  )
  snp[,1] <- 1:nrow(snp)
  snp[,2] <- x@windows[,1]
  snp[,3] <- x@windows[,2]
  vcf <- x@vcf.fix$POS[x@var.info$mask]
  #
  count.snps <- function(x){
    vcf2 <- vcf[vcf >= x[2] & vcf <= x[3]]
    snp[x[1],4] <<- length(vcf2)
  }
  apply(snp, MARGIN=1, count.snps)
  snp[,5] <- snp[,4]/(snp[,3]-snp[,2]+1)
  x@snpden.w <- as.data.frame(snp)
  return(x)
}

var.win <- function(x, win.size=1000){
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  if(is.matrix(as.character(x@seq))){
    seq <- as.character(x@seq)[1:length(x@seq)]    
  }
  if(is.list(as.character(x@seq))){
    seq <- as.character(x@seq)[[1]]
  }
  #
  genic_sites <- rep(0, times=x@len)
  genic_sites[unlist(apply(x@ann[, 4:5], MARGIN=1, function(x){seq(from=x[1], to=x[2], by=1)}))] <- 1
  #
  win.info <- seq(1,x@len, by=win.size)
  win.info <- cbind(win.info, c(win.info[-1]-1, x@len))
  win.info <- cbind(1:nrow(win.info), win.info)
  win.info <- cbind(win.info, win.info[,3]-win.info[,2]+1)
  #  win.info <- cbind(win.info, matrix(ncol=7, nrow=nrow(win.info)))
  #
  win.proc <- function(y, seq){
    seq <- seq[y[2]:y[3]]
    a <- length(grep("[aA]", seq, perl=TRUE))
    c <- length(grep("[cC]", seq, perl=TRUE))
    g <- length(grep("[gG]", seq, perl=TRUE))
    t <- length(grep("[tT]", seq, perl=TRUE))
    n <- length(grep("[nN]", seq, perl=TRUE))
    o <- length(grep("[^aAcCgGtTnN]", seq, perl=TRUE))
    count <- sum(x@vcf.fix$POS[x@var.info$mask] >= y[2] & x@vcf.fix$POS[x@var.info$mask] <= y[3])
    genic <- sum(genic_sites[y[2]:y[3]])
    #
    c(a,c,g,t,n,o, count, genic)
  }
  #
  win.info <- cbind(win.info, t(apply(win.info, MARGIN=1, win.proc, seq=seq)))
  win.info <- as.data.frame(win.info)
  names(win.info) <- c('window','start','end','length','A','C','G','T','N','other','variants', 'genic')
  win.info
}

##### ##### vcf functions #####

vcf.fix2gt.m <- function(x){
  snames <- names(x@vcf.gt)[-1]
  pos <- paste(x@vcf.fix[,1], x@vcf.fix[,2], sep="_")
  #  pos <- x@vcf.fix[,2]
  x1 <- as.matrix(x@vcf.gt)
  nsamp <- ncol(x1) - 1
  #
  x1 <- cbind(unlist(lapply(strsplit(x1[,1], ":"), function(x){grep("GT", x)})),x1)
  #
  get.gt <- function(x){
    cell <- as.numeric(x[1])
    x <- lapply(strsplit(x[3:length(x)], ":"), function(x){x[cell]})
    unlist(x)
  }
  x1 <- apply(x1, MARGIN=1, get.gt)
  x1[x1=="0/0"] <- 0
  x1[x1=="0/1"] <- 1
  x1[x1=="1/0"] <- 1
  x1[x1=="1/1"] <- 2
  x1 <- as.numeric(x1)
  x1 <- matrix(data=x1, ncol=nsamp, byrow=TRUE,
               dimnames=list(pos, snames))
  x@gt.m <- x1
  return(x)
}

##### ##### gt.m2sfs #####

gt.m2sfs <- function(x){
  #  cat(x@pop1)
  #  cat(length(x@pop1))
  #  cat('\n')
  #  if(length(x@pop1) < 1 | length(x@pop2) < 1 | is.na(x@pop1) | is.na(x@pop2)){
  #    cat("One or both populations are not defined\n")
  #    cat("Creating arbitrary populations\n")
  #    x@pop1 <- 1:floor(ncol(x@vcf.gt[,-1])/2)
  #    x@pop2 <- c(1+max(1:floor(ncol(x@vcf.gt[,-1])/2))):ncol(x@vcf.gt)
  #  }
  pop1 <- x@gt.m[x@mask, x@pop1]
  pop2 <- x@gt.m[x@mask, x@pop2]
  sfs <- matrix(ncol=ncol(pop1)*2+1, nrow=ncol(pop2)*2+1)
  sfs1d <- cbind(rowSums(pop2)+1, rowSums(pop1)+1)
  sfs1d[,1] <- nrow(sfs) + 1 - sfs1d[,1]
  apply(sfs1d, MARGIN=1, function(x){
    if(is.na(sfs[x[1],x[2]])){
      sfs[x[1],x[2]] <<- 1
    }else{
      sfs[x[1],x[2]] <<- sfs[x[1],x[2]] +1
    }}
  )
  x@sfs <- sfs
  return(x)
}

gt2popsum <- function(x){
  if(class(x) != "Chrom"){stop("Object is not of class Chrom")}
  #  stopifnot(class(x) == "Chrom")
  #  gt <- extract.gt(x, element = "GT", mask = x@var.info$mask)
  #  stopifnot(length(grep("(1/1|0/0|0/1)", unique(as.vector(gt)))) == 3)
  #  gt <- x@gt.m
  #
  hwe <- function(x){
    # Genotype counts
    n11 <- x[1]
    n1i <- x[2]
    nii <- x[3]
    n   <- sum(n11, n1i, nii)
    #
    # Allele count and frequency
    n1 <- 2*n11 + n1i
    p1 <- n1/(2*n)
    #
    # Probability
    num <- (factorial(n) * factorial(n1) * factorial(2*n - n1) * 2^(n1i))
    den <- (factorial((n1 - n1i)/2) * factorial(n1i) * factorial(n-(n1+n1i)/2) * factorial(2*n))
    prob <- num/den
    #
    # Disequilibrium
    Da <- n11/n - (n1/(2*n))^2
    # Chi-square
    chisq <- ((n*Da)^2)/(n*p1^2) + ((-2*n*Da)^2)/(2*n*p1*(1-p1)) + ((n*Da)^2)/(n*(1-p1)^2)
    p <- 1 - pchisq(chisq, df=1)
    return(c(prob, Da, chisq, p))
  }
  #  tmp[gt == "0/0"] <- 0
  #  tmp[gt == "0/1"] <- 1
  #  tmp[gt == "1/0"] <- 1
  #  tmp[gt == "1/1"] <- 2
  gt <- extract.gt(x, element = "GT", mask = rep(TRUE, times=nrow(x@var.info)))
  tmp <- matrix(ncol=ncol(gt), nrow=nrow(gt))
  tmp[gt == "0/0" | gt == "0|0"] <- 0
  tmp[gt == "0/1" | gt == "0|1"] <- 1
  tmp[gt == "1/0" | gt == "1|0"] <- 1
  tmp[gt == "1/1" | gt == "1|1"] <- 2
  #
  gt <- tmp
  rm(tmp)
  #
  mask <- x@var.info$mask
  summ <- matrix(ncol=19, nrow=nrow(gt), 
                 dimnames=list(c(),
                               c('n', 'RR','RA','AA','nAllele','nREF','nALT','Ho','He',
                                 'hwe.prob', 'hwe.Da', 'hwe.chisq', 'hwe.p',
                                 'Ne','theta_pi','theta_w','theta_h','tajimas_d', 'faywu_h'))
  )
  #
  # Homozygous for reference allele.
  summ[mask,'RR'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==0)}))
  # Heterozygote.
  summ[mask,'RA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1)}))
  # Homozygous for alternate allele.
  summ[mask,'AA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==2)}))
  #
  summ[mask, 'n'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(!is.na(x))}))
  #
  summ[mask,'nREF'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                    function(x){sum(2*length(na.omit(x))-sum(x))})
  )
  summ[mask,'nALT'] <- rowSums(gt[mask, , drop=FALSE])
  summ[,'nAllele'] <- summ[,'nREF']+summ[,'nALT']
  #
  # Observed heterozygosity
  summ[mask,'Ho'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1)/length(na.omit(x))}))
  #
  # Expected heterozygosity
  summ[,'He'] <- 1 - ((summ[,'nREF']/summ[,'nAllele'])^2 + (summ[,'nALT']/summ[,'nAllele'])^2)
  #
  summ[,'Ne'] <- 1/(1-summ[,'He'])
  #
  # Hardy-Weinberg Disequilibrium
  summ[mask,c('hwe.prob', 'hwe.Da', 'hwe.chisq', 'hwe.p')] <- t(apply(summ[mask,c('RR', 'RA', 'AA'), drop=FALSE], MARGIN=1, FUN=hwe))
  #
  # Thetas.
  summ[,c('theta_pi','theta_w','theta_h')] <- t(apply(summ[,c('nREF','nALT'), drop=FALSE], MARGIN=1,thetas))
  #summ[,7:9] <- t(apply(summ[,c('nREF','nALT'), drop=FALSE], MARGIN=1,thetas))
  #
  summ[,'tajimas_d'] <- summ[,'theta_pi'] - summ[,'theta_w']
  summ[,'faywu_h'] <- summ[,'theta_pi'] - summ[,'theta_h']
  #  summ[,10] <- summ[,7] - summ[,8]
  #  summ[,11] <- summ[,7] - summ[,9]
  #
  #  print(head(summ))
  #  x@vcf.stat <- as.data.frame(summ)
  x@var.info <- cbind(x@var.info, as.data.frame(summ))
  return(x)
}

thetas <- function(x){
  #  print(x)
  rnum <- x[1]
  anum <- x[2]
  if(is.na(rnum)){return(c(NA,NA,NA))}
  n <- rnum + anum
  Si <- vector(mode="numeric", length=n)
  Si[anum] <- 1
  theta_w <- sum(1/1:(rnum+anum-1))^-1 * 1
  theta_pi <- (2*anum*rnum)/(n*(n-1))
  theta_h <- (2*1*anum^2)/(n*(n-1))
  return(c(theta_pi, theta_w, theta_h))
}

linkage <- function(x){
  gt <- x@gt.m
  mask <- x@mask
  link.m <- matrix(ncol=8, nrow=nrow(gt)-1,
                   dimnames=list(c(), c('pos', 'len', 'bigD', 'Delta', 'Dprime', 'delta', 'd', 'Q'))
  )
  link <- function(x){
    n1 <- length(!is.na(gt[x,]))
    n2 <- length(!is.na(gt[x+1,]))
    #    print(x)
  }
  lapply(1:nrow(link.m), link)
  #  print(head(gt))
  
  return(x)
}


#' @rdname Chrom-methods
#' @export
#' @aliases proc.chrom
#' 
#' @param verbose logical stating whether to produce verbose output
#'
#' 
#proc.chrom <- function(x, pop1=NA, pop2=NA, win.size=1000, max.win=10000, verbose=TRUE){
proc.chrom <- function(x, verbose=TRUE, ...){
  stopifnot(class(x) == "Chrom")
  #  x <- set.pop1(x, pop1)
  #  x <- set.pop2(x, pop2)
  #  ptime <- system.time(x@acgt.w <- regex.win(x))
  ptime <- system.time(x@seq.info$nuc.win <- regex.win(x))
  if(verbose==TRUE){
    message("Nucleotide regions complete.\n")
    print(ptime)
  }
  ptime <- system.time(x@seq.info$N.win <- regex.win(x, regex="[n]"))
  #  ptime <- system.time(x@n.w <- acgt.win(x, regex="[n]"))
  #  ptime <- system.time(x <- n.win(x))
  if(verbose==TRUE){
    message("N regions complete.\n")
    print(ptime)
  }
  ptime <- system.time(x@win.info <- var.win(x, ...))
  if(verbose==TRUE){
    message("Window analysis complete.\n")
    print(ptime)
  }
  #
  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
    ptime <- system.time(x <- gt2popsum(x))
    if(verbose==TRUE){
      message("Population summary complete.\n")
      print(ptime)
    }
  }
  #  ptime <- system.time(x <- windowize(x, win.size=win.size, max.win=max.win))
  #  ptime <- system.time(x <- windowize(x))
  #  if(verbose==TRUE){
  #    cat("Sliding windows created.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- gc.win(x))
  #  if(verbose==TRUE){
  #    cat("Sliding GC windows complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- snp.win(x))
  #  if(verbose==TRUE){
  #    cat("Sliding SNP windows complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- vcf.fix2gt.m(x))
  #  if(verbose==TRUE){
  #    cat("Genotype matrix complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- gt.m2sfs(x))
  #  cat("gt.m2sfs is commented out\n")
  #  if(verbose==TRUE){
  #    cat("SFS complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- linkage(x))
  #  if(verbose==TRUE){
  #    cat("Linkage calculation complete.\n")
  #    print(ptime)
  #  }
  return(x)
}

#### Graphic functions ####


plot.sfs <- function(x, log10=TRUE, ...){
  sfs <- x@sfs
  if(log10){sfs <- log10(sfs)}
  #
  layout(matrix(c(1,2), nrow=1), widths=c(4,1))
  image(t(sfs)[,nrow(sfs):1], col=rainbow(100, end=0.85),
        axes=FALSE, frame.plot=TRUE)
  #  axis(side=1, at=seq(1,ncol(sfs), by=1)/ncol(sfs), labels=NA)
  axis(side=1, at=seq(0, ncol(sfs)-1, by=1)/(ncol(sfs)-1), labels=NA)
  axis(side=1, at=seq(0, ncol(sfs)-1, by=5)/(ncol(sfs)-1), labels=seq(0, ncol(sfs)-1, by=5), las=1, tcl=-0.7)
  axis(side=3, at=seq(0, ncol(sfs)-1, by=1)/(ncol(sfs)-1), labels=NA)
  axis(side=2, at=seq(0, nrow(sfs)-1, by=1)/(nrow(sfs)-1), labels=NA)
  axis(side=2, at=seq(0, nrow(sfs)-1, by=5)/(nrow(sfs)-1), labels=seq(0, nrow(sfs)-1, by=5), las=1, tcl=-0.7)
  axis(side=4, at=seq(0, nrow(sfs)-1, by=1)/(nrow(sfs)-1), labels=NA)
  abline(a=0, b=1)
  title(main=paste("SFS for", x@name))
  #
  par(mar=c(5,0,4,3))
  barplot(height=rep(1, times=100), width=1, space=0,
          col=rainbow(100, start=0, end=0.85), border=NA, horiz=TRUE, axes=FALSE)
  axis(side=4, at=seq(0,100, length.out=2),
       labels=format(seq(0, 10^max(sfs, na.rm=TRUE), length.out=2), digits=3),
       las=1)
  axis(side=4, at=seq(1, max(sfs, na.rm=TRUE), by=1)*(100/max(sfs, na.rm=TRUE)),
       labels=10^seq(1, max(sfs, na.rm=TRUE), by=1), las=1
  )
  #
  par(mar=c(5,4,4,2), mfrow=c(1,1))
}








#### EOF ####
