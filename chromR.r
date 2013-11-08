# chromR.
##### ##### ##### ##### #####
# Class definition.

setOldClass("DNAbin")

setClass(
  Class="Chrom",
  representation=representation(
    name = "character",
    seq = "DNAbin",
    len = "integer",
    vcf.meta = "character",
    vcf.fix = "data.frame",
    vcf.gt = "data.frame",
    vcf.info = "data.frame",
    ann = "data.frame",
    #
    acgt.w = "matrix",
    n.w = "matrix",
    windows = "matrix",
    nuccomp.w = "data.frame",
    snpden.w = "data.frame",
    #
    mask = "logical"
  ),
  prototype=prototype(
    vcf.fix = data.frame(matrix(ncol=8, nrow=0, 
                              dimnames=list(c(),
c('chrom','pos','id','ref','alt','qual','filter','info'))),
                       stringsAsFactors=FALSE)
#, dimnames=list(c(), c('chrom','pos','id','ref','alt','qual','filter','info')))
#    vcf.fix = matrix(ncol=8, nrow=0, dimnames=list(c(), c('chrom','pos','id','ref','alt','qual','filter','info')))
  )
)

##### ##### ##### ##### #####
# Generic methods.

setMethod(
  f="show",
  signature = "Chrom",
  definition=function(object){
    cat("*** Class Chrom, method Show *** \n")
    cat(paste("Name: ", object@name, "\n"))
    cat(paste("Length: ", object@len, "\n"))
    cat("Use print(object) for more details.\n")
    cat("******* End Show (Chrom) ******* \n")
  }
)

setMethod(
  f="print",
  signature="Chrom",
  definition=function (x,y,...){
    cat("***** Object of class 'Chrom' *****\n")
    cat(paste("Name: ", x@name, "\n"))
    cat(paste("Length: ", x@len, "\n"))
    cat("\nVCF fixed data:\n")
    print(head(x@vcf.fix[,1:7], n=4))
    cat("Last column (info) omitted.\n")
    cat("\nVCF variable data:\n")
    cat(paste("Columns: ", ncol(x@vcf.gt), "\n"))
    cat(paste("Rows: ", nrow(x@vcf.gt), "\n"))
    cat("(First column is format.)\n")
    cat("\nAnnotation data:\n")
    if(length(x@ann)>0){
      print(head(x@ann[,1:8], n=4))
      cat("Last column (attributes) omitted.\n")
    } else {
      cat("Empty slot.\n")
    }
    cat("***** End print (Chrom) ***** \n")
  }
)

setMethod(
  f= "plot",
  signature= "Chrom",
  definition=function (x,y,...){
#    cat("***** Object of class 'Chrom' *****\n")
#    cat("***** Plot not yet implemented *****\n")
    par(mfrow=c(2,2))
    if(length(x@vcf.info)>0){
      hist(x@vcf.info$dp, col=3, main="Depth (DP)", xlab="")
      rug(x@vcf.info$dp)
    } else {
      plot(1:2,1:2, type='n')
    }
    if(length(x@vcf.info)>0){
      hist(x@vcf.info$mq, col=4, main="Mapping quality (MQ)", xlab="")
      rug(x@vcf.info$mq)
    } else {
      plot(1:2,1:2, type='n')
    }
    if(length(x@vcf.fix)>0){
      hist(x@vcf.fix$qual, col=5, main="Quality (QUAL)", xlab="")
      rug(x@vcf.fix$qual)
    } else {
      plot(1:2,1:2, type='n')
    }
    if(length(x@snpden.w)>0){
      hist(x@snpden.w$count, col=6, main="SNP count (per window)", xlab="")
      rug(x@snpden.w$count)
    } else {
      plot(1:2,1:2, type='n')
    }
    par(mfrow=c(1,1))
  }
)

##### ##### ##### ##### #####
# Accessors.

### Getter for "names"
setGeneric("getName",function(object){standardGeneric ("getName")})

setMethod("getName","Chrom",
  function(object){
    return(object@name)
  }
)

# Setter for name.

setGeneric("setName<-",function(object,value){standardGeneric("setName<-")})

setReplaceMethod(
  f="setName",
  signature="Chrom",
  definition=function(object,value){
    object@name <-value
    return (object)
  }
)

# Setter for seq.

setGeneric("seq2chrom<-",function(object,value){standardGeneric("seq2chrom<-")})

setReplaceMethod(
  f="seq2chrom",
  signature="Chrom",
  definition=function(object,value){
    object@seq <-value
    object@len <-length(value[[1]])
    return (object)
  }
)

##### ##### ##### ##### #####
# Data loading functions.

vcf2chrom <- function(x,y,...){
  x@vcf.fix <- as.data.frame(y@fix)
#  x@vcf.fix <- as.data.frame(y[,1:8])
  colnames(x@vcf.fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
  x@vcf.fix[,2] <- as.numeric(x@vcf.fix[,2])
  x@vcf.fix[,6] <- as.numeric(x@vcf.fix[,6])
  #
  x@vcf.gt <- y@gt
#  x@vcf.gt <- y[,9:ncol(y)]
  #
  x@vcf.meta <- y@meta
  #
  info <- matrix(ncol=2, nrow=nrow(y@fix))
  colnames(info) <- c('dp','mq')
  info[,1] <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^DP=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  info[,2] <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^MQ=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  x@vcf.info <- as.data.frame(info)
  #
  x@mask <- rep(TRUE, times=nrow(x@vcf.fix))
  # assign may be more efficient.
  return(x)
}

ann2chrom <- function(x,y,...){
  x@ann <- as.data.frame(y)
  colnames(x@ann) <- c('seqid','source','type','start','end','score','strand','phase','attributes')
  x@ann$start <- as.numeric(as.character(x@ann$start))
  x@ann$end   <- as.numeric(as.character(x@ann$end))
  return(x)
}

create.chrom <- function(name, seq, vcf=NULL, ann=NULL){
  x <- new(Class="Chrom")
  setName(x) <- name
  if(class(seq)=="DNAbin"){
    seq2chrom(x) <- seq
  } else {
    cat("** Error: seq is not of class DNAbin** \n")
    break
  }
  if(length(vcf)>0){
    x <- vcf2chrom(x, vcf)
  }
  if(length(ann)>0){
    x <- ann2chrom(x, ann)
  }
  return(x)
}

##### ##### Set a mask ##### #####

masker <- function(x, QUAL=999, mindp=0.25, maxdp=0.75, minmq=0.25, maxmq=0.75, ...){
  x@mask <- as.numeric(x@vcf.fix[,6]) >= QUAL
  #
  x@mask[x@vcf.info$dp[x@vcf.info$dp <= quantile(x@vcf.info$dp, probs=c(mindp))]] <- FALSE
  x@mask[x@vcf.info$dp[x@vcf.info$dp >= quantile(x@vcf.info$dp, probs=c(maxdp))]] <- FALSE
  #
  x@mask[x@vcf.info$mq[x@vcf.info$mq <= quantile(x@vcf.info$mq, probs=c(mindp))]] <- FALSE
  x@mask[x@vcf.info$mq[x@vcf.info$mq >= quantile(x@vcf.info$mq, probs=c(maxdp))]] <- FALSE
  #
  return(x)
}

##### ##### Window functions ##### #####

acgt.win <- function(x, max.win=1000, regex="[acgtwsmkrybdhv]"){
  seq <- as.character(x@seq)[[1]]
  seq <- grep(regex, seq, ignore.case=T, perl=TRUE)
  bp.windows <- matrix(NA, ncol=2, nrow=max.win)
  bp.windows[1,1] <- seq[1]
  i <- 1
  for(j in 2:length(seq)){
    if(seq[j]-seq[j-1] > 1){
      bp.windows[i,2] <- seq[j-1]
      i <- i+1
      bp.windows[i,1] <- seq[j]
    }
  }
  bp.windows[i,2] <- seq[j]
  bp.windows <- bp.windows[1:i,]
  x@acgt.w <- bp.windows
  return(x)
}

n.win <- function(x, max.win=1000, regex="[n]"){
  seq <- as.character(x@seq)[[1]]
  seq <- grep(regex, seq, ignore.case=T, perl=TRUE)
  bp.windows <- matrix(NA, ncol=2, nrow=max.win)
  bp.windows[1,1] <- seq[1]
  i <- 1
  for(j in 2:length(seq)){
    if(seq[j]-seq[j-1] > 1){
      bp.windows[i,2] <- seq[j-1]
      i <- i+1
      bp.windows[i,1] <- seq[j]
    }
  }
  bp.windows[i,2] <- seq[j]
  bp.windows <- bp.windows[1:i,]
  x@n.w <- bp.windows
  return(x)
}

windowize <- function(x, win.size=1000, max.win=5000){
  acgt.w <- x@acgt.w
  windows <- matrix(NA, ncol=2, nrow=max.win)
  i <- 1
  for(j in 1:nrow(acgt.w)){
    while(acgt.w[j,2]-acgt.w[j,1] > win.size){
      windows[i,1] <- acgt.w[j,1]
      windows[i,2] <- acgt.w[j,1] + win.size - 1
      acgt.w[j,1] <- acgt.w[j,1] + win.size + 0
      i <- i+1
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
  vcf <- x@vcf.fix[x@mask,]$pos
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

proc.chrom <- function(x, win.size=1000, max.win=5000){
  x <- acgt.win(x)
  x <- n.win(x)
  x <- windowize(x, win.size=win.size, max.win=max.win)
  x <- gc.win(x)
  x <- snp.win(x)
  return(x)
}

##### ##### ##### ##### #####
# Graphic function.

chromo <- function(x, verbose=TRUE, nsum=TRUE, DP=TRUE, QUAL=TRUE, MQ=TRUE, SNPDEN=TRUE, NUC=TRUE, ANN=TRUE, ...){
  brows <- 0
  srows <- 0
  #
  if(length(x@vcf.info)>0 & DP){brows <- brows+1} # dp
  if(length(x@vcf.info)>0 & MQ){brows <- brows+1} # mq
  if(length(x@vcf.fix)>0 & QUAL){brows <- brows+1} # qual
  if(length(x@snpden.w)>0 & SNPDEN){brows <- brows+1}
  if(length(x@nuccomp.w)>0 & NUC){brows <- brows+1}
  #
  if(length(x@ann)>0 & ANN){srows <- srows+1}
  if(length(x@acgt.w)>0){srows <- srows+1}
  #
  if(verbose){
    cat('  Chromo\n')
    cat(paste("  brows: ", brows, "\n"))
    cat(paste("  srows: ", srows, "\n"))
  }
  #
  layout(matrix(c(1:((brows+srows)*2)), ncol=2, byrow=T),
    widths=rep(c(1,0.1), times=c(brows+srows)),
    heights=c(rep(1,times=brows),rep(0.4, times=srows)))
  par(mar=c(0,0,0,0))
  par(oma=c(4,4,3,1))
  #
  if(length(x@vcf.info)>0 & DP){ # dp
    plot(x@vcf.fix[x@mask,2], x@vcf.info[x@mask,1], pch=20, col="#0080ff22", axes=F, frame.plot=T, ylab="", ...)
    axis(side=2, las=2)
    title(main="Read depth (DP)", line=-1)
    boxplot(x@vcf.info[x@mask,1], axes=FALSE, frame.plot=T, col="#0080ff")
  }
  #
  if(length(x@vcf.info)>0 & MQ){ # dp
    plot(x@vcf.fix[x@mask,2], x@vcf.info[x@mask,2], pch=20, col="#3CB37122", axes=F, frame.plot=T, ylab="", ...)
    axis(side=2, las=2)
    title(main="Mapping quality (MQ)", line=-1)
    boxplot(x@vcf.info[x@mask,2], axes=FALSE, frame.plot=T, col="#3CB371")
  }
  #
  if(length(x@vcf.fix)>0 & QUAL){ # qual
    plot(x@vcf.fix[x@mask,2], x@vcf.fix[x@mask,6], pch=20, col="#80008022", axes=F, frame.plot=T, ylab="", ...)
    axis(side=2, las=2)
    title(main="QUAL", line=-1)
    boxplot(as.numeric(x@vcf.fix[x@mask,6]), axes=FALSE, frame.plot=T, col="#800080")
  }
  #
  if(length(x@snpden.w)>0 & SNPDEN){
    # SNP density.
    plot(c(0,x@len), c(0,max(x@snpden.w$density)), type='n', xlab="", ylab="", axes=F, frame.plot=T, ...)
    abline(h=seq(0.1, 1, by=0.1), col="#a0a0a0")
    abline(h=seq(0.02, 0.08, by=0.02), lty=3, col="#a0a0a0")
    rect(x@snpden.w$start, 0, x@snpden.w$stop, x@snpden.w$density, col="#cc0000", border=NA)
    axis(side=2, las=2)
    title(main="SNPs per site", line=-1)
    title(main=paste(sum(x@snpden.w$count), "total SNPs"), line=-2)
    boxplot(x@snpden.w$density, axes=FALSE, frame.plot=T, col="#cc0000", ylim=c(0,max(x@snpden.w$density)))
  }
  if(length(x@nuccomp.w)>0 & NUC){
    # GC and AT content.
    plot(c(0,x@len), c(0,1), type='n', xlab="", ylab="", axes=F, frame.plot=T, ...)
    rect(x@nuccomp.w[,2], 0, x@nuccomp.w[,3], x@nuccomp.w[,6], col="#0000cc", border=NA)
    rect(x@nuccomp.w[,2], x@nuccomp.w[,6], x@nuccomp.w[,3], x@nuccomp.w[,6]+x@nuccomp.w[,7], col="#ffd700", border=NA)
    axis(side=2, las=2)
    title(main="GC and AT content", line=-1)
    boxplot(x@nuccomp.w[,6:7], axes=FALSE, frame.plot=T, col=c("#0000cc", "#ffd700"), ylim=c(0,1), border=c("#0000cc", "#ffd700"))
    text(1,0.1, "G/C")
    text(2,0.1, "A/T")
  }
  #
  if(length(x@ann)>0 & ANN){
    # Annotations.
    plot(c(0,x@len), c(-1,1), type='n', xlab="", ylab="", las=1, axes=FALSE, frame.plot=TRUE, ...)
    lines(c(0,x@len),c(0, 0), lwd=2)
    rect(as.numeric(as.character(x@ann[,4])), -1, as.numeric(as.character(x@ann[,5])), 1, col="#b22222", border=NA)
    #
    plot(1:10,1:10,type='n', axes=FALSE, xlab="", ylab="", ...)
    if(nsum){
      text(5,10,"Genic bases:")
      text(5,9,sum(as.numeric(as.character(x@ann[,5]))-as.numeric(as.character(x@ann[,4]))))
      text(5,8,format(sum(as.numeric(as.character(x@ann[,5]))-as.numeric(as.character(x@ann[,4])))/x@len, digits=3))
    }
  }
  #
  if(length(x@acgt.w)>0){
    # Chromosome.
    plot(c(0,x@len), c(-1,1), type='n', xlab="", ylab="", las=1, axes=FALSE, frame.plot=TRUE, ...)
    lines(c(0,x@len),c(0, 0), lwd=2)
    rect(x@acgt.w[,1], -0.7, x@acgt.w[,2], 0.7, col="#00cc00", border=NA)
    rect(x@n.w[,1], -0.4, x@n.w[,2], 0.4, col="#ff6666", border=NA)
    axis(side=1)
    title(xlab="Base pairs", line=2, outer=T)
    #
    plot(1:10,1:10,type='n', axes=FALSE, xlab="", ylab="")
    if(nsum){
      text(5,9,"Length (bp)")
      text(5,8, x@len)
      text(5,6,"GC (0.25; 0.5; 0.75)")
      text(5,5, paste(format(quantile(x@nuccomp.w$gcf, probs=c(0.25,0.5,0.75)), digits=3), collapse="; "))
      text(5,3, "SNPs (0.25; 0.5; 0.75)")
      text(5,2, paste(format(quantile(x@snpden.w$density, probs=c(0.25,0.5,0.75)), digits=3), collapse="; "))
    }
  }
  title(main=x@name, outer=TRUE, line=1)
  #
  # Return to defaults
  layout(matrix(c(1), ncol=1, byrow=T))
  par(mar=c(5,4,4,2))
  par(oma=c(0,0,0,0))
}


##### ##### ##### ##### #####
# EOF.
