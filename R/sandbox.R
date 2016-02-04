

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

##### ##### Set populations #####

# @rdname Chrom-methods
# @export
# @aliases set.pop1 
# 
# @param pop1 a numeric vector indicating the samples in population 1
# 
set.pop1 <- function(x, pop1){
  x@pop1 <- pop1
  return(x)  
}

# @rdname Chrom-methods
# @export
# @aliases set.pop2
# 
# @param pop2 a numeric vector indicating the samples in population 2
# 
set.pop2 <- function(x, pop2){
  x@pop2 <- pop2
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

#### Graphic functions ####


plot.sfs <- function(x, log10=TRUE, ...){
  sfs <- x@sfs
  if(log10){sfs <- log10(sfs)}
  #
  graphics::layout(matrix(c(1,2), nrow=1), widths=c(4,1))
  graphics::image(t(sfs)[,nrow(sfs):1], col=grDevices::rainbow(100, end=0.85),
        axes=FALSE, frame.plot=TRUE)
  #  axis(side=1, at=seq(1,ncol(sfs), by=1)/ncol(sfs), labels=NA)
  graphics::axis(side=1, at=seq(0, ncol(sfs)-1, by=1)/(ncol(sfs)-1), labels=NA)
  graphics::axis(side=1, at=seq(0, ncol(sfs)-1, by=5)/(ncol(sfs)-1), labels=seq(0, ncol(sfs)-1, by=5), las=1, tcl=-0.7)
  graphics::axis(side=3, at=seq(0, ncol(sfs)-1, by=1)/(ncol(sfs)-1), labels=NA)
  graphics::axis(side=2, at=seq(0, nrow(sfs)-1, by=1)/(nrow(sfs)-1), labels=NA)
  graphics::axis(side=2, at=seq(0, nrow(sfs)-1, by=5)/(nrow(sfs)-1), labels=seq(0, nrow(sfs)-1, by=5), las=1, tcl=-0.7)
  graphics::axis(side=4, at=seq(0, nrow(sfs)-1, by=1)/(nrow(sfs)-1), labels=NA)
  graphics::abline(a=0, b=1)
  graphics::title(main=paste("SFS for", x@name))
  #
  graphics::par(mar=c(5,0,4,3))
  graphics::barplot(height=rep(1, times=100), width=1, space=0,
          col=grDevices::rainbow(100, start=0, end=0.85), border=NA, horiz=TRUE, axes=FALSE)
  graphics::axis(side=4, at=seq(0,100, length.out=2),
       labels=format(seq(0, 10^max(sfs, na.rm=TRUE), length.out=2), digits=3),
       las=1)
  graphics::axis(side=4, at=seq(1, max(sfs, na.rm=TRUE), by=1)*(100/max(sfs, na.rm=TRUE)),
       labels=10^seq(1, max(sfs, na.rm=TRUE), by=1), las=1
  )
  #
  graphics::par(mar=c(5,4,4,2), mfrow=c(1,1))
}



