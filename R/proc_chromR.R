
#' @title Process chromR object
#' @name Process chromR objects
#' @rdname proc_chromR
#' @description Functions which process chromR objects 
#' 
#' @param x object of class chromR
#' @param win.size integer indicating size for windowing processes
#' @param verbose logical indicating whether verbose output should be reported
#' @param ... arguments to be passed to methods
#' @param max.win maximum window size
#' @param regex a regular expression to indicate nucleotides to be searched for
#' 
#' @details
#' The function \strong{proc_chromR()} calls helper functions to process the data present in a chromR object into summaries statistics.
#' 
#' The function \strong{regex.win()} is used to generate coordinates to define rectangles to represent regions of the chromosome containing called nucleotides (acgtwsmkrybdhv).
#' It is then called a second time to generate coordinates to define rectangles to represent regions called as uncalled nucleotides (n, but not gaps).
#' 
#' The function \strong{gt2popsum} is called to create summaries of the variant data.
#' 
#' The function \strong{var.win} is called to create windowized summaries of the chromR object.
#' 
#' 


#' @rdname proc_chromR
#' @export
#' @aliases proc.chromR
#'
proc.chromR <- function(x, win.size = 1e3, verbose=TRUE){
  stopifnot(class(x) == "chromR")
  
  if( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL." )
  }
  if( nrow(x@ann) == 0 & verbose == TRUE ){
    warning( "annotation slot has no rows." )
  }
  
  
#  ptime <- system.time(x@seq.info$nuc.win <- regex.win(x))
  if(class(x@seq) == "DNAbin"){
    ptime <- system.time(x@seq.info$nuc.win <- seq2rects(x)) 
    if(verbose==TRUE){
#      print("Nucleotide regions complete.")
#      print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("Nucleotide regions complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, chromosome representation not made (seq2rects)." )
  }
  
  if(class(x@seq) == "DNAbin"){
#  ptime <- system.time(x@seq.info$N.win <- regex.win(x, regex="[n]"))
    ptime <- system.time(x@seq.info$N.win <- seq2rects(x, chars="n")) 
    if(verbose==TRUE){
#      print("N regions complete.")
#      print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("N regions complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))      
    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, chromosome representation not made (seq2rects, chars=n)." )
  }

    
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  if(nrow(x@vcf@gt[x@var.info$mask,])>0){
#    ptime <- system.time(x <- gt2popsum(x))
    ptime <- system.time(x <- gt.to.popsum(x))
    if(verbose==TRUE){
#      print("Population summary complete.")
#      print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("Population summary complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  }
  
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  ptime <- system.time(x@win.info <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=win.size, max_bp=x@len))
  x@win.info <- cbind(rep(x@var.info$CHROM[1], times=nrow(x@win.info)), x@win.info)
  names(x@win.info)[1] <- "CHROM"
  if(verbose==TRUE){
#    print("window_init complete.")
#    print(paste("  elapsed time: ", round(ptime[3], digits=4)))
    message("window_init complete.")
    message(paste("  elapsed time: ", round(ptime[3], digits=4)))
  }
#  }

  if(class(x@seq) == "DNAbin"){
    if(nrow(x@vcf@gt[x@var.info$mask,])>0){
      ptime <- system.time(x@win.info <- .Call('vcfR_windowize_fasta', 
                                               PACKAGE = 'vcfR',
                                               wins=x@win.info,
                                               seq=as.character(x@seq)[1,]
                                               ))
      if(verbose==TRUE){
#        print("windowize_fasta complete.")
#        print(paste("  elapsed time: ", round(ptime[3], digits=4)))
        message("windowize_fasta complete.")
        message(paste("  elapsed time: ", round(ptime[3], digits=4)))
      }
    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, windowize_fasta not run." )
  }


#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  if( nrow(x@ann) > 0 ){
    if( nrow(x@vcf@gt[x@var.info$mask,]) > 0 ){
      ptime <- system.time(x@win.info <- .Call('vcfR_windowize_annotations', PACKAGE = 'vcfR', wins=x@win.info,
                                               ann_starts=as.numeric(as.character(x@ann[,4])), 
                                               ann_ends=as.numeric(as.character(x@ann[,5])),
                                               chrom_length=x@len)
      )
      if(verbose==TRUE){
#        print("windowize_annotations complete.")
#        print(paste("  elapsed time: ", round(ptime[3], digits=4)))
        message("windowize_annotations complete.")
        message(paste("  elapsed time: ", round(ptime[3], digits=4)))
      }
    }
  } else if ( nrow(x@ann) == 0 & verbose == TRUE ){
    warning( "ann slot has zero rows, windowize_annotations not run." )
  }

  
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  if( nrow(x@vcf@gt[x@var.info$mask,]) > 0 ){
    ptime <- system.time(x@win.info <- .Call('vcfR_windowize_variants', PACKAGE = 'vcfR', windows=x@win.info, variants=x@var.info[c('POS','mask')]))
    if(verbose==TRUE){
#      print("windowize_variants complete.")
#      print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("windowize_variants complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  }
  
  return(x)
}



##### ##### seq.info functions #####

#' @rdname proc_chromR
#' @export
#' @aliases regex.win
#' 
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


#' @rdname proc_chromR
#' @aliases seq2rects
#' 
#' @description
#' Create representation of a sequence.
#' Begining and end points are determined for stretches of nucleotides.
#' Stretches are determined by querying each nucleotides in a sequence to determine if it is represented in the database of characters (chars).
#' 
#' 
#' @param chars a vector of characters to be used as a database for inclusion in rectangles
#' @param lower converts the sequence and database to lower case, making the search case insensitive
#' 
#'   
#' @export
#' 
seq2rects <- function(x, chars="acgtwsmkrybdhv", lower=TRUE){

  if(is.matrix(as.character(x@seq))){
#    seq <- as.character(x@seq)[1:length(x@seq)]
    seq <- as.character(x@seq)[1,]
  }

  if(lower == TRUE){
    seq <- tolower(seq)
    chars <- tolower(chars)
  }

  rects <- .Call('vcfR_seq_to_rects', PACKAGE = 'vcfR', seq, targets=chars)
  return(rects)
}


#' @rdname proc_chromR
#' @export
#' @aliases var.win
#' 
#var.win <- function(x, win.size=1e3){
var.win <- function(x, win.size=1e3){
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  
  # Convert DNAbin to string of chars.
  if(class(x@seq) == "DNAbin"){
    if(is.matrix(as.character(x@seq))){
      seq <- as.character(x@seq)[1:length(x@seq)]    
    } else if(is.list(as.character(x@seq))){
      seq <- as.character(x@seq)[[1]]
    }
  }

  # Create a vector of 0 and 1 marking genic sites.
  if(nrow(x@ann) > 0){
    genic_sites <- rep(0, times=x@len)
    genic_sites[unlist(apply(x@ann[, 4:5], MARGIN=1, function(x){seq(from=x[1], to=x[2], by=1)}))] <- 1
  }
  
  # Initialize data.frame of windows.
  win.info <- seq(1, x@len, by=win.size)
  win.info <- cbind(win.info, c(win.info[-1]-1, x@len))
  win.info <- cbind(1:nrow(win.info), win.info)
  win.info <- cbind(win.info, win.info[,3]-win.info[,2]+1)
  #  win.info <- cbind(win.info, matrix(ncol=7, nrow=nrow(win.info)))

  # Declare a function to count nucleotide classes.
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
  
  # Implement function to count nucleotide classes.
  if(class(x@seq) == "DNAbin"){
    win.info <- cbind(win.info, t(apply(win.info, MARGIN=1, win.proc, seq=seq)))
    win.info <- as.data.frame(win.info)
    names(win.info) <- c('window','start','end','length','A','C','G','T','N','other','variants', 'genic')
  }
  
  win.info
}






#' @rdname proc_chromR
#' @export
#' @aliases gt2popsum
#' 
gt2popsum <- function(x){
  if(class(x) != "chromR"){stop("Object is not of class chromR")}
  #  stopifnot(class(x) == "chromR")
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
    p <- 1 - stats::pchisq(chisq, df=1)
    return(c(prob, Da, chisq, p))
  }
  #  tmp[gt == "0/0"] <- 0
  #  tmp[gt == "0/1"] <- 1
  #  tmp[gt == "1/0"] <- 1
  #  tmp[gt == "1/1"] <- 2
#  gt <- extract.gt(x, element = "GT", mask = rep(TRUE, times=nrow(x@var.info)))
  gt <- extract.gt(x, element = "GT")
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
                                  function(x){sum(x==0, na.rm=TRUE)}))
  # Heterozygote.
  summ[mask,'RA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1, na.rm=TRUE)}))
  # Homozygous for alternate allele.
  summ[mask,'AA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==2, na.rm=TRUE)}))
  #
  summ[mask, 'n'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(!is.na(x))}))
  #
  summ[mask,'nREF'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                    function(x){sum(2*length(stats::na.omit(x))-sum(x), na.rm=TRUE)})
  )
  summ[mask,'nALT'] <- rowSums(gt[mask, , drop=FALSE], na.rm=TRUE)
  summ[,'nAllele'] <- summ[,'nREF']+summ[,'nALT']
  #
  # Observed heterozygosity
  summ[mask,'Ho'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1, na.rm=TRUE)/length(stats::na.omit(x))}))
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




#' @rdname proc_chromR
#' @aliases gt.to.popsum
#' 
#' @export
gt.to.popsum <- function(x){
  if(class(x) != "chromR"){stop("Object is not of class chromR")}
  
  # Extract genotypes from vcf.gt
  gt <- extract.gt(x, element="GT")
  
  # If summaries already exist, we'll remove them.
  x@var.info <- x@var.info[,grep("^n$|^Allele_counts$|^He$|^Ne$", colnames(x@var.info), invert = TRUE)]
  
  x@var.info <- .Call('vcfR_gt_to_popsum', PACKAGE = 'vcfR', var_info=x@var.info, gt=gt)

  return(x)
}



