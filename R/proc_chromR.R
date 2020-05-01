
#' @title Process chromR object
#' @name Process chromR objects
#' @rdname proc_chromR
#' @description Functions which process chromR objects 
#' 
#' @param x object of class chromR
#' @param win.size integer indicating size for windowing processes
#' @param verbose logical indicating whether verbose output should be reported
# ' @param ... arguments to be passed to methods
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
#' Each \strong{window} receives a \strong{name} and its coordinates.
#' Several attempts are made to name the windows appropriately.
#' First, the CHROM column of vcfR@fix is queried for a name.
#' Next, the label of the sequence is queried for a name.
#' Next, the first cell of the annotation matrix is queried.
#' If an appropriate name was not found in the above locations the chromR object's 'name' slot is used.
#' Note that the 'name' slot has a default value.
#' If this default value is not updated then all of your windows may receive the same name.
#' 
#' 


# ' @rdname proc_chromR
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
  
  if(class(x@seq) == "DNAbin"){
    ptime <- system.time(x@seq.info$nuc.win <- seq2rects(x)) 
    if(verbose==TRUE){
      message("Nucleotide regions complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, chromosome representation not made (seq2rects)." )
  }
  
  if(class(x@seq) == "DNAbin"){
    ptime <- system.time(x@seq.info$N.win <- seq2rects(x, chars="n")) 
    if(verbose==TRUE){
      message("N regions complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))      
    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, chromosome representation not made (seq2rects, chars=n)." )
  }


  # Population summary
  if( nrow(x@vcf@gt) > 0 ){
    if( nrow( x@vcf@gt[ x@var.info$mask, , drop = FALSE ] ) > 0 ){
      ptime <- system.time(x <- gt.to.popsum(x))
      if(verbose==TRUE){
        message("Population summary complete.")
        message(paste("  elapsed time: ", round(ptime[3], digits=4)))
      }
    }
  }
  
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  
  # Initialize windows.
  if( length(x@len) > 0 ){
    ptime <- system.time(x@win.info <- .window_init(window_size=win.size, max_bp=x@len))

    # Name of windows based on chromosome name.
    if( !is.na(x@var.info$CHROM[1]) ){
      x@win.info <- cbind(rep(x@var.info$CHROM[1], times=nrow(x@win.info)), x@win.info)
      names(x@win.info)[1] <- "CHROM"
    } else if( !is.null(x@seq) ){
      x@win.info <- cbind(rep( labels(x@seq)[1], times=nrow(x@win.info)), x@win.info)
      names(x@win.info)[1] <- "CHROM"
    } else if( nrow(x@ann) > 0 ){
      x@win.info <- cbind(rep( x@ann[1,1], times=nrow(x@win.info)), x@win.info)
      names(x@win.info)[1] <- "CHROM"
    } else {
      x@win.info <- cbind(rep( x@name, times=nrow(x@win.info)), x@win.info)
      names(x@win.info)[1] <- "CHROM"
    }
    if(verbose==TRUE){
#    print("window_init complete.")
#    print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("window_init complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  }
#  }

  if(class(x@seq) == "DNAbin"){
#    if( nrow( x@vcf@gt[x@var.info$mask, , drop = FALSE ] ) > 0 ){
      ptime <- system.time(x@win.info <- .windowize_fasta(wins=x@win.info,
                                                          seq=as.character(x@seq)[1,]
                                                          ))
      if(verbose==TRUE){
#        print("windowize_fasta complete.")
#        print(paste("  elapsed time: ", round(ptime[3], digits=4)))
        message("windowize_fasta complete.")
        message(paste("  elapsed time: ", round(ptime[3], digits=4)))
      }
#    }
  } else if ( is.null( x@seq ) & verbose == TRUE ){
    warning( "seq slot is NULL, windowize_fasta not run." )
  }
  
  # Windowize annotations.
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  if( nrow(x@ann) > 0 ){
    #if( nrow( x@vcf@gt[x@var.info$mask, , drop = FALSE] ) > 0 ){
      ptime <- system.time(x@win.info <- .windowize_annotations(wins=x@win.info,
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
    #}
  } else if ( nrow(x@ann) == 0 ){
    if ( verbose == TRUE ){
      warning( "ann slot has zero rows." )
    }
    if( nrow(x@win.info) > 0 ){
      x@win.info$genic <- 0
    }
  }

  # Windowize variants.
#  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
  if( nrow( x@vcf@fix[x@var.info$mask, , drop = FALSE ] ) > 0 ){
    ptime <- system.time(x@win.info <- .windowize_variants(windows=x@win.info, variants=x@var.info[c('POS','mask')]))
    if(verbose==TRUE){
#      print("windowize_variants complete.")
#      print(paste("  elapsed time: ", round(ptime[3], digits=4)))
      message("windowize_variants complete.")
      message(paste("  elapsed time: ", round(ptime[3], digits=4)))
    }
  } else {
    if( nrow(x@win.info) > 0 ){
      x@win.info$variants <- 0
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

  rects <- .seq_to_rects(seq, targets=chars)
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





