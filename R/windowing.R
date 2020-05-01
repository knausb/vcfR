#' @title Create window summaries of data
#' @name Windowing
#' @rdname windowing
#' 
#' 
#' @description
#' Create windows of non-overlapping data and summarize.
#' 
#' @param x A NumericMatrix
#' @param pos A vector of chromosomal positions for each row of data (variants) 
#' @param maxbp Length of chromosome
#' @param winsize Size (in bp) for windows
#' @param depr logical (T/F), this function has been deprecated, set to FALSE to override.
#' 
#' @details
#' The numeric matrix where samples are in columns and variant data are in rows.
#' The windowing process therefore occurs along columns of data.
#' This matrix could be created with \code{\link{extract.gt}}.
#' 
#' The chromosome is expected to contain positions 1 though maxbp.
#' If maxbp is not specified this can be inferred from the last element in pos.
#' 
#'
#' @param starts integer vector of starting positions for windows
#' @param ends integer vector of ending positions for windows
#' @param summary string indicating type of summary (mean, median, sum)
#' 


# ' @rdname windowing
# @aliases windowing alias NM2winNM
#' 
#' @export
#' 
NM2winNM <- function(x, pos, maxbp, winsize = 100L, depr = TRUE) {
  
  if( depr ){
    myMsg <- "The function NM2winNM was deprecated in vcfR version 1.6.0. If you use this function and would like to advocate for its persistence, please contact the maintainer of vcfR. The maintainer can be contacted at maintainer('vcfR')"
    stop(myMsg)
  }
  
  .NM2winNM(x, pos, maxbp, winsize)
}



#' @rdname windowing
#' @export
#' 
z.score <- function(x){
  winave <- apply(x, MARGIN=2, mean, na.rm=TRUE)
  winsd  <- apply(x, MARGIN=2, stats::sd, na.rm=TRUE)
  zsc <- sweep(x, MARGIN=2, STATS=winave, FUN="-")
  zsc <- sweep(zsc, MARGIN=2, STATS=winsd, FUN="/")
  zsc
}


#' @rdname windowing
#' 
#' 
#' @export
#' 
windowize.NM <- function(x, pos, starts, ends, summary="mean", depr = TRUE){
  
  if( depr ){
    myMsg <- "The function windowizeNM was deprecated in vcfR version 1.6.0. If you use this function and would like to advocate for its persistence, please contact the maintainer of vcfR. The maintainer can be contacted at maintainer('vcfR')"
    stop(myMsg)
  }
  
  .windowize_NM(x, pos, starts, ends, summary=summary)  
}

