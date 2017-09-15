#' @title Ranking variants within windows
#' @name Ranking
#' @rdname ranking
#' 
#' @description
#' Rank variants within windows.
#' 
#' @param x an object of class Crhom or a data.frame containing...
# @param ends a vector containing the position of the end of each window
#' @param scores a vector of scores for each variant to be used to rank the data
#' 
#' 
#' 


#' @rdname ranking
#' @aliases rank.variants.chromR
#' 
#' @export
rank.variants.chromR <- function(x, scores){
  if( class(x) != "chromR" ){
    stop("expecting object of class chromR or data.frame")
  }
  stopifnot(class(x@var.info) == 'data.frame')
  stopifnot(is.vector(x@win.info$end))
  stopifnot(class(x@win.info$end) == 'numeric')
#  stopifnot(is.vector(x@win.info['end']))
#  stopifnot(class(x@win.info['end']) == 'numeric')
  
  stopifnot(is.vector(scores))
  stopifnot(class(scores) == 'numeric')
  
  if( nrow(x@vcf) != length(scores) ){
    msg <- "The number of variants and scores do not match."
    msg <- paste(msg, " nrow(x@vcf): ", nrow(x@vcf), sep = "")
    msg <- paste(msg, ", length(scores): ", length(scores), sep = "")
    stop(msg)
  }
  
  x@var.info <- .rank_variants(x@var.info, x@win.info$end, scores)
  
  return(x)
}


