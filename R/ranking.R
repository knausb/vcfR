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
  
  x@var.info <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', x@var.info, x@win.info$end, scores)
  #  vars <- .Call('vcfR_rank_variants', PACKAGE = 'vcfR', pinf_mt@var.info, pinf_mt@win.info$end, testv)
  
  return(x)
}


