#' @title Query genotypes for heterozygotes
#' @name is.het
#' @rdname is_het
#' 
#' @description Query a matrix of genotypes for heterozygotes
#' 
#' 
#' @aliases is.het
#' 
#' @param x a matrix of genotypes
#' @param na_is_false should missing data be returned as NA (FALSE) or FALSE (TRUE)
#' 
#' @details
#' 
#' This function was designed to identify heterozygous positions in a matrix of genotypes.
#' The matrix of genotypes can be created with \code{\link{extract.gt}}.
#' Because the goal was to identify heterozygotes it may be reasonable to ignore missing values by setting na_is_false to TRUE so that the resulting matrix will consist of only TRUE and FALSE.
#' In order to preserve missing data as missing na_is_false can be set to FALSE where if at least one allele is missing NA is returned. 
#' 
#' 
#' @seealso
#' \code{\link{extract.gt}}
#' 
#' @examples 
#' data(vcfR_test)
#' gt <- extract.gt(vcfR_test)
#' hets <- is_het(gt)
#' # Censor non-heterozygous positions.
#' is.na(vcfR_test@gt[,-1][!hets]) <- TRUE
#' 
#' @export
is.het <- function(x, na_is_false = TRUE){
  if( class(x) != 'matrix' ){
    stop( paste( "Expecting a matrix, received a",  class(x) ) )
  }
  
  test_gt <- function(x, na_is_false = na_is_false){
    is.na( x[ x=="." ] ) <- TRUE
    
    if( sum( is.na(x) ) > 0 &  na_is_false == FALSE ){
      return(NA)
    } else {
      x <- unique(x)
      if( length(x) > 1 ){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
  
  proc_gt <- function(x, na_is_false = na_is_false){
    x <- strsplit(x, split="[/\\|]")
    x <- lapply(x, test_gt, na_is_false)
    unlist(x)
  }
  
  x2 <- apply( x, MARGIN=2, proc_gt, na_is_false = na_is_false )

  return(x2)
}
