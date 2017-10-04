
#' @title Plot freq_peak object
#' @name freq_peak_plot
#' @rdname freq_peak_plot
#' 
#' @description
#' Converts allele balance data produced by \code{freq_peak()} to a copy number by assinging the allele balance data (frequencies) to its closest expected ratio.
#'  
#' @param pos chromosomal position of variants
#' @param ab1 matrix of allele balances for allele 1
#' @param ab2 matrix of allele balances for allele 2
#' @param fp1 feq_peak object for allele 1
#' @param fp2 feq_peak object for allele 2
#' @param mySamp sample indicator
#' @param col1 color 1
#' @param col2 color 2
#' @param mhist logical indicating to include a marginal histogram
#' @param ... parameters passed on to other functions
#' 
#' @details 
#' Converts allele balance data produced by \code{freq_peak()} to a copy number.
#' 
#' @return an invisible NULL.
#' 
#' 
#' @examples
#' data(vcfR_example)
#' freq_peak_plot(pos=1:40)
#' 
#' @export
freq_peak_plot <- function(pos, 
                           ab1 = NULL, 
                           ab2 = NULL, 
                           fp1 = NULL, 
                           fp2 = NULL, 
                           mySamp = 1, 
                           col1 = "#A6CEE344",
                           col2 = "#1F78B444",
                           mhist = TRUE,
                           ...){
  
  if( !inherits(fp1, "freq_peak") & !is.null(fp1) ){
    msg <- "fp1 does not appear to be a freq_peak object"
    stop(msg)
  }
  if( !inherits(fp2, "freq_peak") & !is.null(fp2) ){
    msg <- "fp2 does not appear to be a freq_peak object"
    stop(msg)
  }
  
  # Store original options.
  orig_opts <- options()
  
  # Determine plot geometry.
  if( mhist == TRUE ){
    graphics::layout(matrix(1:2, nrow=1), widths = c(4,1))
    graphics::par(mar=c(5,4,4,0))
  }
  
  # Initialize plot
  plot( range(pos, na.rm = TRUE), c(0,1), ylim=c(0,1), type="n", yaxt='n', 
       main = "", xlab = "Position", ylab = "Allele balance")
  graphics::axis(side=2, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
                 labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
  graphics::abline(h=c(0.2,0.25,0.333,0.5,0.666,0.75,0.8), col=8)  

  
  
  # Null marginal histogram
  if( mhist == TRUE & is.null(ab1) & is.null(ab2) ){
    graphics::par(mar=c(5,1,4,2))
    graphics::barplot(height=0.01, width=0.02,  space = 0, horiz = T, add = FALSE, col="#000000", xlim = c(0,1.0))
#    barplot(height=0.01, width=0.02,  space = 0, horiz = T, add = TRUE, col="#1F78B4")
  }
  
  
  options(orig_opts)
  return( invisible(NULL) )
}


