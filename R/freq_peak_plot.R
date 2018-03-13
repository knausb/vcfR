
#' @title Plot freq_peak object
#' @name freq_peak_plot
#' @rdname freq_peak_plot
#' 
#' @description
#' Converts allele balance data produced by \code{freq_peak()} to a copy number by assinging the allele balance data (frequencies) to its closest expected ratio.
#'  
#' @param pos chromosomal position of variants
#' @param posUnits units ('bp', 'Kbp', 'Mbp', 'Gbp') for `pos` to be converted to in the main plot
#' @param ab1 matrix of allele balances for allele 1
#' @param ab2 matrix of allele balances for allele 2
#' @param fp1 feq_peak object for allele 1
#' @param fp2 feq_peak object for allele 2
#' @param mySamp sample indicator
#' @param col1 color 1
#' @param col2 color 2
#' @param alpha sets the transparency for dot plot (0-255)
#' @param main main plot title.
#' @param mhist logical indicating to include a marginal histogram
#' @param layout call layout
#' @param ... parameters passed on to other functions
#' 
#' @details 
#' 
#' Creates a visualization of allele balance data consisting of a dot plot with position as the x-axis and frequency on the y-axis and an optional marginal histogram.
#' The only required information is a vector of chromosomal positions, however this is probably not going to create an interesting plot.
#' 
#' 
#' @return An invisible NULL.
#' 
#' @seealso 
#' freq_peak,
#' peak_to_ploid
#' 
#' @examples
#' 
#' # An empty plot.
#' freq_peak_plot(pos=1:40)
#' 
#' data(vcfR_example)
#' gt <- extract.gt(vcf)
#' hets <- is_het(gt)
#' # Censor non-heterozygous positions.
#' is.na(vcf@gt[,-1][!hets]) <- TRUE
#' # Extract allele depths.
#' ad <- extract.gt(vcf, element = "AD")
#' ad1 <- masplit(ad, record = 1)
#' ad2 <- masplit(ad, record = 2)
#' freq1 <- ad1/(ad1+ad2)
#' freq2 <- ad2/(ad1+ad2)
#' myPeaks1 <- freq_peak(freq1, getPOS(vcf))
#' is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
#' myPeaks2 <- freq_peak(freq2, getPOS(vcf), lhs = FALSE)
#' is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
#' freq_peak_plot(pos = getPOS(vcf), ab1 = freq1, ab2 = freq2, fp1 = myPeaks1, fp2=myPeaks2)
#' 
#' 
#' 
#' @export
freq_peak_plot <- function(pos, 
                           posUnits = 'bp',
                           ab1 = NULL, 
                           ab2 = NULL, 
                           fp1 = NULL, 
                           fp2 = NULL, 
                           mySamp = 1, 
                           col1 = "#A6CEE3",
                           col2 = "#1F78B4",
                           alpha = 44,
                           main = NULL,
                           mhist = TRUE,
                           layout = TRUE,
                           ...){
  
  if( !inherits(fp1, "freq_peak") & !is.null(fp1) ){
    msg <- "fp1 does not appear to be a freq_peak object"
    stop(msg)
  }
  if( !inherits(fp2, "freq_peak") & !is.null(fp2) ){
    msg <- "fp2 does not appear to be a freq_peak object"
    stop(msg)
  }
  
  # Handle x-axis units.
  if( posUnits == 'bp'){
    # Don't need to do anything.
  } else if(posUnits == 'Kbp'){
    pos <- pos/1e3
    fp1$wins[,'START_pos'] <- fp1$wins[,'START_pos']/1e3
    fp1$wins[,'END_pos']   <- fp1$wins[,'END_pos']/1e3
    fp2$wins[,'START_pos'] <- fp2$wins[,'START_pos']/1e3
    fp2$wins[,'END_pos']   <- fp2$wins[,'END_pos']/1e3
  } else if(posUnits == 'Mbp'){
    pos <- pos/1e6
    fp1$wins[,'START_pos'] <- fp1$wins[,'START_pos']/1e6
    fp1$wins[,'END_pos']   <- fp1$wins[,'END_pos']/1e6
    fp2$wins[,'START_pos'] <- fp2$wins[,'START_pos']/1e6
    fp2$wins[,'END_pos']   <- fp2$wins[,'END_pos']/1e6
  } else if(posUnits == 'Gbp'){
    pos <- pos/1e9
    fp1$wins[,'START_pos'] <- fp1$wins[,'START_pos']/1e9
    fp1$wins[,'END_pos'] <-   fp1$wins[,'END_pos']/1e9
    fp2$wins[,'START_pos'] <- fp2$wins[,'START_pos']/1e9
    fp2$wins[,'END_pos'] <-   fp2$wins[,'END_pos']/1e9
  }
  xlabel <- paste('Position (', posUnits, ')', sep = '')  
  
  # Handle color
  col1 <- grDevices::col2rgb(col1, alpha = FALSE)
  col2 <- grDevices::col2rgb(col2, alpha = FALSE)
  col1 <- grDevices::rgb(col1[1,1], col1[2,1], col1[3,1], maxColorValue=255)
  col2 <- grDevices::rgb(col2[1,1], col2[2,1], col2[3,1], maxColorValue=255)
  col1d <- paste(col1, alpha, sep = "")
  col2d <- paste(col2, alpha, sep = "")
  
  # Store original par.
  orig_par <- graphics::par(no.readonly = TRUE)
  
  # Determine plot geometry.
  if( mhist == TRUE & layout == TRUE){
    graphics::layout(matrix(1:2, nrow=1), widths = c(4,1))
    
  }
  graphics::par(mar=c(5,4,4,0))
  
  # Initialize plot
  plot( range(pos, na.rm = TRUE), c(0,1), ylim=c(0,1), type="n", yaxt='n', 
       main = "", xlab = xlabel, ylab = "Allele balance")
  graphics::axis(side=2, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
                 labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
  graphics::abline(h=c(0.2,0.25,0.333,0.5,0.666,0.75,0.8), col=8)  

  # Add dot plots
  if( !is.null(ab1) ){
    graphics::points(pos, ab1[,mySamp], pch = 20, col= col1d)
  }
  if( !is.null(ab2) ){
    graphics::points(pos, ab2[,mySamp], pch = 20, col= col2d)
  }
  
  # Add window peak indicators
  if( !is.null(fp1) ){
    graphics::segments(x0=fp1$wins[,'START_pos'], y0=fp1$peaks[,mySamp],
                       x1=fp1$wins[,'END_pos'], lwd=3)
  }
  if( !is.null(fp2) ){
    graphics::segments(x0=fp2$wins[,'START_pos'], y0=fp2$peaks[,mySamp],
                       x1=fp2$wins[,'END_pos'], lwd=3)
  }
  
  # Title
  if( !is.null(main) ){
    graphics::title(main = main)
  }
#  if( !is.null(ab1) ){
#    graphics::title(main = colnames(ab1[, mySamp, drop = F]))
#  } else if( !is.null(ab2) ){
#    graphics::title(main = colnames(ab2[, mySamp, drop = F]))
#  }
    
  # Marginal histogram
  if( mhist == TRUE){
    graphics::par(mar=c(5,1,4,2))
    
    if( !is.null(fp1$bin_width) ){
      hsbrk <- seq(0,1,by=fp1$bin_width)
    } else if( !is.null(fp2$bin_width) ){
      hsbrk <- seq(0,1,by=fp2$bin_width)
    } else {
      hsbrk <- seq(0,1,by=0.02)
    }
    # Ensure floating point comparison don't get us
    hsbrk[1] <- -0.001
    hsbrk[length(hsbrk)] <- 1.001
    
    if( is.null(ab1) & is.null(ab2) ){
      # Null marginal histogram
      graphics::barplot(height=0.01, width=0.02,  space = 0, horiz = T, add = FALSE, col="#000000", xlim = c(0,1.0))
    }
    
    if ( !is.null(ab1) & is.null(ab2) ){
      bp1 <- graphics::hist(ab1[,mySamp], breaks = hsbrk, plot = FALSE)
      graphics::barplot(height=bp1$counts, width=fp1$bin_width,  space = 0, horiz = T, add = FALSE, col=col1)
    }
    if ( is.null(ab1) & !is.null(ab2) ){
      bp2 <- graphics::hist(ab2[,mySamp], breaks = hsbrk, plot = FALSE)
      graphics::barplot(height=bp2$counts, width=fp2$bin_width,  space = 0, horiz = T, add = FALSE, col=col2)
    }
    if ( !is.null(ab1) & !is.null(ab2) ){
      bp1 <- graphics::hist(ab1[,mySamp], breaks = hsbrk, plot = FALSE)
      graphics::barplot(height=bp1$counts, width=fp1$bin_width,  space = 0, horiz = T, add = FALSE, col=col1)
      bp2 <- graphics::hist(ab2[,mySamp], breaks = hsbrk, plot = FALSE)
      graphics::barplot(height=bp2$counts, width=fp2$bin_width,  space = 0, horiz = T, add = TRUE, col=col2)
    }
    graphics::title(xlab="Count")
    graphics::par(mar=c(5,4,4,2))
  }
  
  if(layout == TRUE){
    graphics::par(orig_par)
  }
  return( invisible(NULL) )
}


