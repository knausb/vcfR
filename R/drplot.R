
#' @title dr plot elements
#' @name  dr plot elements
#'
#' @description Plot chromR objects and their components
#' @rdname drplot
#' 
#' @param dmat a numeric matrix for dot plots where the first column is position (POS) and subsequent columns are y-values.
#' @param chrom.s start position for the chromosome
#' @param chrom.e end position for the chromosome
#' @param rlst a list containing numeric matrices containing rectangle coordinates.
#' @param title optional string to be used for the plot title.
#' @param hline vector of positions to be used for horizontal lines.
#' @param dcol vector of colors to be used for dot plots.
#' @param rcol vector of colors to be used for rectangle plots.
#' @param rbcol vector of colors to be used for rectangle borders.
#' @param ... arguments to be passed to other methods.
#' 
#' 
#' @details Plot details
#' The parameter \strong{rlist} is list of numeric matrices containing rectangle coordinates.
#' The first column of each matrix is the left positions, the second column is the bottom coordinates, the third column is the right coordinates and the fourth column is the top coordinates.
#' 
#' @return Returns the y-axis minimum and maximum values invisibly.
#' 
#' @seealso 
#' \code{\link[graphics]{rect}}
#' \code{\link{chromo}}
#' 
#' @export
dr.plot <- function( dmat = NULL, rlst = NULL,
                     chrom.s = 1, chrom.e = NULL,
                     title = NULL, hline = NULL,
                     dcol = NULL, 
                     rcol = NULL, rbcol = NULL,
                     ... ){
  
  # Attempt to handle rlst
  if( class(rlst) != "list" & class(rlst) == "matrix" ){
    rlst <- list( rlst )
  }
  if( class(rlst) != "list" & !is.null(rlst) ){
    stop( paste("parameter rlst is of type", class(rlst), "instead of type list.") )
  }
  
  
  # Determine x max.
  if( is.null(chrom.e) ){
    stop("chrom.e (end chromosome position) must be specified.")
  }
    
  # Determine y min and max.
  if( !is.null(dmat) ){
    ymin <- min( dmat[,-1], na.rm = TRUE )
    ymax <- max( dmat[,-1], na.rm = TRUE )
  } else {
    ymin <-  0.0
    ymax <-  0.1
  }
  if( !is.null(rlst) ){
    rmin <- min( unlist( lapply( rlst, function(x){ min(x[,c(2,4)], na.rm = TRUE) } ) ) )
    rmax <- max( unlist( lapply( rlst, function(x){ max(x[,c(2,4)], na.rm = TRUE) } ) ) )
  } else {
    rmin <- 0
    rmax <- 0    
  }
  ymin <- min( c(ymin, rmin), na.rm = TRUE )
  ymax <- max( c(ymax, rmax), na.rm = TRUE )
  ymin <-  ymin * 1.1
  ymax <-  ymax * 1.1
#  if( ymin == 0 ){ ymin <- -0.9 }
  
  # Color palettes.
  if( is.null(dcol) ){
    dcol <- 1:8
  }
  if( is.null(rcol) ){
    rcol <- 1:8
  }
  if( is.null(rbcol) ){
    rbcol <- 1:8
  }
  
  # Initialize the plot.
  plot( c(chrom.s, chrom.e), c(0,0), type="n",
        xaxt = "n", xlab="", 
        ylab="", ylim = c(ymin, ymax),
        las=1,
        ... )
  
  # Horizontal lines
  if( !is.null(hline) ){
    graphics::abline( h = hline, lty = 2, col = "#808080" )
  }
  
  # Rect plot.
  if( length(rlst) > 1 ){
    rcol  <- rep( rcol, times=length(rlst))
    rbcol <- rep(rbcol, times=length(rlst))
  }
  if( !is.null(rlst) ){
    for( i in 1:length(rlst) ){
      rmat <- rlst[[i]]
      graphics::rect( xleft = rmat[,1], ybottom = rmat[,2], 
            xright = rmat[,3], ytop = rmat[,4],
            col = rcol[i], border = rbcol[i],
            ... )
    }
  }
#  palette("default")
  
  # Dot plot.
  if( !is.null(ncol(dmat)) ){
    dcol  <- rep( dcol, times=ncol(dmat) )
  }
  if( !is.null(dmat) ){
    POS <- dmat[ ,1 ]
    dmat <- dmat[ ,-1 , drop=FALSE ]
    for( i in 1:ncol(dmat) ){
      graphics::points( POS, dmat[,i], pch = 20, col = dcol[i] )
    }
  }
  
  graphics::title( main = title, line = -1.2 )
  
  return( invisible( c(ymin, ymax) ) )
}


#' @rdname drplot
#' 
#' 
#' @export
null.plot <- function(){
  org.mar <- graphics::par("mar")
  graphics::par( mar=c(0,4,0,0) )
  
  plot( 1:10, 1:10, type = "n", axes = FALSE, frame.plot = FALSE, xlab="", ylab="" )
  
  graphics::par( mar=org.mar)
}


##### ##### ##### ##### #####
# EOF.