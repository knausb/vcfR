#'
#'
#' @rdname chromo_plot2
#' @title chromo plot v2
#' @name  chromo plot v2
# @aliases dot.plot
#'
#' @description Plot chromR objects and their components
#' 
#' @param mat a numeric matrix where column one represents the chromosomal position (POS) and subsequent columns represent y-values.
#' @param lst a list containing numeric matrices where each first column is the start and each second column is the end coordinates for rectangles.
#' @param title a title for the plot.
#' @param scale logical specifying whether to scale bars to one.
#' @param col a vector of colors for each column of y-values, recycled as necessary.
#' @param hist logical specifying whether marginal histograms should be plotted.
#' @param layout logical specifying whether to call layout or not.
#' @param mwidth numeric specifying the relative width of the main panel when using layout.
#' @param hline numeric vector specifying positions for horizontal lines to be drawn.
#' @param heights a numeric vector of heights for rectangles, recycled as necessary.
#' @param xmin minimum value for rectangle plots.
#' @param xmax maximum value for rectangle plots.
#' @param ... arguments to be passed to/from other methods.
#' 
#' 
#' @details Plot chromR objects and their components.
#' 
#' The function \strong{dot.plot}
#'
#' The function \strong{bar.plot}
#'
#' The function \strong{rect.plot}
#' 


#' @rdname chromo_plot2
#' 
#' 
#' @export
dot.plot <- function( mat, title = NULL, hline = NULL, col = NULL, hist = TRUE, layout = TRUE, mwidth = 4, ... ){
  
  xmin <- min( mat[, 1], na.rm =TRUE )
  xmax <- max( mat[, 1], na.rm =TRUE )
  ymin <- min( mat[,-1], na.rm =TRUE )
  ymax <- max( mat[,-1], na.rm =TRUE )
  
  if( is.null(col) ){
    col <- 1:8
  }
  
  if( layout == TRUE & hist == TRUE ){
    layout( matrix( 1:2, nrow=1, ncol=2), widths = c(mwidth,1) )
  }
  
  org.mar <- par("mar")
  par( mar=c(0,4,0,0) )
  
  plot( x=c(xmin, xmax), y=c(ymin, ymax), type = "n", xlab="", xaxt="n", ylab="", las=2)

  if( !is.null(hline) ){
    abline( h = hline, lty=2, col="#808080" )
  }
  
  for( i  in 2:ncol(mat) ){
    points( mat[,1], mat[,i], pch=20, col=col[i])
  }
  title( main = title, line = -1)
  
  if( hist == TRUE ){
    par( mar=c(0,0,0,0) )
    boxplot( mat[,-1], xlab="", xaxt="n", ylab="", yaxt="n", col= c(col + 1) )
  }

  if( layout == TRUE & hist == TRUE ){
    par( mfrow = c(1,1) )
  }
  
  par( mar=org.mar)
}


#' @rdname chromo_plot2
#' 
#' 
#' @export
bar.plot <- function( mat, scale = FALSE, title = NULL, col = NULL, hist = TRUE, layout = TRUE, mwidth = 4, ... ){
  
  if( scale == TRUE ){
    sums <- rowSums(mat[,-1], na.rm = TRUE)
    mat[,-1] <- sweep(mat[,-1], MARGIN=1, STATS=sums, FUN="/")
  }
#  ymin <- min( mat[,-1], na.rm =TRUE )
#  ymin <- 0
  if( scale == TRUE ){
    ymax <- 1
  } else {
    ymax <- sum(apply(mat[,-1], MARGIN=2, max, na.rm = TRUE), na.rm = TRUE)
  }

  if( is.null(col) ){
    col <- 1:8
  }

  if( layout == TRUE & hist == TRUE ){
    layout( matrix( 1:2, nrow=1, ncol=2), widths = c(mwidth,1) )
  }
  
  org.mar <- par("mar")
  par( mar=c(0,4,0,0) )

  barplot( t(mat[,-1]), xlim=c(0,nrow(mat)), 
           space=0, col=c(col+1), border=NA,
           las=2,
           pty="m" )
  title( main = title, line = -1)
  
  if( hist == TRUE ){
    par( mar=c(0,0,0,0) )
    boxplot( mat[,-1], xlab="", xaxt="n", ylab="", yaxt="n", col= c(col + 1), ylim=c(0,ymax) )
  }

  if( layout == TRUE & hist == TRUE ){
    par( mfrow = c(1,1) )
  }
  par( mar=org.mar)
}


#' @rdname chromo_plot2
#' 
#' 
#' @export
rect.plot <- function( lst, heights = 1, xmin = 0, 
                       xmax = NULL, title = NULL, 
                       col = NULL,
                       ... ){
  
  if( is.null(xmax) ){
    stop("xmax not specified.")
  }
  if( is.null(col) ){
    col <- 1:8
  }

  ymax <- max(heights)
  ymin <- -1 * ymax
  
  org.mar <- par("mar")
  par( mar=c(0,4,0,0) )
  
  plot( x=c(xmin, xmax), y=c(ymin, ymax),
        type = "n", xlab="", xaxt="n",
        ylab="", yaxt="n", las=2, frame.plot = FALSE )

  lines( x = c(xmin, xmax), y = c(0,0))
  
  for( i in 1:length(lst) ){
    rect( xleft=lst[[i]][,1], 
          ybottom= c(-1 * heights[i]), 
          xright=lst[[i]][,2], 
          ytop=heights[i],
          col = col[i],
          border = NA )
  }
  
  title( main = title, line = -1)
  
  par( mar=org.mar)
}


#' @rdname chromo_plot2
#' 
#' 
#' @export
null.plot <- function(){
  org.mar <- par("mar")
  par( mar=c(0,4,0,0) )
  
  plot( 1:10, 1:10, type = "n", axes = FALSE, frame.plot = FALSE, xlab="", ylab="" )
  
  par( mar=org.mar)
}

# EOF.