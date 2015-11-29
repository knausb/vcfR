#'
#'
#' @rdname chromo_elements
#' @title chromo plot elements
#' @name  chromo plot elements
# @aliases dot.plot
#'
#' @description Plot chromR objects and their components
#' 
#' @param mat a numeric matrix where column one represents the chromosomal position (POS) and subsequent columns represent y-values.
#' @param lst a list containing numeric matrices where each first column is the start and each second column is the end coordinates for rectangles.
#' @param title a title for the plot.
#' @param scale logical specifying whether to scale bars to one.
#' @param col a vector of colors for each column of y-values, recycled as necessary.
#' @param boxp logical specifying whether marginal boxplots should be plotted.
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


#' @rdname chromo_elements
#' 
#' 
#' @export
dot.plot <- function( mat, title = NULL, hline = NULL, col = NULL, boxp = TRUE, layout = TRUE, mwidth = 4, ... ){
#  org.par <- par( no.readonly = TRUE )
  org.mar <- par("mar")
  
  xmin <- min( mat[, 1], na.rm =TRUE )
  xmax <- max( mat[, 1], na.rm =TRUE )
  ymin <- min( mat[,-1], na.rm =TRUE )
  ymax <- max( mat[,-1], na.rm =TRUE )
  
  POS <- mat[ ,  1 ]
  mat <- mat[ , -1, drop = FALSE ]
  
  if( is.null(col) ){
    col <- 1:8
  }
  
  if( layout == TRUE & boxp == TRUE ){
    layout( matrix( 1:2, nrow=1, ncol=2), widths = c(mwidth,1) )
  }
  
  par( mar=c(0,4,0,0) )
  
  plot( x=c(xmin, xmax), 
        y=c(ymin, ymax), 
        type = "n", xlab="",
        xaxt="n", ylab="", las=2)

  if( !is.null(hline) ){
    abline( h = hline, lty=2, col="#808080" )
  }
  
  for( i  in 1:ncol(mat) ){
    points( POS, mat[,i], pch=20, col=col[i])
  }
  title( main = title, line = -1)
  
  if( boxp == TRUE ){
    par( mar=c(0,0,0,0) )
    col2 <- substr( col, start = 1, stop = 7)
    boxplot( mat, xlab="", xaxt="n", ylab="", yaxt="n", col= col2 )
  }

  if( layout == TRUE & boxp == TRUE ){
    par( mfrow = c(1,1) )
  }
  
#  par( org.par )
  par( mar = org.mar )
}


#' @rdname chromo_elements
#' 
#' 
#' @export
bar.plot <- function( mat, scale = FALSE, title = NULL, 
                      hline = NULL, col = NULL, boxp = TRUE,
                      layout = TRUE, mwidth = 4, ... ){
#  org.par <- par( no.readonly = TRUE )
  org.mar <- par( "mar" )
  POS <- mat[ ,  1 ]
  mat <- mat[ , -1, drop = FALSE ]
  
  if( scale == TRUE ){
    sums <- rowSums(mat, na.rm = TRUE)
    mat <- sweep(mat, MARGIN=1, STATS=sums, FUN="/")
  }
  if( scale == TRUE ){
    ymax <- 1
  } else {
    ymax <- sum(apply(mat, MARGIN=2, max, na.rm = TRUE), na.rm = TRUE)
  }

  if( is.null(col) ){
    col <- 1:8
  }

  if( layout == TRUE & boxp == TRUE ){
    layout( matrix( 1:2, nrow=1, ncol=2), widths = c(mwidth,1) )
  }

  par( mar=c(0,4,0,0) )

  barplot( t(mat), xlim=c(0,nrow(mat)), 
           space=0, col=col, border=NA,
           las=2,
           pty="m" )
  title( main = title, line = -1)
  if( !is.null(hline) ){
    abline( h = hline, lty=2, col="#808080" )
  }
  
  if( boxp == TRUE ){
    par( mar=c(0,0,0,0) )
    boxplot( mat, xlab="", xaxt="n", ylab="", yaxt="n", col= col, ylim=c(0,ymax) )
  }

  if( layout == TRUE & boxp == TRUE ){
    par( mfrow = c(1,1) )
  }
#  par( org.par )
  par( mar = org.mar )
}


#' @rdname chromo_elements
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


#' @rdname chromo_elements
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