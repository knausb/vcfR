

#' @rdname chromo_plot
#' @title Plot chromR object
#' @name chromo_plot
#' @export
#' @aliases chromo
#' 
#' @param chrom an object of class chrom.
#' @param boxp logical specifying whether marginal boxplots should be plotted.
#' @param alpha degree of transparency applied to points in dot plots.
#' @param ... arguments to be passed to other methods.
#' 
chromo <- function( chrom,
#                   verbose=TRUE,
#                   nsum=TRUE,
                   ...){
}


#' @rdname chromo_plot
#' @export
#' @aliases chromoqc
#'
chromoqc <- function( chrom, boxp = TRUE, alpha = 255, ...){
  org.par <- par( no.readonly = TRUE )

  # Set parameters.
  mwidth <- 8
  ncols <- 1
  nrows <- 3
  heights <- c(1,1,1)
  par( oma = c(3,0,1,0) )

  
  # Add to parameters.
  if( boxp == TRUE ){
    ncols <- ncols + 1
    mwidth <- c(mwidth, 1)
  }
  if( !is.null(chrom@win.info$variants) ){
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }
  if( !is.null(chrom@seq) ){ 
    nrows <- nrows + 2
    heights <- c( heights, 1, 0.25)
  }
  if( nrow(chrom@ann) > 0 ){
    nrows <- nrows + 1
    heights <- c( heights, 0.25)    
  }
  
  # Plot
  chrom.par <- par( no.readonly = TRUE )
  layout( matrix( 1:c( ncols * nrows ), nrow=nrows, ncol=ncols, byrow = TRUE ),
          widths = mwidth, heights = heights
  )
  dmat <- as.matrix( cbind(chrom@var.info[,"POS"], chrom@var.info[,"DP"]) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  dot.plot( dmat,
            title="Read Depth (DP)", mwidth=mwidth,
            layout = FALSE, boxp = boxp,
            col=rgb( red=30, green=144, blue=255, alpha=alpha, maxColorValue = 255),
            ... )
  dmat <- as.matrix( cbind(chrom@var.info[,"POS"], chrom@var.info[,"MQ"]) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  dot.plot( dmat,
            title="Mapping Quality (MQ)", mwidth=mwidth,
            layout=F, boxp = boxp,
            col=rgb( red=46, green=139, blue=87, alpha=alpha, maxColorValue = 255),
            ... )
  dmat <- as.matrix( cbind(chrom@var.info[,"POS"], 
                           as.numeric( chrom@vcf@fix[,"QUAL"] ) ) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  dot.plot( dmat,
            title="Phred-Scaled Quality (QUAL)", mwidth=mwidth,
            layout=F, boxp = boxp,
            col=rgb( red=139, green=0, blue=139, alpha=alpha, maxColorValue = 255),
            ... )
  if( !is.null(chrom@win.info$variants) ){
    bmat <- cbind( chrom@win.info[,"window"], 
                   chrom@win.info[,"variants"]/c(chrom@win.info[,"end"]-chrom@win.info[,"start"]+1) )
    bar.plot( bmat,
              title="Variants", hline = seq(0,0.1,by=0.02),
              layout = FALSE, scale = FALSE,
              col=rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255),
              boxp = boxp,
              ... )
  }
#  if( !is.null(chrom@win.info[,"A"]) ){
  if( length( grep("^A$", colnames(chrom@win.info)) ) == 1 ){  
    bmat <- cbind( chrom@win.info[,"window"], 
                   chrom@win.info[,"A"] + chrom@win.info[,"T"],
                   chrom@win.info[,"C"] + chrom@win.info[,"G"])
    bmat[,2] <- bmat[,2] / c(chrom@win.info[,"end"]-chrom@win.info[,"start"]+1)
    bmat[,3] <- bmat[,3] / c(chrom@win.info[,"end"]-chrom@win.info[,"start"]+1)
    bar.plot( bmat,
              title="Nucleotide content",
              layout = FALSE, scale = FALSE,
              col=c(rgb( red=000, green=034, blue=205, alpha=255, maxColorValue = 255),
                    rgb( red=255, green=235, blue=000, alpha=255, maxColorValue = 255)
                    ),
              boxp = boxp,
              ... )
    rect.plot(chrom@seq.info, xmax=chrom@len,
              title="Nucleotides", 
              heights = c(1, 0.5), 
              col=c('green', 'red'),
              ... )
    if( boxp == TRUE){
      null.plot()
    }
  }
  if( nrow(chrom@ann) > 0 ){
    rect.plot( list(chrom@ann[,4:5]),
               xmax=chrom@len,
               title = "Annotations",
               col = rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255),
               ... )
  }

#  axis( side = 1, line = 5.2 )
#  axis( side = 1, line = 0 )
#  title( xlab = "Base pairs", line = 1.6, outer = TRUE )

  title( main = chrom@name, line = 0.2, outer = TRUE )
  # identify, locator

  par( org.par )
}


#' @rdname chromo_plot
#' @export
#' @aliases chromoqc
#'
chromohwe <- function(chrom, ...){
  stop("Function not implemented.")
}


#' @rdname chromo_plot
#' @export
#' @aliases chromodot
#'
chromodot <- function( chrom, ...){
  message("Function not implemented.")
}


#' @rdname chromo_plot
#' @export
#' @aliases chromopop
#'
chromopop <- function( chrom, ...){
  stop("Function not implemented.")
}





##### ##### ##### ##### #####
# EOF.