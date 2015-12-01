

#' @rdname chromo_plot
#' @title Plot chromR object
#' @name chromo_plot
#' @export
#' @aliases chromo
#' 
#' @param chrom an object of class chrom.
#' @param boxp logical specifying whether marginal boxplots should be plotted.
#' @param dp.alpha degree of transparency applied to points in dot plots.
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
chromoqc <- function( chrom, boxp = TRUE, dp.alpha = 255,
                      ...){
  # Save original parameters for later.
#  org.par <- par( no.readonly = TRUE )
  
  # Set parameters.
  mwidth <- 8
  ncols <- 1
  nrows <- 3
  mheight <- 0.3 # Minor plot height
  heights <- c(1,1,1)
  if( length(chrom@name) > 0 ){
    par( oma = c(3,0,1,0) )
  } else {
    par( oma = c(3,0,0,0) )
  }

  # Add to parameters.
  if( boxp == TRUE ){
    ncols <- ncols + 1
    mwidth <- c(mwidth, 1)
  }
  if( !is.null(chrom@win.info$variants) ){
    # Add variant plot.
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }
  if( !is.null(chrom@seq) ){
    # Add nucleotide content and sequence plots.
    nrows <- nrows + 2
    heights <- c( heights, 1, mheight)
  }
  if( nrow(chrom@ann) > 0 ){
    # Add annotation plot.
    nrows <- nrows + 1
    heights <- c( heights, mheight)
  }

  # Plot
  layout( matrix( 1:c( ncols * nrows ), nrow=nrows, ncol=ncols, byrow = TRUE ),
          widths = mwidth, heights = heights
  )
  
  # DP plot
  dmat <- as.matrix( cbind(chrom@var.info[,c("POS", "DP")]) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  par( mar = c(0,4,0,0) )
  bdim <- dr.plot( dmat, rlst = NULL, chrom.s = 1, chrom.e = chrom@len,
                   title = "Read Depth (DP)", 
                   #hline = seq(0,1e5, by=1e3),
                   dcol = rgb( red=30, green=144, blue=255, alpha=dp.alpha, maxColorValue = 255),
                   rcol = NULL, rbcol = NULL,
                   ... )
  par( mar = c(5,4,4,2) + 0.1 )
  
  if( boxp == TRUE ){
    par( mar = c(0,0,0,0) )
    boxplot( chrom@var.info[,"DP"], ylim=bdim, yaxt = "n",
             col= rgb( red=30, green=144, blue=255, maxColorValue = 255))
    par( mar = c(5,4,4,2) + 0.1 )
  }
  
  # MQ plot
  dmat <- as.matrix( cbind(chrom@var.info[,c("POS", "MQ")]) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  par( mar = c(0,4,0,0) )
  bdim <- dr.plot( dmat, rlst = NULL, chrom.s = 1, chrom.e = chrom@len,
                   title = "Mapping Quality (MQ)", hline = seq(0,1e5, by=1e3),
                   dcol = rgb( red=46, green=139, blue=87, alpha=dp.alpha, maxColorValue = 255),
                   rcol = NULL, rbcol = NULL,
                   ... )
  par( mar = c(5,4,4,2) + 0.1 )
  
  if( boxp == TRUE ){
    par( mar = c(0,0,0,0) )
    boxplot( chrom@var.info[,"MQ"], ylim=bdim, yaxt = "n",
             col = rgb( red=46, green=139, blue=87, maxColorValue = 255) )
    par( mar = c(5,4,4,2) + 0.1 )
  }
  
  # QUAL plot
  dmat <- as.matrix( cbind(chrom@var.info[,"POS"], 
                           as.numeric( chrom@vcf@fix[,"QUAL"] ) ) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  par( mar = c(0,4,0,0) )
  bdim <- dr.plot( dmat, rlst = NULL, chrom.s = 1, chrom.e = chrom@len,
                   title = "Phred-Scaled Quality (QUAL)", hline = NULL,
                   dcol = rgb(red=139, green=0, blue=139, alpha=dp.alpha, maxColorValue = 255),
                   rcol = NULL, rbcol = NULL,
                   ... )
  par( mar = c(5,4,4,2) + 0.1 )
  
  if( boxp == TRUE ){
    par( mar = c(0,0,0,0) )
    boxplot( as.numeric(chrom@vcf@fix[,"QUAL"]), ylim=bdim, yaxt = "n",
             col = rgb(red=139, green=0, blue=139, maxColorValue = 255) )
    par( mar = c(5,4,4,2) + 0.1 )
  }
  
  # Variant plot
  if( !is.null(chrom@win.info$variants) ){
    rmat <- cbind(chrom@win.info[,'start'] ,
                  0,
                  chrom@win.info[,'end'],
                  chrom@win.info[,'variants'] / c(chrom@win.info[,'end'] - chrom@win.info[,'start'])
    )
    
    par( mar = c(0,4,0,0) )
    bdim <- dr.plot( dmat = NULL, rlst = list( rmat ), chrom.s = 1, chrom.e = chrom@len,
                     title = "Variants per Site", hline = NULL,
                     dcol = NULL,
                     rcol  = rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ), 
                     rbcol = rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ),
                     ... )
    if( length( grep("^A$", colnames(chrom@win.info)) ) == 0 & nrow(chrom@ann) == 0 ){
      axis( side = 1, line = 0 )      
    }
    par( mar = c(5,4,4,2) + 0.1 )
    
    if( boxp == TRUE ){
      par( mar = c(0,0,0,0) )
      boxplot( rmat[,4], ylim=bdim, yaxt = "n",
               col = rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ) )
      par( mar = c(5,4,4,2) + 0.1 )
    }
  }
  
  # Nucleotide plot
  if( length( grep("^A$", colnames(chrom@win.info)) ) == 1 ){
    rmat1 <- cbind(chrom@win.info[,'start'],
                   0,
                   chrom@win.info[,'end'],
                   rowSums(chrom@win.info[,c('A', 'T')]) / c(chrom@win.info[,'end'] - chrom@win.info[,'start'] )
    )
    rmat2 <- cbind(chrom@win.info[,'start'],
                   rmat1[,4],
                   chrom@win.info[,'end'],
                   rmat1[,4] + rowSums(chrom@win.info[,c('C', 'G')]) / c(chrom@win.info[,'end'] - chrom@win.info[,'start'] )
    )
    
    par( mar = c(0,4,0,0) )
    bdim <- dr.plot( dmat = NULL, rlst = list(rmat1, rmat2), chrom.s = 1, chrom.e = chrom@len,
                     title = "Nucleotide Content", hline = NULL,
                     dcol = NULL,
                     rcol = c(rgb( red=000, green=034, blue=205, maxColorValue = 255), 
                              rgb( red=255, green=235, blue=000, maxColorValue = 255)),
                     rbcol = c(rgb( red=000, green=034, blue=205, maxColorValue = 255), 
                              rgb( red=255, green=235, blue=000, maxColorValue = 255)),
                     ... )
    par( mar = c(5,4,4,2) + 0.1 )

    if( boxp == TRUE ){
      par( mar = c(0,0,0,0) )
      rmat1 <- cbind(rmat1[,4], rmat2[,4])
      rmat1[,2] <- rmat1[,2] - rmat1[,1]
      boxplot( rmat1, ylim=bdim, yaxt = "n",
               col = c(rgb( red=000, green=034, blue=205, maxColorValue = 255),
                       rgb( red=255, green=235, blue=000, maxColorValue = 255)
                      ),
               border = c(rgb( red=000, green=034, blue=205, maxColorValue = 255),
                          rgb( red=255, green=235, blue=000, maxColorValue = 255)
                          ),
               xaxt = "n"
      )
      par( mar = c(5,4,4,2) + 0.1 )
    }
    
    # Sequence plot.
    rmat1 <- cbind(chrom@seq.info$nuc.win[,1], -1, chrom@seq.info$nuc.win[,2], 1)
    rmat2 <- cbind(chrom@seq.info$N.win[,1], -0.5, chrom@seq.info$N.win[,2], 0.5)
    par( mar = c(0,4,0,0) )
    dr.plot( rlst = list( rmat1, rmat2 ), chrom.s = 1, chrom.e = chrom@len,
                    title = "Nucleotides", hline = NULL,
                    dcol = NULL,
                    rcol  = c('green', 'red'),
                    rbcol = c('green', 'red'),
                    yaxt = "n",
                    #frame.plot = FALSE,
                    ... )
    
    if( nrow(chrom@ann) == 0 ){
      axis( side = 1, line = 0 )
    }
    par( mar = c(5,4,4,2) + 0.1 )
    
    if( boxp == TRUE){
      null.plot()
    }
  }
  
  # Annotation plot.
  if( nrow(chrom@ann) > 0 ){
    rmat <- cbind( chrom@ann[,4], -1, chrom@ann[,5], 1)
    par( mar=c(0,4,0,0) )
    dr.plot( rlst = list( rmat ), chrom.e = chrom@len, title = "Annotations",
             rcol = rgb(178,34,34, maxColorValue = 255),
             rbcol = rgb(178,34,34, maxColorValue = 255),
             hline = 0,
             yaxt = "n",
             ...)
    axis( side = 1, line = 0 )
    par( mar = c(5,4,4,2) + 0.1 )
    

#    title( xlab = "Base pairs", line = 0.2, outer = TRUE )
    
    if( boxp == TRUE){
      null.plot()
    }
    
  }

  title( xlab = "Base pairs", line = 1.6, outer = TRUE )
  if( length(chrom@name) > 0 ){
    title( main = chrom@name, line = 0.2, outer = TRUE )
  }
  
  par( mar = c(5,4,4,2) + 0.1 )
  par( oma = c(0,0,0,0) )
  par( mfrow = c(1,1) )
#  par( org.par )
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