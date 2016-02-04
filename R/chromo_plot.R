

#' @rdname chromo_plot
#' @title Plot chromR object
#' @name chromo_plot
#' @export
#' @aliases chromo
#' 
#' @param chrom an object of class chrom.
#' @param boxp logical specifying whether marginal boxplots should be plotted [T/F].
#' @param dp.alpha degree of transparency applied to points in dot plots [0-255].
#' @param chrom.s	start position for the chromosome.
#' @param chrom.e	end position for the chromosome.
#' @param drlist1 a named list containing elements to create a drplot
#' @param drlist2 a named list containing elements to create a drplot
#' @param drlist3 a named list containing elements to create a drplot
#' @param ... arguments to be passed to other methods.
#' 
#' 
#' @details 
#' Each \strong{drlist} parameter is a list containing elements necessarry to plot a dr.plot.
#' This list should contain up to seven elements named title, dmat, rlist, dcol, rcol, rbcol and bwcol.
#' These elements are documented in the dr.plot page where they are presented as individual parameters.
#' The one exception is bwcol which is a vector of colors for the marginal box and whisker plot.
#' This is provided so that different colors may be used in the dot plot and the box and whisker plot.
#' For example, transparency may be desired in the dot plot but not the box and whisker plot.
#' When one (or more) of these elements is omitted an attempt to use default values is made.
#' 
#' 
#' 
#' @seealso \code{\link{dr.plot}}
#' 
#' 
chromo <- function( chrom,
                    boxp = TRUE,
                    dp.alpha = TRUE,
                    chrom.s = 1,
                    chrom.e = NULL,
                    drlist1 = NULL, drlist2 = NULL, drlist3 = NULL,
#                    title1, dmat1, rlist1, dcol1, rcol1, rbcol1,
#                    title2, dmat2, rlist2, dcol2, rcol2, rbcol2,
#                    title3, dmat3, rlist3, dcol3, rcol3, rbcol3,
#                   verbose=TRUE,
#                   nsum=TRUE,
                   ...){
  
  if( class(chrom) != "chromR" ){
    stop("Expecting object of class chromR")
  }
  
  # Save original parameters.
  orig.oma <- graphics::par('oma')
  orig.mar <- graphics::par('mar')
    
  # Initialize parameters.
  mwidth <- 8
  ncols <- 1
  nrows <- 0
  mheight <- 0.3 # Minor plot height
  heights <- c()

  ##### ##### ##### ##### #####
  #  
  # Determine the layout of the plot.
  #
  ##### ##### ##### ##### #####
  
  # Plot title
  if( length(chrom@name) > 0 ){
    graphics::par( oma = c(3,0,1,0) )
  } else {
    graphics::par( oma = c(3,0,0,0) )
  }
  
  # Marginal boxplots.
  if( boxp == TRUE ){
    ncols <- ncols + 1
    mwidth <- c(mwidth, 1)
  }

  # drlist1
  if( !is.null(drlist1) ){
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }
  # drlist2
  if( !is.null(drlist2) ){
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }
  # drlist3
  if( !is.null(drlist3) ){
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }

  # Variant plot.
  if( !is.null(chrom@win.info$variants) ){
    nrows <- nrows + 1
    heights <- c( heights, 1)
  }
  
  # Nucleotide content and sequence plots.
  if( length( grep("^A$", colnames(chrom@win.info) ) ) == 1 ){
    nrows <- nrows + 2
    heights <- c( heights, 1, mheight)
  }
  
  # Annotation plot.  
  if( nrow(chrom@ann) > 0 ){
    nrows <- nrows + 1
    heights <- c( heights, mheight)
  }
  

  if( nrows == 0 ){
    stop("no data has been included!")
  }

  
  ##### ##### ##### ##### #####
  #
  # Establish layout for plot
  #
  ##### ##### ##### ##### #####
  
  graphics::layout( matrix( 1:c( ncols * nrows ),
                  nrow=nrows, 
                  ncol=ncols, 
                  byrow = TRUE ),
          widths = mwidth,
          heights = heights
  )


  ##### ##### ##### ##### #####
  #
  # Plot
  #
  ##### ##### ##### ##### #####  

  ##### ##### ##### ##### #####
  #
  # drplots
  #
  ##### ##### ##### ##### #####  

  # drplot1
  if( !is.null(drlist1) ){
    graphics::par( mar = c(0,4,0,0) )
  
    bdim <- dr.plot( dmat    = drlist1$dmat,
                     rlst    = drlist1$rlst,
                     chrom.s = chrom.s,
                     chrom.e = chrom.e,
                     title   = drlist1$title,
                     dcol    = drlist1$dcol,
                     rcol    = drlist1$rcol,
                     rbcol   = drlist1$rbcol,
                     ... )
    graphics::par( mar = orig.mar )
  
    if( boxp == TRUE ){
      graphics::par( mar = c(0,0,0,0) )
      if( is.null(drlist1$bwcol) ){
        drlist1$bwcol <- drlist1$dcol
      }
      graphics::boxplot( x    = drlist1$dmat[,-1],
               ylim = bdim,
               yaxt = "n",
               col  = drlist1$bwcol
              )
      graphics::par( mar = orig.mar )
    }
  }
  
  # drplot2
  if( !is.null(drlist2) ){
    graphics::par( mar = c(0,4,0,0) )
  
    bdim <- dr.plot( dmat    = drlist2$dmat,
                     rlst    = drlist2$rlst,
                     chrom.s = chrom.s,
                     chrom.e = chrom.e,
                     title   = drlist2$title,
                     dcol    = drlist2$dcol,
                     rcol    = drlist2$rcol,
                     rbcol   = drlist2$rbcol,
                     ... )
    graphics::par( mar = orig.mar )
  
    if( boxp == TRUE ){
      graphics::par( mar = c(0,0,0,0) )
      if( is.null(drlist2$bwcol) ){
        drlist2$bwcol <- drlist2$dcol
      }
      graphics::boxplot( x    = drlist2$dmat[,-1],
               ylim = bdim,
               yaxt = "n",
               col  = drlist2$bwcol
              )
      graphics::par( mar = orig.mar )
    }
  }

  # drplot3
  if( !is.null(drlist3) ){
    graphics::par( mar = c(0,4,0,0) )
  
    bdim <- dr.plot( dmat    = drlist3$dmat,
                     rlst    = drlist3$rlst,
                     chrom.s = chrom.s,
                     chrom.e = chrom.e,
                     title   = drlist3$title,
                     dcol    = drlist3$dcol,
                     rcol    = drlist3$rcol,
                     rbcol   = drlist3$rbcol,
                     ... )
    graphics::par( mar = orig.mar )
  
    if( boxp == TRUE ){
      graphics::par( mar = c(0,0,0,0) )
      if( is.null(drlist3$bwcol) ){
        drlist3$bwcol <- drlist3$dcol
      }
      graphics::boxplot( x    = drlist3$dmat[,-1],
               ylim = bdim,
               yaxt = "n",
               col  = drlist3$bwcol
              )
      graphics::par( mar = orig.mar )
    }
  }

  ##### ##### ##### ##### #####
  #
  # chromR plots
  #
  ##### ##### ##### ##### #####  

  # Variant plot
  if( !is.null(chrom@win.info$variants) ){
    rmat <- cbind(chrom@win.info[,'start'] ,
                  0,
                  chrom@win.info[,'end'],
                  chrom@win.info[,'variants'] / c(chrom@win.info[,'end'] - chrom@win.info[,'start'])
    )
    
    graphics::par( mar = c(0,4,0,0) )
    bdim <- dr.plot( dmat = NULL, rlst = list( rmat ), chrom.s = 1, chrom.e = chrom@len,
                     title = "Variants per Site", hline = NULL,
                     dcol = NULL,
                     rcol  = grDevices::rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ), 
                     rbcol = grDevices::rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ),
                     ... )
    if( length( grep("^A$", colnames(chrom@win.info)) ) == 0 & nrow(chrom@ann) == 0 ){
      graphics::axis( side = 1, line = 0 )      
    }
    graphics::par( mar = c(5,4,4,2) + 0.1 )
    
    if( boxp == TRUE ){
      graphics::par( mar = c(0,0,0,0) )
      graphics::boxplot( rmat[,4], ylim=bdim, yaxt = "n",
               col = grDevices::rgb( red=178, green=34, blue=34, alpha=255, maxColorValue = 255 ) )
      graphics::par( mar = c(5,4,4,2) + 0.1 )
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
    
    graphics::par( mar = c(0,4,0,0) )
    bdim <- dr.plot( dmat = NULL, rlst = list(rmat1, rmat2), chrom.s = 1, chrom.e = chrom@len,
                     title = "Nucleotide Content", hline = NULL,
                     dcol = NULL,
                     rcol = c(grDevices::rgb( red=000, green=034, blue=205, maxColorValue = 255), 
                              grDevices::rgb( red=255, green=235, blue=000, maxColorValue = 255)),
                     rbcol = c(grDevices::rgb( red=000, green=034, blue=205, maxColorValue = 255), 
                              grDevices::rgb( red=255, green=235, blue=000, maxColorValue = 255)),
                     ... )
    graphics::par( mar = c(5,4,4,2) + 0.1 )

    if( boxp == TRUE ){
      graphics::par( mar = c(0,0,0,0) )
      rmat1 <- cbind(rmat1[,4], rmat2[,4])
      rmat1[,2] <- rmat1[,2] - rmat1[,1]
      graphics::boxplot( rmat1, ylim=bdim, yaxt = "n",
               col = c(grDevices::rgb( red=000, green=034, blue=205, maxColorValue = 255),
                       grDevices::rgb( red=255, green=235, blue=000, maxColorValue = 255)
                      ),
               border = c(grDevices::rgb( red=000, green=034, blue=205, maxColorValue = 255),
                          grDevices::rgb( red=255, green=235, blue=000, maxColorValue = 255)
                          ),
               xaxt = "n"
      )
      graphics::par( mar = c(5,4,4,2) + 0.1 )
    }
    
    # Sequence plot.
    rmat1 <- cbind(chrom@seq.info$nuc.win[,1], -1, chrom@seq.info$nuc.win[,2], 1)
    rmat2 <- cbind(chrom@seq.info$N.win[,1], -0.5, chrom@seq.info$N.win[,2], 0.5)
    graphics::par( mar = c(0,4,0,0) )
    dr.plot( rlst = list( rmat1, rmat2 ), chrom.s = 1, chrom.e = chrom@len,
                    title = "Nucleotides", hline = NULL,
                    dcol = NULL,
                    rcol  = c('green', 'red'),
                    rbcol = c('green', 'red'),
                    yaxt = "n",
                    #frame.plot = FALSE,
                    ... )
    
    if( nrow(chrom@ann) == 0 ){
      graphics::axis( side = 1, line = 0 )
    }
    graphics::par( mar = c(5,4,4,2) + 0.1 )
    
    if( boxp == TRUE){
      null.plot()
    }
  }
  
  # Annotation plot.
  if( nrow(chrom@ann) > 0 ){
    rmat <- cbind( chrom@ann[,4], -1, chrom@ann[,5], 1)
    graphics::par( mar=c(0,4,0,0) )
    dr.plot( rlst = list( rmat ), chrom.e = chrom@len, title = "Annotations",
             rcol = grDevices::rgb(178,34,34, maxColorValue = 255),
             rbcol = grDevices::rgb(178,34,34, maxColorValue = 255),
             hline = 0,
             yaxt = "n",
             ...)
    graphics::axis( side = 1, line = 0 )
    graphics::par( mar = c(5,4,4,2) + 0.1 )

    if( boxp == TRUE){
      null.plot()
    }
  }

  graphics::title( xlab = "Base pairs", line = 1.6, outer = TRUE )
  if( length(chrom@name) > 0 ){
    graphics::title( main = chrom@name, line = 0.2, outer = TRUE )
  }

  
  ##### ##### ##### ##### #####
  #  
  # Reset graphics parameters to defaults.
  #
  ##### ##### ##### ##### #####

  graphics::par( mar = orig.mar )
  graphics::par( oma = orig.oma )
  graphics::par( mfrow = c(1,1) )
}

##### ##### ##### ##### #####
#  
# End chromo
#
##### ##### ##### ##### #####

##### ##### ##### ##### #####
#  
# Begin chromoqc
#
##### ##### ##### ##### #####

#' @rdname chromo_plot
#' @export
#' @aliases chromoqc
#'
chromoqc <- function( chrom, 
                      boxp = TRUE, 
                      dp.alpha = 255,
                      ...){
  
  if( class(chrom) != "chromR" ){
    stop( paste("expecting an object of class chromR, got", class(chrom), "instead.") )
  }
  
  # Read depth
  myList1 <- list(title = "Read Depth (DP)",
                  dmat  = chrom@var.info[ chrom@var.info[,"mask"] , c("POS","DP") ],
                  dcol  = grDevices::rgb( red=30, green=144, blue=255, alpha=dp.alpha, maxColorValue = 255),
                  bwcol = grDevices::rgb( red=30, green=144, blue=255, maxColorValue = 255)
  )

  # Mapping Quality (MQ)
  myList2 <- list(title = "Mapping Quality (MQ)",
                  dmat  = chrom@var.info[ chrom@var.info[,"mask"] , c("POS","MQ") ],
                  dcol  = grDevices::rgb( red=46, green=139, blue=87, alpha=dp.alpha, maxColorValue = 255),
                  bwcol = grDevices::rgb( red=46, green=139, blue=87, maxColorValue = 255)
  )
  
  # Phred-Scaled Quality (QUAL)
  dmat <- as.matrix( cbind(chrom@var.info[,"POS"], 
                           as.numeric( chrom@vcf@fix[,"QUAL"] ) ) )
  dmat <- dmat[ chrom@var.info[,"mask"], , drop = FALSE]
  myList3 <- list(title = "Phred-Scaled Quality (QUAL)",
                  dmat  = dmat,
                  dcol  = grDevices::rgb(red=139, green=0, blue=139, alpha=dp.alpha, maxColorValue = 255),
                  bwcol = grDevices::rgb(red=139, green=0, blue=139, maxColorValue = 255)
  )
  
  chromo( chrom, boxp = boxp, 
          chrom.e = chrom@len, 
          drlist1 = myList1,
          drlist2 = myList2,
          drlist3 = myList3,
          ...
  )
}


##### ##### ##### ##### #####
#  
# End chromoqc
#
##### ##### ##### ##### #####




# ' @rdname chromo_plot
# ' @export
# ' @aliases chromoqc
# '
#chromohwe <- function(chrom, ...){
#  stop("Function not implemented.")
#}

# ' @rdname chromo_plot
# ' @export
# ' @aliases chromodot
# '
#chromodot <- function( chrom, ...){
#  message("Function not implemented.")
#}

# ' @rdname chromo_plot
# ' @export
# ' @aliases chromopop
# '
#chromopop <- function( chrom, ...){
#  stop("Function not implemented.")
#}

##### ##### ##### ##### #####
# EOF.