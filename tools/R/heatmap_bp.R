#### Misc functions ####

#' @title Heatmap with barplots
#'
#' @name heatmap.bp
#' @rdname heatmap_bp
#' @aliases heatmap.bp
#' @export
#' 
#' @description
#' Heatmap of a numeric matrix with barplots summarizing columns and rows.
#'
#' 
#' @param x a numeric matrix.
#' @param cbarplot a logical indicating whether the columns should be summarized with a barplot.
#' @param rbarplot a logical indicating whether the rows should be summarized with a barplot.
#' @param legend a logical indicating whether a legend should be plotted.
#' @param clabels a logical indicating whether column labels should be included.
#' @param rlabels a logical indicating whether row labels should be included.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "none".
#' @param col.ramp vector of colors to be used for the color ramp.
#' @param ...  additional arguments to be passed on.
#' 
#' @details The function heatmap.bp creates a heatmap from a numeric matrix with optional barplots to summarize the rows and columns.
#' 
#' @seealso \code{\link[stats]{heatmap}}, \code{\link[graphics]{image}}, heatmap2 in \href{https://cran.r-project.org/package=gplots}{gplots}, \href{https://cran.r-project.org/package=pheatmap}{pheatmap}.
#' 
#' @examples
#' library(vcfR)
#' 
#' x  <- as.matrix(mtcars)
#' 
#' heatmap.bp(x)
#' heatmap.bp(x, scale="col")
#' # Use an alternate color ramp
#' heatmap.bp(x, col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100))
# library(viridis)
#' heatmap.bp(x)
#' 
#' \dontrun{
#' heatmap.bp(x, cbarplot = FALSE, rbarplot = FALSE, legend = FALSE)
#' heatmap.bp(x, cbarplot = FALSE, rbarplot = TRUE, legend = FALSE)
#' heatmap.bp(x, cbarplot = FALSE, rbarplot = FALSE, legend = TRUE)
#' heatmap.bp(x, cbarplot = FALSE, rbarplot = TRUE, legend = TRUE)
#' 
#' heatmap.bp(x, cbarplot = TRUE, rbarplot = FALSE, legend = FALSE)
#' heatmap.bp(x, cbarplot = TRUE, rbarplot = TRUE, legend = FALSE)
#' heatmap.bp(x, cbarplot = TRUE, rbarplot = FALSE, legend = TRUE)
#' heatmap.bp(x, cbarplot = TRUE, rbarplot = TRUE, legend = TRUE)
#' }
#'
#'  
# data(vcfR_example)
# pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
# pinf_mt <- masker(pinf_mt)
# pinf_gq <- extract.gt(pinf_mt, element="GQ", as.numeric=TRUE)
# heatmap.bp(pinf_gq)
# heatmap.bp(pinf_gq, scale="col")
# heatmap.bp(pinf_gq, col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100))
# heatmap.bp(pinf_gq, col.ramp = colorRampPalette(c("#D55E00", "#F0E442", "#009E73"))(100))
#' @importFrom viridisLite viridis
heatmap.bp <- function(x, cbarplot = TRUE, rbarplot = TRUE,
                       legend = TRUE, clabels = TRUE, rlabels = TRUE,
                       na.rm = TRUE, scale = c("row", "column", "none"),
#                       col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100),
#                       col.ramp = viridis::viridis(n = 100, alpha=1),
                       col.ramp = viridisLite::viridis(n = 100, alpha=1),
                       ...){
#  require(viridis)
#  viridisLite::viridis(n=4)
  stopifnot(class(x) == 'matrix')
  scale <- if(missing(scale))
    "none"
  else match.arg(scale)
  #
  
  # Determine the geometry of the plot.
  nrows <- 1
  ncols <- 1
  if(cbarplot == TRUE){ nrows <- nrows + 1 }
  if(rbarplot == TRUE){ ncols <- ncols + 1 }
  if(legend == TRUE){ ncols <- ncols + 1 }
  
  # Scale the data as appropriate.
  if (scale == "row") {
    x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 1L, stats::sd, na.rm = na.rm)
    x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
  }
  else if (scale == "column") {
    x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
    sx <- apply(x, 2L, stats::sd, na.rm = na.rm)
    x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
  }
  
  # Handle column names.
  if(is.null(colnames(x))){colnames(x) <- 1:ncol(x)}
  if(is.null(rownames(x))){rownames(x) <- 1:nrow(x)}
  
  
  # Get user's par(), ignoring the read-only variables.
  userpar <- graphics::par(no.readonly = TRUE)
  # Promise to reset graphics device
  on.exit({
    graphics::par(userpar)
  })


  # Set plot geometry.
  if( cbarplot == FALSE & rbarplot == FALSE & legend == FALSE ){
    # One panel.
  }
  if( cbarplot == FALSE & rbarplot == TRUE & legend == FALSE ){
     graphics::layout(matrix(1:2, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6))
  }
  if( cbarplot == FALSE & rbarplot == FALSE & legend == TRUE ){
     graphics::layout(matrix(1:2, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.2))
  }
  if( cbarplot == FALSE & rbarplot == TRUE & legend == TRUE ){
     graphics::layout(matrix(1:3, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6, 0.2))
  }
  if( cbarplot == TRUE & rbarplot == FALSE & legend == FALSE ){
    graphics::layout(matrix(1:2, nrow=nrows, ncol=ncols, byrow = TRUE),
           heights=c(1, 4))    
  }
  if( cbarplot == TRUE & rbarplot == TRUE & legend == FALSE ){
    graphics::layout(matrix(1:4, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6), heights=c(1, 4))
  }
  if( cbarplot == TRUE & rbarplot == FALSE & legend == TRUE ){
    graphics::layout(matrix(1:4, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.2), heights=c(1, 4))
  }
  if( cbarplot == TRUE & rbarplot == TRUE & legend == TRUE ){
    graphics::layout(matrix(1:6, nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6, 0.2), heights=c(1, 4))
  }
  
  # Global parameters.
  graphics::par(mar=c(0,0,0,0))
  graphics::par(oma=c(1,1,1,1))
  
  if( cbarplot == TRUE ){
    graphics::barplot(colSums(x, na.rm=na.rm),
            space=0, border=NA, axes=FALSE,
            names.arg="",
            col=c("#808080", "#c0c0c0"), xaxs="i")
    if(clabels == TRUE & scale == 'none'){
      graphics::text(c(1:ncol(x))-0.5, 0.0, colnames(x), adj=c(0.0,0.5), srt=90)
    } else if (clabels == TRUE & scale != 'none'){
      graphics::text(c(1:ncol(x))-0.5, min(colSums(x, na.rm=na.rm), na.rm=na.rm), colnames(x), adj=c(0.0,0.5), srt=90)
    }
    if( rbarplot == TRUE ){
      plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
    }
    if( legend == TRUE ){
      plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
    }
  }
  
  # Plot image matrix.
  graphics::image(t(x), col = col.ramp,
                  axes=FALSE, frame.plot=TRUE)
  
  # Row barplot.
  if( rbarplot == TRUE ){
    graphics::barplot(rowSums(x, na.rm=na.rm),
            space=0, border=NA,
            horiz=TRUE, axes=FALSE, names.arg="",
            col=c("#808080", "#c0c0c0"),
            yaxs="i")
    if(rlabels == TRUE & scale == 'none'){
      graphics::text(0, c(1:nrow(x))-0.5, rownames(x), adj=c(0.0, 0.5), srt=0)
    } else if(rlabels == TRUE & scale != 'none'){
      graphics::text(min(rowSums(x, na.rm=na.rm), na.rm=na.rm), c(1:nrow(x))-0.5, rownames(x), adj=c(0.0,0.5), srt=0)
    }
  }
  
  # Legend.
  if( legend == TRUE ){
    mp <- graphics::barplot(rep(1, times=length(col.ramp)), space=0, border=NA, horiz = TRUE,
                            col = col.ramp, axes=FALSE)
#    graphics::text(0.5, 5, "Low", col="#FFFFFF")
#    graphics::text(0.5, 95, "High", col="#FFFFFF")
    if ( mp[nrow(mp),1] - mp[1,1] >= 1 ){
      graphics::text(0.5, mp[1,1], "Low", col="#FFFFFF", adj=c(0.5,0))
      graphics::text(0.5, mp[nrow(mp),1], "High", col="#FFFFFF", adj=c(0.5,1))
    }
  }
  invisible(NULL)
}





##### ##### ##### ##### #####
# EOF
