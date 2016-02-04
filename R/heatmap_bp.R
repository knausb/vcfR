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
#' @seealso \code{\link[stats]{heatmap}}, \code{\link[graphics]{image}}, heatmap2 in \href{https://cran.r-project.org/package=gplots}{gplots}.
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
#'  
# data(vcfR_example)
# pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
# pinf_mt <- masker(pinf_mt)
# pinf_gq <- extract.gt(pinf_mt, element="GQ", as.numeric=TRUE)
# heatmap.bp(pinf_gq)
# heatmap.bp(pinf_gq, scale="col")
# heatmap.bp(pinf_gq, col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100))
# heatmap.bp(pinf_gq, col.ramp = colorRampPalette(c("#D55E00", "#F0E442", "#009E73"))(100))
#' 
heatmap.bp <- function(x, cbarplot = TRUE, rbarplot = TRUE,
                       legend = TRUE, clabels = TRUE, rlabels = TRUE,
                       na.rm = TRUE, scale = c("row", "column", "none"),
#                       col.ramp = colorRampPalette(c("red", "yellow", "#008000"))(100),
#                       col.ramp = viridis::viridis(n = 100, alpha=1),
                       col.ramp = viridisLite::viridis(n = 100, alpha=1),
                       ...){
#  require(viridis)
  stopifnot(class(x) == 'matrix')
  scale <- if(missing(scale))
    "none"
  else match.arg(scale)
  #
  nrows <- 1
  ncols <- 1
  if(cbarplot == TRUE){ nrows <- nrows + 1 }
  if(rbarplot == TRUE){ ncols <- ncols + 1 }
  if(rbarplot == TRUE){ ncols <- ncols + 1 }
  #
  #
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
  #
  if(is.null(colnames(x))){colnames(x) <- 1:ncol(x)}
  if(is.null(rownames(x))){rownames(x) <- 1:nrow(x)}
  #
  if(ncols == 1 & nrows == 1){
    cat("1 and 1\n")
    graphics::image(t(x),
          col = col.ramp,          
          axes=FALSE, frame.plot=TRUE)
  }
  if(nrows == 1 & ncols == 2){}
  if(nrows == 1 & ncols == 3){}
  if(nrows == 2 & ncols == 1){}
  if(nrows == 2 & ncols == 2){}
  if(nrows == 2 & ncols == 3){
    graphics::layout(matrix(1:6,nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6, 0.2), heights=c(1, 4))
    graphics::par(mar=c(0,0,0,0))
    graphics::par(oma=c(1,1,1,1))
    graphics::barplot(colSums(x, na.rm=na.rm),
            space=0, border=NA, axes=FALSE,
            names.arg="",
            col=c("#808080", "#c0c0c0"), xaxs="i")
    if(clabels == TRUE & scale == 'none'){
      graphics::text(c(1:ncol(x))-0.5, 0.0, colnames(x), adj=c(0.0,0.5), srt=90)
    } else if (clabels == TRUE & scale != 'none'){
      graphics::text(c(1:ncol(x))-0.5, min(colSums(x, na.rm=na.rm), na.rm=na.rm), colnames(x), adj=c(0.0,0.5), srt=90)
    }
    #
    plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
    plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
    graphics::image(t(x),
          col = col.ramp,
          axes=FALSE, frame.plot=TRUE)
    #
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
    #
    graphics::barplot(rep(1, times=length(col.ramp)), space=0, border=NA, horiz = TRUE,
            col = col.ramp,
            axes=FALSE)
    graphics::text(0.5, 5, "Low", col="#FFFFFF")
    graphics::text(0.5, 95, "High", col="#FFFFFF")
  }
  graphics::par(mfrow=c(1,1))
  graphics::par(mar=c(5,4,4,2))
  #
}





##### ##### ##### ##### #####
# EOF