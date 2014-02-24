#### Misc functions ####

#' @title Heatmap with barplots
#'
#' @name heatmap-barplots
#' @rdname heatmap.bp
#' @aliases heatmap.bp
#' @export
#' 
#' @description
#' Heatmap of a numeric matrix with barplots summarizing columns and rows.
#'
#' 
#' @param x a numeric matrix
#' @param rbarplot logical indicating whether the rows should be summarized with a barplot
#' @param cbarplot logical indicating whether the columns should be summarized with a barplot
#' @param legend logical indicating whether a legend should be plotted
#' @param ...  additional arguments to be passed on
#' 
#' @details The function heatmap.bp creates a heatmap from a numeric matrix with optional barplots to summarize the rows and columns.
#' 
#' @seealso heatmap, image, heatmap2 in library(gplots)
#' 
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' pinf_mt <- create.chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
#' pinf_gq <- extract.gt(pinf_mt, element="GQ", as.matrix=TRUE)
#' heatmap.bp(pinf_gq)
#' 
heatmap.bp <- function(x, cbarplot = TRUE, rbarplot = TRUE, legend = TRUE, ...){
  stopifnot(class(x) == 'matrix')
  #
  nrows <- 1
  ncols <- 1
  if(cbarplot == TRUE){ nrows <- nrows + 1 }
  if(rbarplot == TRUE){ ncols <- ncols + 1 }
  if(rbarplot == TRUE){ ncols <- ncols + 1 }
  #
  cat("nrows = ", nrows, "\n")
  cat("ncols = ", ncols, "\n")
  #
  if(ncols == 1 & nrows == 1){
    cat("1 and 1\n")
    image(t(x), col = colorRampPalette(c("red", "yellow", "#008000"))(100),
          axes=FALSE, frame.plot=TRUE)
  }
  if(nrows == 1 & ncols == 2){}
  if(nrows == 1 & ncols == 3){}
  if(nrows == 2 & ncols == 1){}
  if(nrows == 2 & ncols == 2){}
  if(nrows == 2 & ncols == 3){
    cat("2 and 3\n")
    layout(matrix(1:6,nrow=nrows, ncol=ncols, byrow = TRUE),
           widths=c(4, 0.6, 0.2), heights=c(1, 4))
    par(mar=c(0,0,0,0))
    par(oma=c(1,1,1,1))
#    barplot(colSums(x), space=0, border=NA, axes=FALSE, col=c("#808080", "#c0c0c0"))
#    plot(colSums(x), type="h", axes=FALSE, col=c("#808080", "#c0c0c0"), lwd=4)
    par(mar=c(2,0,0,0))
    plot(0:ncol(x), seq(0, max(colSums(x)), length=c(ncol(x)+1)), type="n",
         axes=FALSE,
         plt=c(0,0,0,0))
    rect(xleft=0:c(ncol(x)-1), ybottom=0, xright=1:ncol(x), ytop=colSums(x), 
         col=c("#808080", "#c0c0c0"),
         border=NA, xlim=c(0,ncol(x)+1))
    axis(side=1)
    #
    plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
    plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
#    plot(1, 1, type="n")
#    plot(1:10,1:10)
    image(t(x), col = colorRampPalette(c("red", "yellow", "#008000"))(100),
          axes=FALSE, frame.plot=TRUE)
    #
    barplot(rowSums(x), space=0, border=NA, horiz=TRUE, axes=FALSE, col=c("#808080", "#c0c0c0"))
    #
    barplot(rep(1, times=100), space=0, border=NA, horiz = TRUE,
            col=colorRampPalette(c("red", "yellow", "#008000"))(100),
            axes=FALSE)
    text(0.5, 5, "Low", col="#FFFFFF")
    text(0.5, 95, "High", col="#FFFFFF")
  }
  par(mfrow=c(1,1))
  #

  #
#  image(t(x), col = colorRampPalette(c("red", "yellow", "#008000"))(100),
#        axes=FALSE, frame.plot=TRUE)
}




#### EOF ####