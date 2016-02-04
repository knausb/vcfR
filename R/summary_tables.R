#'
#' 
#' @title Write summary tables from chromR objects
#' @rdname summary_tables
#' @export
#'
#' @description
#' Write summary tables from chromR objects.
#' 
#' @param file A filename for the output file
#' @param x An object of class chromR
#' @param mask logical vector indicating rows to use
#' @param APPEND logical indicating whether to append to existing file (omitting the header) or write a new file
#'
#' @details
#' The function \strong{write.var.info} takes the variant information table from a chromR object and writes it as a comma delimited file. 
#' 
#' The function \strong{write.win.info} takes the window information table from a chromR object and writes it as a comma delimited file.
#' 
#'
#' @seealso
#' \code{\link{write.vcf}}
#'
# CRAN:
# \href{http://cran.r-project.org/web/packages/pegas/index.html}{pegas}::read.vcf,
# \href{http://cran.r-project.org/web/packages/PopGenome/index.html}{PopGenome}::readVCF,
# \href{http://cran.r-project.org/web/packages/data.table/index.html}{data.table}::fread
# 
# Bioconductor:
# \href{http://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation}::readVcf
#


#' @rdname Summary.tables
#' @aliases write.var.info
#' 
#' @export
#' 
write.var.info <- function(x, file = "", mask = FALSE, APPEND = FALSE){
  if(class(x) == "vcfR"){
    stop("Unexpected class! Detected class vcfR. This class does not contain variant summaries.")
  }
  if(class(x) != "chromR"){
    stop("Unexpected class! Expecting an object of class chromR.")
  }
  
  if(mask == FALSE){
    utils::write.table(x@var.info, file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
  } else if(mask == TRUE){
    utils::write.table(x@var.info[x@var.info$mask,], file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
  }
}



#' @rdname summary_tables
#' @aliases write.win.info
#' 
#' @export
#' 
write.win.info <- function(x, file = "", APPEND = FALSE){
  if(class(x) == "vcfR"){
    stop("Unexpected class! Detected class vcfR. This class does not contain window summaries.")
  }
  if(class(x) != "chromR"){
    stop("Unexpected class! Expecting an object of class chromR.")
  }
  
  utils::write.table(x@win.info, file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
}

