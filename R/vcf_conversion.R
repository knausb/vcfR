#' @title Convert vcf data to other formats
#' @name Format conversion
#' @rdname vcf_conversion
#' @description
#' Convert vcfR object to other formats
#'  
#' @param x A matrix of genotypes


#' @rdname vcf_conversion
#' @export
gt_to_df <- function(x) {
  fmat.col <- grep("FORMAT", names(x))
  gt.pos <- grep("^GT$", unlist(strsplit(as.character(x[1, fmat.col]), ":")))
  
  gt2snpbin <- function(x) {
    gt.a <- unlist(lapply(strsplit(as.character(x), ":"), function(x) {
      x[gt.pos]
    }))
    gt.a <- unlist(lapply(strsplit(gt.a, "/"), function(x) {
      paste(as.numeric(x[1]), as.numeric(x[2]), sep = "/")
    }))
    gt.a <- gsub("0", "2", gt.a)
  }
  gt.df <- data.frame(matrix(nrow = ncol(x) - fmat.col, ncol = nrow(x)))
  for (i in c(fmat.col + 1):ncol(x)) {
    gt.df[i - fmat.col, ] <- gt2snpbin(x[, i])
  }
  rownames(gt.df) <- colnames(x)[c(fmat.col + 1):ncol(x)]
  gt.df
}

