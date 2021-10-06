#' 
#' @title Convert a vcfR object to hapmap
#' 
#' @description Converts a vcfR object to hapmap
#'
#' @param vcf a vcfR object.
# ' @param out_file name of output file.
# ' @param method should 'N' or 'H' format data be generated?
#'
#' @details
#' Converts a vcfR object to a hapmap format.
#' 
#' @return a data.frame that can be used as an input for GAPIT.
#' 
#' @author Brian J. Knaus
#' 
# ' @seealso \href{http://popgen.sc.fsu.edu/Migrate/Migrate-n.html}{Migrate-N} website.
#' @examples
#' data(vcfR_test)
#' myHapMap <- vcfR2hapmap(vcfR_test)
#' class(myHapMap)
#' \dontrun{
#' # Example of how to create a (GAPIT compliant) HapMap file.
#' write.table(myHapMap, 
#'             file = "myHapMap.hmp.txt",
#'             sep = "\t", 
#'             row.names = FALSE,
#'             col.names = FALSE)
#' }
#' 
#' @export
vcfR2hapmap <- function(vcf) {
#  print("vcfR2hapmap works!")
  
  vcf <- vcf[!is.indel(vcf), ]
  vcf <- vcf[is.biallelic(vcf), ]
  gt <- extract.gt(vcf, return.alleles = TRUE)
  gt <- sub("/|\\|", "", gt, fixed = FALSE)
  gt[ is.na(gt) ] <- "NN"
  gt[ gt == "." ] <- "NN"

  hapMap <- matrix(data = NA, nrow = nrow(gt), ncol = ncol(gt) + 11)
  hapMap <- as.data.frame(hapMap)
  colnames(hapMap) <- c(
    c("rs", "alleles", "chrom", "pos", "strand", "assembly", "center",
      "protLSID", "assayLSID", "panel", "QCcode"),
    colnames(gt)
    )
  hapMap[,1] <- rownames(gt)
  hapMap[,3] <- getCHROM(vcf)
  hapMap[,4] <- getPOS(vcf)
  hapMap[, 12:ncol(hapMap)] <- gt
  class(hapMap) <- c("hapMap", class(hapMap))
  
  # GAPIT compatibility
  hapMap <- rbind(colnames(hapMap), hapMap)
  
  return(hapMap)
}
  

