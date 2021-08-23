#' 
#' @title Convert a vcfR object to hapmap
#' @description The function converts a vcfR object to hapmap
#'
#' @param vcf a vcfR object.
# ' @param out_file name of output file.
# ' @param method should 'N' or 'H' format data be generated?
#'
#' @details
#' This function converts a vcfR object to a hapmat format.
#' 
#' @return a data.frame that can be used as an input for GAPIT.
#' 
#' @author Brian J. Knaus
#' 
# ' @seealso \href{http://popgen.sc.fsu.edu/Migrate/Migrate-n.html}{Migrate-N} website.
# ' @examples
#' 
#' @export
vcfR2hapmap <- function(vcf) {
#  print("vcfR2hapmap works!")
  
  vcf <- vcf[!is.indel(vcf), ]
  vcf <- vcf[is.biallelic(vcf), ]
  gt <- extract.gt(vcf, return.alleles = TRUE)
  gt <- sub("/|\\|", "", gt, fixed = FALSE)

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
  hapMap <- rbind(hapMap[1, ], hapMap)
  
  return(hapMap)
}
  

