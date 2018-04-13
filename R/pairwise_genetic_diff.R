#' @title Pairwise genetic differentiation across populations
#' 
#' @aliases pairwise_genetic_diff
#' 
#' @description
#' \code{pairwise_genetic_diff} Calculate measures of genetic differentiation across all population pairs. 
#'
#' @param vcf a vcfR object
#' @param pops factor indicating populations
#' @param method the method to measure differentiation
#' 
#' @author Javier F. Tabima
#' 
#' @return a data frame containing the pairwise population differentiation indices of interest across all pairs of populations in the population factor.
#' 
#' @examples
#' data(vcfR_example)
#' pops <- as.factor(rep(c('a','b'), each = 9))
#' myDiff <- pairwise_genetic_diff(vcf, pops, method = "nei")
#' colMeans(myDiff[,c(4:ncol(myDiff))], na.rm = TRUE)
#' pops <- as.factor(rep(c('a','b','c'), each = 6))
#' myDiff <- pairwise_genetic_diff(vcf, pops, method = "nei")
#' colMeans(myDiff[,c(4:ncol(myDiff))], na.rm = TRUE)
#' 
#' @seealso  \code{\link{genetic_diff}} in  \code{\link{vcfR}}
#' 
#' @export
pairwise_genetic_diff <- function (vcf, pops, method="nei"){
  var_info <- as.data.frame(vcf@fix[, 1:2, drop = FALSE])
  if (is.null(var_info$mask)) {
    var_info$mask <- TRUE
  }
  # Create a list of pairwise comparisons.
  combination.df <- utils::combn(x = as.character(unique(pops)), m = 2, simplify = FALSE)
  
  # Function to make pairwise comparisons
  pwDiff <- function (x) {
    # x contains the names of two populations to be compared.
    # pops is a factor of population designations for each sample.
    # vcf is the vcfR object.
    # method is the method to be used in the comparison.
    vcf.temp <- vcf
    pop.tem <- as.factor(as.character(pops[pops %in% x]))
    samples.temp <- colnames(vcf.temp@gt)[-1][pops %in% x]
    vcf.temp@gt <- vcf.temp@gt[, c(TRUE, colnames(vcf.temp@gt)[-1] %in% samples.temp)]
    temp.gendif <- genetic_diff(vcf.temp, pop.tem, method = method)
    if (method == "nei") {
      temp.genind <- temp.gendif[,colnames(temp.gendif) %in% c("Gst","Gprimest")]
      colnames(temp.genind) <- paste0(colnames(temp.genind),"_", paste0(levels(pop.tem), collapse = "_"))
    } else if (method == "jost") {
      temp.genind <- temp.gendif[,colnames(temp.gendif) %in% c("Dest_Chao","Db")]
      colnames(temp.genind) <- paste0(colnames(temp.genind),"_", paste0(levels(pop.tem), collapse = "_"))
    }
    return(temp.genind)
  }
  
  test <- lapply(combination.df, pwDiff)
  pop.diff <- as.data.frame(do.call(cbind,test))
  pop.diff <- cbind(var_info, pop.diff)
  return(pop.diff)
}

