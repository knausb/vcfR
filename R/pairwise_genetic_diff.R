#' @title Pairwise genetic differentiation across populations
#'
#' \code{pairwise_genetic_diff} Calculate measures of genetic differentiation across all population pairs. 
#'
#' @param vcf a vcfR object
#' @param pops factor indicating populations
#' @param method the method to measure differentiation
#' @return a data frame containing the pairwise population differentiation indices of interest across all pairs of populations in the population factor.
#' @examples
#' data(vcfR_example)
#' myPops <- as.factor(rep(c('a','b'), each = 9))
#' myDiff <- pairwise_genetic_diff(vcf, myPops, method = "nei")
#' colMeans(myDiff[,c(4:ncol(myDiff))], na.rm = TRUE)
#' 
#' @seealso  \code{\link{pairwise_genetic_diff}} in  \code{\link{vcfR}}

pairwise_genetic_diff <- function (vcf, pops, method="nei"){
  var_info <- as.data.frame(vcf@fix[, 1:2, drop = FALSE])
  if (is.null(var_info$mask)) {
    var_info$mask <- TRUE
  }
  combination.df <- combn(as.character(unique(myPops)), 2, simplify = F)
  test <- lapply(combination.df, function (x) {
    vcf.temp <- vcf
    pop.tem <- as.factor(as.character(myPops[myPops %in% x]))
    samples.temp <- colnames(vcf.temp@gt)[-1][myPops %in% x]
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
  })
  pop.diff <- as.data.frame(do.call(cbind,test))
  pop.diff <- cbind(var_info, pop.diff)
  return(pop.diff)
}

