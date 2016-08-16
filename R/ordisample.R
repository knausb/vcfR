#' 
#' 
#' @title Ordinate a sample's data
#' @name ordisample
#' @rdname ordisample
#' 
#' @description
#' Ordinate information from a sample's GT region and INFO column.
#'
#' @param x an object of class vcfR or chromR.
#' @param sample a sample number where the first sample (column) is 2
#' @param distance metric to be used for ordination, options are in \code{\link[vegan]{vegdist}}
#' @param plot logical specifying whether to plot the ordination
#' @param verbose logical specifying whether to produce verbose output
#' @param ... parameters to be passed to child processes
#' 
#' 
#' @details 
#' The INFO column of VCF data contains descriptors for each variant.
#' Each sample typically includes several descriptors of each variant in the GT region as well.
#' This can present an overwhelming amount of information.
#' Ordination is used in this function to reduce this complexity.
#' 
#' The ordination procedure can be rather time consuming depending on how much data is used.
#' I good recommendation  is to always start with a small subset of your full dataset and slowly scale up.
#' There are several steps in this function that attempt to eliminate variants or characters that have missing values in them.
#' This that while starting with a small number is good, you will need to have a large enough number so that a substantial amount of the data make it to the ordination step.
#' In the example I use 100 variants which appears to be a reasonable compromise.
#' 
#' 
#' Several of the steps of this ordination make distributional assumptions.
#' That is, they assume the data to be normally distributed.
#' There is no real reason to assume this assumption to be valid with VCF data.
#' It has been my experience that this assumption is frequently violated with VCF data.
#' It is therefore suggested to use this funciton as an exploratory tool that may help inform other decisions.
#' These analyst may be able to address these issues through data transformation or other topics beyond the scope of this function.
#' Interpretation of the results of this function need to take into account that assumptions may hav been violated.
#' 
#' 
#' @return 
#' A list consisting of two objects.
#' The first object is an object of class 'metaMDS' created by the function vegan::metaMDS.
#' The second is an object of class 'envfit' created by the function vegan::envfit.
#' 
#' 
#' @seealso
#' \code{\link[vegan]{metaMDS}}
#' \code{\link[vegan]{vegdist}}
#' 
#' @examples
#' \dontrun{
# data(vcfR_test)
#' data(vcfR_example)
#' ordisample(vcf[1:100,], sample = 2)
#' 
#' myOrd <- ordisample(vcf[1:100,], sample = 2, plot = FALSE)
#' names(myOrd)
#' plot(myOrd$metaMDS, type = "n")
#' points(myOrd$metaMDS, display = "sites", pch=20, col="#8B4513")
#' text(myOrd$metaMDS, display = "spec", col="blue")
#' plot(myOrd$envfit, col = "#008000", add = TRUE)
#' }
#' 
#' 
#' @import vegan
#' @export
#' 
ordisample <- function(x, sample, distance = "bray", plot = TRUE, verbose = TRUE, ...){
#  require(vegan, quietly = verbose)
  if( class(sample) == "character" ){
    sample <- grep( sample, colnames(x@gt), fixed = TRUE )
  }
  if( length(sample) !=1 ){
    stop( "Invalid specification of 'sample.' Please use either an integer or a character." )
  }
  
  x <- x[,c(1,sample)]
  
  # INFO data
  myINFO <- INFO2df(x)
  myMETA <- metaINFO2df(x)

  for(i in 1:ncol(myINFO)){
    tmp <- myINFO[,i]
    if( class(myINFO[,i]) == "character" ){
      tmp <- strsplit(tmp, split = ",")
      tmp <- lapply(tmp, function(x){x[1]})
      tmp <- unlist(tmp)

      if( myMETA$Type[i] == "Integer"){
        tmp <- as.integer(tmp)
      }
      if( myMETA$Type[i] == "Float"){
        tmp <- as.numeric(tmp)
      }
    }
    myINFO[,i] <- tmp
  }
  
  # Get FORMAT fields
  myFORMAT <- metaINFO2df(x, field = "FORMAT")
  myFORMAT <- myFORMAT[grep("^GT$", myFORMAT$ID, invert = TRUE),]
  
  myGT <- data.frame( matrix( nrow=nrow(x), ncol=nrow(myFORMAT) ) )
  names(myGT) <- myFORMAT$ID
  
  for(i in 1:ncol(myGT)){
    tmp <- extract.gt( x, element = colnames(myGT)[i] )
      tmp <- strsplit(tmp, split = ",")
      tmp <- lapply(tmp, function(x){x[1]})
      tmp <- unlist(tmp)
    if( myFORMAT$Type[i] == "Integer"){
      tmp <- as.integer(tmp)
    }
    if( myFORMAT$Type[i] == "Float"){
      tmp <- as.numeric(tmp)
    }
    myGT[,i] <- tmp
  }
  
  # Manage NAs.
  # GT missingness.
  badVars <- apply(myGT, MARGIN=1, function(x){ sum( is.na(x) ) > 0 })
  if( verbose == TRUE & sum(badVars) > 0 ){
    message(paste( sum(badVars), "variants containing missing values removed." ))
  }
  myGT   <- myGT[!badVars,]
  myINFO <- myINFO[!badVars,]
  
  
  badChars <- apply(myINFO, MARGIN=2, function(x){ sum( is.na(x) ) == length(x) })
  if( verbose == TRUE & sum(badChars) > 0 ){
    message(paste("INFO character: ", names(myINFO)[badChars], " omitted due to missingness.\n" ))
  }
  myINFO <- myINFO[,!badChars]
  
  badChars <- apply(myINFO, MARGIN=2, function(x){ length( unique(x) ) == 1 })
  if( verbose == TRUE & sum(badChars) > 0 ){
    message(paste("INFO character: ", names(myINFO)[badChars], " omitted as monomorphic.\n" ))
  }
  myINFO <- myINFO[,!badChars]

  # INFO missingness.
  badVars <- apply(myINFO, MARGIN=1, function(x){ sum( is.na(x) ) > 0 })
  if( verbose == TRUE & sum(badVars) > 0 ){
    message(paste( sum(badVars), "variants containing missing values removed." ))
  }
  myGT   <- myGT[!badVars,]
  myINFO <- myINFO[!badVars,]
  
  # Ordination
  mds1 <- vegan::metaMDS(myGT, distance = distance, k = 2)

  # Covariates
  ord.fit <- vegan::envfit(mds1, env=myINFO, perm=999, na.rm = TRUE)

  # Plot
  if( plot == TRUE ){
    graphics::plot(mds1, type = "n")
    graphics::points(mds1, display = "sites", cex = 0.8, pch=20, col=grDevices::rgb(139,69,19, alpha=50, maxColorValue = 255))
    vegan::ordiellipse(mds1, groups = factor( rep(1, times=nrow(myGT)) ), kind = "sd", conf = 0.68, col = "#808080" )  
    graphics::text(mds1, display = "spec", col="blue", ...)
    graphics::title( main = colnames(x@gt)[2] )
  
    graphics::plot(ord.fit, choices = c(1,2), at = c(0,0),
         axis = FALSE,
    #     p.max = 0.05,
         col = "#008000", add = TRUE, ... )
  }
  
  invisible( list( metaMDS = mds1, envfit = ord.fit ) )
}

