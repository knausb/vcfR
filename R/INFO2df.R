


#' @title Reformat INFO data as a data.frame
#' @name INFO2df
#' @rdname INFO2df
#' 
#' @description
#' Reformat INFO data as a data.frame and handle class when possible.
#'  
#' @param x an object of class vcfR or chromR.
#' 
#' @details 
#' The INFO column of VCF data contains descriptors for each variant.
#' Because this column may contain many comma delimited descriptors it may be difficult to interpret.
#' The function INFO2df converts the data into a data.frame.
#' The function metaINFO2df extracts the information in the meta section that describes the INFO descriptors.
#' This function is called by INFO2df to help it handle the class of the data.
#' 
#' @return 
#' A data.frame
#' 
#' 
#' @examples
#' data(vcfR_test)
#' metaINFO2df(vcfR_test)
#' getINFO(vcfR_test)
#' INFO2df(vcfR_test)
#' 
#' 
#' @export
#' 
INFO2df <- function(x){
  
  if( inherits(x, "chromR") ){
    x <- x@vcfR
  }  
  
  INFO.Type <- metaINFO2df(x)

  # Initialize a data.frame for the INFO data
  INFOdf <- data.frame( matrix( nrow=nrow(x@fix),  ncol=nrow(INFO.Type) ) )
  names(INFOdf) <- INFO.Type[,'ID']
  
  for( i in 1:nrow(INFO.Type) ){
    tmp <- extract.info(x, element = INFO.Type[,'ID'][i])
    if( INFO.Type[,'Type'][i] == "Integer" & INFO.Type[,'Number'][i] == "1" ){
      tmp <- as.integer( tmp )
    }
    if( INFO.Type[,'Type'][i] == "Float" & INFO.Type[,'Number'][i] == "1" ){
      tmp <- as.numeric( tmp )
    }
    INFOdf[,INFO.Type[,'ID'][i]] <- tmp
  }
  return(INFOdf)
}


#' @rdname INFO2df
#'
#' @export
#'
metaINFO2df <- function(x){

  if( inherits(x, "chromR") ){
    x <- x@vcfR
  }
  
  # Isolate INFO from meta
  INFO <- x@meta[grep("##INFO=", x@meta)]
  
  # Clean things up a bit.
  INFO <- sub("##INFO=<", "", INFO)
  INFO <- sub(">$", "", INFO)
  INFO <- sub('\"', "", INFO)
  INFO <- sub('\"$', "", INFO)
  
  INFO <- strsplit(INFO, split = ",")
  
  ID          <- unlist( lapply(INFO, function(x){ grep("^ID=", x, value=TRUE) }) )
  Number      <- unlist( lapply(INFO, function(x){ grep("^Number=", x, value=TRUE) }) )
  Type        <- unlist( lapply(INFO, function(x){ grep("^Type=", x, value=TRUE) }) )
  Description <- unlist( lapply(INFO, function(x){ grep("^Description=", x, value=TRUE) }) )
  Source      <- unlist( lapply(INFO, function(x){ grep("^Source=", x, value=TRUE) }) )
  Version     <- unlist( lapply(INFO, function(x){ grep("^Version=", x, value=TRUE) }) )
  
  ID          <- unlist(lapply(strsplit(ID, split = "=" ), function(x){x[2]}))
  Number      <- unlist(lapply(strsplit(Number, split = "=" ), function(x){x[2]}))
  Type        <- unlist(lapply(strsplit(Type, split = "=" ), function(x){x[2]}))
  Description <- unlist(lapply(strsplit(Description, split = "=" ), function(x){x[2]}))
  Source      <- unlist(lapply(strsplit(Source, split = "=" ), function(x){x[2]}))
  Version     <- unlist(lapply(strsplit(Version, split = "=" ), function(x){x[2]}))

#  INFO.Type <- cbind(ID, Number, Type, Description, Source, Version)
  INFO.Type <- data.frame( ID = ID, Number = Number, Type = Type, 
                           Description = Description, stringsAsFactors = FALSE)
  if( !is.null(Source) ) { INFO.Type$Source = Source }
  if( !is.null(Version) ){ INFO.Type$Version = Version }
  
  return( INFO.Type )  
}

