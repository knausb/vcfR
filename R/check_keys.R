

#### check_keys ####
#' @rdname check_keys
#' @aliases check_keys
#' 
#' @title Check that INFO and FORMAT keys are unique
#' 
#' @param x an oblect of class vcfR
#' 
#' @description 
#' The INFO and FORMAT columns contain information in key-value pairs.
#' If for some reason a key is not unique it will create issues in retrieving this information.
#' This function checks the keys defined in the meta section to make sure they are unique.
#' Note that it does not actually check the INFO and FORMAT columns, just their definitions in the meta section.
#' This is because each variant can have different information in their INFO and META cells.
#' Checking these on large files will tehrefore come with a performance cost.
#' 
#' @seealso queryMETA()
#' 
#' @examples 
#' data(vcfR_test)
#' check_keys(vcfR_test)
#' queryMETA(vcfR_test)
#' queryMETA(vcfR_test, element = 'DP')
#' # Note that DP occurs as unique in INFO and FORMAT but they may be different.
#' 
#' 
#' @export
check_keys <- function(x) {
  if(class(x) != 'vcfR'){
    stop( paste('Expecting a vcfR object, instead received:', class(x)) )
  }
  
  # First check INFO.
  myKeys <- grep('INFO', x@meta, value = TRUE)
  myKeys <- sub('##INFO=<ID=','',myKeys)
  myKeys <- unlist(lapply(strsplit(myKeys, ','), function(x){x[1]}))
  myKeys <- table(myKeys)
  myKeys <- myKeys[myKeys > 1]
  if( length(myKeys) > 0){
    warning(paste("The following INFO key occurred more than once:", names(myKeys), '\n'))
  }
  
  # Check FORMAT.
  myKeys <- grep('FORMAT', x@meta, value = TRUE)
  myKeys <- sub('##FORMAT=<ID=','',myKeys)
  myKeys <- unlist(lapply(strsplit(myKeys, ','), function(x){x[1]}))
  myKeys <- table(myKeys)
  myKeys <- myKeys[myKeys > 1]
  if( length(myKeys) > 0){
    warning(paste("The following FORMAT key occurred more than once:", names(myKeys), '\n'))
  }

}



