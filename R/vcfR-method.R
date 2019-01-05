

##### Method show #####
#' 
#' @rdname vcfR-method
# ' @aliases show.vcfR,show,vcfR-method
#' @aliases show,vcfR-method
#' @title show
#' 
#' @description 
#' Display a summary of a vcfR object.
#' 
#' @param object a vcfR object
#' 
#' @details 
#' The method \strong{show} is used to display an object.
#' Because vcf data are relatively large, this has been abbreviated.
#' Here we display the first four lines of the meta section, and truncate them to no more than 80 characters.
#' The first eight columns and six rows of the fix section are also displayed. 
#' 
setMethod(
  f="show",
  signature = "vcfR",
  definition=function(object){
    
    if( ncol(object@gt) > 1 ){
      nsamp <- ncol(object@gt) - 1
    } else {
      nsamp <- 0
    }
    nchrom <- length( unique( getCHROM( object ) ) )
    nvar <- nrow(object@fix)
    nna <- sum( is.na(object@gt[,-1]) )
    pna <- nna / c( nsamp * nvar )

    cat("***** Object of Class vcfR *****\n")
    cat( paste( nsamp, "samples\n") )
    cat( paste( nchrom, "CHROMs\n") )
    cat( paste( format(nvar, big.mark=","), "variants\n") )
    cat( "Object size: ")
    print(object.size(object), units="MB")
    cat( paste( format(pna * 100, digits = 4), "percent missing data\n") )
    cat("*****        *****         *****\n")
#    message("*****        --*--         *****")
  }
)


#### Method head ####
#' 
#' @name head
#' @rdname vcfR-method
#' @title head
#' @aliases head,vcfR-method
#' @docType methods
#' 
#' @param x object of class vcfR
#' @param n number of rows to print
#' @param maxchar maximum number of characters to print per line
#' @param ... arguments to be passed to other methods
#' 
#' @description \strong{head} returns the first parts of an object of class vcfR.
#' 
#' @details 
#' The method \strong{head} is similar to show, but is more flexible.
#' The number of rows displayed is parameterized by  the variable n.
#' And the maximum number of characters to print per line (row) is also parameterized.
#' In contract to show, head includes a summary of the gt portion of the vcfR object.
#' 
#' 
setMethod(
  f="head",
  signature="vcfR",
  definition=function (x, n=6, maxchar=80){
    print("***** Object of class 'vcfR' *****")
    print("***** Meta section *****")
    
    if(length(x@meta) > n){
      for( i in 1:n ){
        if( nchar(x@meta[i]) <= maxchar ){
          print(x@meta[i])
        } else {
          print( paste( substr(x@meta[i], 1, maxchar-12 ), "[Truncated]" ) )
        }
      }
      print(paste("First", n, "rows."))
    } else {
      print(x@meta)
    }
    
    print("", quote=FALSE)

    print("***** Fixed section *****")
    if(nrow(x@fix) >= n){
      print(x@fix[1:n,1:7])
    } else {
      print(x@fix[,1:7])
    }
    print("", quote=FALSE)

    print("***** Genotype section *****")
    if(nrow(x@gt) >= n){
      if(ncol(x@gt)<6){
        print(x@gt[1:n,])
      } else {
        print(x@gt[1:n,1:6])
        print("First 6 columns only.")
      }
    } else {
      if(ncol(x@gt)<6){
        print(x@gt)
      } else {
        print(x@gt[,1:6])
      }
    }
    print("", quote=FALSE)
    
    print("Unique GT formats:")
    if( nrow(x@gt) == 0 ){
      print("No gt slot present")
      print("", quote=FALSE)      
    } else {
      print(unique(as.character(x@gt[,1])))
      print("", quote=FALSE)      
    }
  }
)


#### Method [] ####
#' 
#' @rdname vcfR-method
#' @title Brackets
#' @description The brackets ('[]') subset objects of class vcfR
#' @details
#' The \strong{square brackets ([])} are used to subset objects of class vcfR.
#' Rows are subset by providing a vector i to specify which rows to use.
#' The columns in the fix slot will not be subset by j.
#' The parameter j is a vector used to subset the columns of the gt slot.
#' Note that it is essential to include the first column here (FORMAT) or downsream processes will encounter trouble.
#' 
#' The \strong{samples} parameter allows another way to select samples.
#' Because the first column of the gt section is the FORMAT column you typically need to include that column and sample numbers therefore begin at two.
#' Use of the samples parameter allows you to select columns by a vector of numerics, logicals or characters.
#' When numerics are used the samples can be selected starting at one.
#' The function will then add one to this vector and include one to select the desired samples and the FORMAT column.
#' When a vector of characters is used it should contain the desired sample names.
#' The function will add the FORMAT column if it is not the first element.
#' When a vector of logicals is used a TRUE will be added to the vector to ensure the FORMAT column is selected.
#' Note that specification of samples will override specification of j.
#' 
#' 
# @export
# @aliases []
#'
#' @aliases [,vcfR-method
#'
#' @param i vector of rows (variants) to include
#' @param j vector of columns (samples) to include
#' @param samples vector (numeric, character or logical) specifying samples, see details
#' @param drop delete the dimensions of an array which only has one level
#'
setMethod(
  f= "[",
  signature(x = "vcfR"),
#  signature(x = "vcfR", i = "ANY", j = "ANY"),
#  signature(x = "vcfR", i = "ANY", j = "ANY", samples = "ANY"),
  definition=function(x, i, j, samples = NULL, ..., drop){
#  definition=function(x, i, j, ..., drop){
    if( !is.null(samples) ){
      if( inherits(samples, what  = c("numeric", "integer") ) ){
        samples <- samples + 1
        j <- c(1, samples)
      } else if( inherits(samples, what  = "character") ){
        if( samples[1] != "FORMAT" ){ 
          j <- c("FORMAT", samples)
        } else {
          j <- samples
        }
      } else if( inherits(samples, what  = "logical") ){
        j <- c(TRUE, samples)
      } else {
        stop(paste("samples specified, expecting a numeric, character or logical but received", class(samples)))
      }
    }

    if(nrow(x@gt) == nrow(x@fix)){
      x@gt <- x@gt[ i, j, drop = FALSE ]
    } else if (nrow(x@gt) == 0){
      # Do nothing.
    } else {
      msg <- paste("The fix slot has", nrow(x@fix), "rows while the gt slot has", nrow(x@gt), "rows, this should never happen.")
      stop(msg)
    }
    x@fix <- x@fix[ i, , drop = FALSE ]
    
    if(nrow(x@gt) > 0){
      if(colnames(x@gt)[1] != 'FORMAT'){
        warning("You have chosen to omit the FORMAT column, this is typically undesireable.")        
      }
    }

    return(x)
  }
)


setGeneric("plot")
#### Method plot ####
#'
#' @rdname vcfR-method
#' @aliases plot,vcfR-method
#' 
#' @title plot.vcfR
#' @description The \strong{plot} method visualizes objects of class vcfR
# @export
# ' @aliases plot.vcfR
# ' @aliases vcfR,vcfR-method
#' 
#' @param y not used
#' 
#' @details 
#' The \strong{plot} method generates a histogram from data found in the 'QUAL' column from the 'fix' slot.
#' 
setMethod(
  f="plot",
  signature= "vcfR",
  definition=function(x, y, ...){
    x <- as.numeric(x@fix[,'QUAL'])
    graphics::hist(x, col=5, main='Histogram of qualities', xlab='QUAL')
    graphics::rug(x)
  }
)


##### ##### ##### ##### #####
#
# rbind
#
##### ##### ##### ##### #####


setMethod("rbind",
  signature( "vcfR" ),
  function (..., deparse.level = 0) 
  {
    ## store arguments
    dots <- list(...)
    
    ## extract arguments which are vcfR objects
    myList <- dots[sapply(dots, inherits, "vcfR")]
    if(!all(sapply(myList, class)=="vcfR")) stop("some objects are not vcfR objects")
    
    ## keep the rest in 'dots'
    dots <- dots[!sapply(dots, inherits, "vcfR")]
    
    # Initialize
    x <- myList[[1]]
    
    browser()
    # Implement
    x@fix <- do.call( rbind, lapply( myList, function(x){ x@fix } ) )
    x@gt  <- do.call( rbind, lapply( myList, function(x){ x@gt  } ) )
    
    return(x)
  }
)



#' @rdname vcfR-method
#' @aliases rbind2.vcfR
#' 
setMethod("rbind2",
  signature(x = "vcfR", y = "missing"),
  function (x, y, ...) 
  {
#    message("y is missing.")
    return(x)
  }
)

#' @rdname vcfR-method
#' @aliases rbind2.vcfR
#' 
setMethod("rbind2",
  signature(x = "vcfR", y = "ANY"),
  function (x, y, ...) 
  {
#    message("y is ANY.")
    return(x)
  }
)

#setGeneric("rbind2")
#' @rdname vcfR-method
#' @aliases rbind2.vcfR
#' @export
#' 
setMethod("rbind2",
  signature( x="vcfR", y="vcfR" ),
  function (x, y, ...)
  {
#    message("rbind2.vcfR")
#    browser()
    x@fix <- rbind( x@fix, y@fix )
    x@gt  <- rbind( x@gt,  y@gt  )
    return(x)
  }
)



##### ##### ##### ##### #####
#
# nrow
#
##### ##### ##### ##### #####


#' @rdname vcfR-method
#' @aliases dim.vcfR
#' @export
#' 
setMethod("dim",
  signature(x = "vcfR"),
  function (x) 
  {
    x <- c( nrow(x@fix), ncol(x@fix), ncol(x@gt) )
    names(x) <- c( 'variants', 'fix_cols', 'gt_cols')
    return(x)
  }
)


#' @rdname vcfR-method
#' @aliases nrow.vcfR
#' @export
#' 
setMethod("nrow",
  signature(x = "vcfR"),
  function (x) 
  {
    rows <- nrow(x@fix)
    return(rows)
  }
)


##### ##### ##### ##### #####
# EOF.