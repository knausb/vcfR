#'
#' @rdname vcfR-methods
#' 
# @export
#' @aliases vcfR-methods
#' @title vcfR methods
#' 
#' @description Methods to show, subset or plot data from objects of class vcfR.
#' 
#' 
#' @param object object of class vcfR
#' @param x object of class vcfR
#' @param n number of rows to print
#' @param maxchar maximum number of characters to print per line
#' 
#' @details 
#' The method \strong{show} is used to display an object.
#' Because vcf data are relatively large, this has been abbreviated.
#' Here we display the first four lines of the meta section, and truncate them to no more than 80 characters.
#' The first eight columns and six rows of the fix section are also displayed. 
#' 
#' The method \strong{head} is similar to show, but is more flexible.
#' The number of rows displayed is parameterized by  the variable n.
#' And the maximum number of characters to print per line (row) is also parameterized.
#' In contract to show, head includes a summary of the gt portion of the vcfR object.
#' 
#' The \strong{square brackets ([])} are used to subset objects of class vcfR.
#' Rows are subset by providing a vector i to specify which rows to use.
#' The columns in the fix slot will not be subset by j.
#' The parameter j is a vector used to subset the columns of the gt slot.
#' Note that it is essential to include the first column here (FORMAT) or downsream processes will encounter trouble.
#' 
#' The \strong{plot} method generates a histogram from data found in the 'QUAL' column from the 'fix' slot.
#' 
#' 


##### Method show #####
#' 
setMethod(
  f="show",
  signature = "vcfR",
  definition=function(object){
    
    nsamp <- ncol(object@gt) - 1
    nvar <- nrow(object@gt)
    nna <- sum( is.na(object@gt[,-1]) )
    pna <- nna / c( nsamp * nvar )

    message("***** Object of Class vcfR *****")
    message( paste( nsamp, "samples") )
    message( paste( format(nvar, big.mark=","), "variants") )
    message( paste( format(pna * 100, digits = 4), "percent missing data") )
    message("*****        *****         *****")
#    message("*****        --*--         *****")
  }
)


#### Method head ####
#' 
#' @rdname vcfR-methods
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
    print(unique(as.character(x@gt[,1])))
    print("", quote=FALSE)
  }
)



#### Method [] ####
#' @rdname vcfR-methods
#' 
# @export
# @aliases []
#'
#' @param i vector of rows (variants) to include
#' @param j vector of columns (samples) to include
#' @param drop delete the dimensions of an array which only has one level
#'
setMethod(
  f= "[",
  signature="vcfR",
  definition=function(x, i, j, drop){
    x@fix <- x@fix[ i, , drop = FALSE ]
    x@gt <- x@gt[ i, j, drop = FALSE ]
    return(x)
  }
)


setGeneric("plot")
#### Method plot ####
#'
#' @rdname vcfR-methods
#' @export
#' @aliases plot.vcfR
#' 
#' @param y not used
#' @param ... Arguments to be passed to methods
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


#setGeneric("rbind", signature="...")
#setGeneric("rbind")
# '
# ' @rdname vcfR-methods
# ' @aliases rbind.vcfR, rbind
# ' 
# ' @param deparse.level integer controlling the construction of labels. see ?rbind.
#  @param ... objects of class vcfR passed to rbind
# ' 
# ' @export
# '
#setMethod(
#  f = "rbind",
#  signature = "vcfR", 
#  definition = function( ..., deparse.level=1 )
#  {
    ## store arguments
#    dots <- list(...)

    ## extract arguments which are vcfR objects
#    myList <- dots[sapply(dots, inherits, "vcfR")]
#    if(!all(sapply(myList, class)=="vcfR")) stop("some objects are not vcfR objects")
    
    ## keep the rest in 'dots'
#    dots <- dots[!sapply(dots, inherits, "vcfR")]
    
    # Initialize
#    x <- myList[[1]]
    
    # Implement
#    x@fix <- do.call( rbind, lapply( myList, function(x){ x@fix } ) )
#    x@gt  <- do.call( rbind, lapply( myList, function(x){ x@gt  } ) )

#    return(x)
#  }
#)



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


#' @rdname vcfR-methods
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

#' @rdname vcfR-methods
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
#' @rdname vcfR-methods
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


#' @rdname vcfR-methods
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


#' @rdname vcfR-methods
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