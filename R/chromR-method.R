#'
#' @rdname chromR-method
#' @title chromR-method
#'
#' @aliases chromR,chromR-method
#'
#' @description Methods that act on objects of class chromR
#'
#'
#' @param x an object of class chromR
#' @param y not currently used
#' @param object an object of class chromR
#' @param value a character containing a name
#' @param n integer indicating the number of elements to be printed from an object
#' @param ... Arguments to be passed to methods
#'
#'
#' @details
#' Methods that act on objects of class chromR.
#'
#' @importFrom utils object.size
#'
#'



##### Generic methods. #####

setMethod( f="show",
  signature = "chromR",
  definition=function(object){
            #1234567890123456789012345678901234567890
    cat( "*****   Class chromR, method Show   *****\n" )
#    cat( "\n" )
    cat( paste("Name:", object@names, "\n") )
#    cat( "\n" )
    cat( paste("Chromosome length:", format(object@len, big.mark=","), "bp\n") )
    cat( "  Chromosome labels: ")
    if( length( labels(object@seq) ) > 0 ){
      cat( paste( labels(object@seq), sep = ",") )
    } else {
      cat( "None" )
    }
    cat( "\n" )
    cat( paste("Annotation (@ann) count:", format(nrow(object@ann), big.mark=","), "\n") )
    cat( "  Annotation chromosome names: " )
    if( length( unique( object@ann[,1] ) ) > 0 ){
      cat( paste( unique( object@ann[,1] ) ), sep = "," )
    } else {
      cat( "None" )
    }
    cat( "\n" )
    cat( paste("Variant (@vcf) count:", format(nrow(object@vcf), big.mark=","), "\n") )
    cat( "  Variant (@vcf) chromosome names: " )
    if( length( unique(getCHROM(object@vcf)) ) > 0 ){
      cat( paste( unique(getCHROM(object@vcf)), sep = "," ) )
    } else {
      cat( "None" )
    }
    cat( "\n" )
#    cat( "\n" )
    cat( "Object size: ")
    print( utils::object.size(object), units="MB" )
    cat( "Use head(object) for more details.\n" )
#    cat( "\n" )
    cat( "*****      End Show (chromR)        *****\n" )
  }
)



#' @rdname chromR-method
# ' @aliases plot
#' @aliases plot,chromR-method
#' @export
#'
setMethod( f="plot",
  signature= "chromR",
  definition=function (x,y,...){
    DP <- x@var.info$DP[x@var.info$mask]
    MQ <- x@var.info$MQ[x@var.info$mask]
    QUAL <- as.numeric(x@vcf@fix[x@var.info$mask, 'QUAL'])

    if( nrow(x@win.info ) > 0 ){
      SNPS <- x@win.info$variants/x@win.info$length
    } else {
      SNPS <- NULL
    }

    graphics::par(mfrow=c(2,2))
    if( length(stats::na.omit(DP)) > 0 ){
      graphics::hist(DP, col=3, main="Read depth (DP)", xlab="")
      graphics::rug(DP)
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No depths found")
    }
    if( length(stats::na.omit(MQ)) > 0 ){
      graphics::hist(MQ, col=4, main="Mapping quality (MQ)", xlab="")
      graphics::rug(MQ)
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No mapping qualities found")
    }
    if( length(stats::na.omit(QUAL)) > 0 ){
      graphics::hist(QUAL, col=5, main="Quality (QUAL)", xlab="")
      graphics::rug(QUAL)
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No qualities found")
    }
    if( length(SNPS) > 0 ){
      graphics::hist( SNPS, col=6, main="Variant count (per window)", xlab="")
      graphics::rug( SNPS )
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No SNP densities found")
    }
    graphics::par(mfrow=c(1,1))
    return(invisible(NULL))
  }
)



##### ##### ##### ##### #####


#' @rdname chromR-method
#' @export
#'
setMethod( f="print",
  signature="chromR",
  definition=function (x,y,...){
    cat("***** Object of class 'chromR' *****\n")
    cat(paste("Name: ", x@names, "\n"))
    cat(paste("Length: ", format(x@len, big.mark=","), "\n"))
    cat("\nVCF fixed data:\n")
    cat("Last column (info) omitted.\n")
    cat("\nVCF variable data:\n")
    cat(paste("Columns: ", ncol(x@vcf@gt), "\n"))
    cat(paste("Rows: ", nrow(x@vcf@gt), "\n"))
    cat("(First column is format.)\n")
    cat("\nAnnotation data:\n")
    if(length(x@ann)>0){
      print(head(x@ann[,1:8], n=4))
      cat("Last column (attributes) omitted.\n")
    } else {
      cat("Empty slot.\n")
    }
    cat("***** End print (chromR) ***** \n")
  }
)


#' @rdname chromR-method
#' @export
#'
setMethod( f="head",
  signature = "chromR",
  definition=function(x, n = 6){
            #1234567890123456789012345678901234567890
    cat("*****   Class chromR, method head   *****")
    cat("\n")
    cat(paste("Name: ", x@names))
    cat("\n")
    cat( paste("Length: ", format(x@len, big.mark=",")) )
    cat("\n")
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****     Sample names (chromR)     *****")
    cat("\n")
    if(ncol(x@vcf@gt) <= 2 * n){
      print(colnames(x@vcf@gt)[-1])
    } else {
      print(head(colnames(x@vcf@gt)[-1]))
      print("...")
      print(utils::tail(colnames(x@vcf@gt)[-1]))
    }
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****    VCF fixed data (chromR)    *****")
    cat("\n")
    if(nrow(x@vcf@gt) <= 2 * n){
      print(x@vcf@fix[,1:7])
    } else {
      print(head(x@vcf@fix[,1:7]))
      print("...")
      print(utils::tail(x@vcf@fix[,1:7]))
    }
    cat("\n")
    cat("INFO column has been suppressed, first INFO record:")
    cat("\n")
    print(unlist(strsplit(as.character(x@vcf@fix[1, 'INFO']), split=";")))
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****  VCF genotype data (chromR)  *****")
    cat("\n")
    if(ncol(x@vcf@gt)>=6){
              #1234567890123456789012345678901234567890
      cat("*****     First 6 columns      *********")
      cat("\n")
      if(nrow(x@vcf@gt) <= 2 * n){
        print(x@vcf@gt[,1:6])
      } else {
        print(head(x@vcf@gt[,1:6]))
#        print("...")
#        print(tail(x@vcf@gt[,1:6]))
      }
    } else {
      print(x@vcf@gt[1:6,])
    }
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****      Var info (chromR)       *****")
    cat("\n")
    if(ncol(x@var.info)>=6){
              #1234567890123456789012345678901234567890
      cat("*****       First 6 columns        *****")
      cat("\n")
      print(x@var.info[1:n,1:6])
    } else {
      print(x@var.info[1:n,])
    }
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****      VCF mask (chromR)        *****")
    cat("\n")
    cat( paste("Percent unmasked:", round(100*(sum(x@var.info$mask)/length(x@var.info$mask)), digits=2 ) ) )
    cat("\n")
    cat("\n")
            #1234567890123456789012345678901234567890
    cat("*****      End head (chromR)        *****")
    cat("\n")

  }
)


#' @rdname chromR-method
#'
setMethod(f="names<-",
          signature( x = "chromR", value = "character" ),
          function(x, value){
            if( length(value) >=1 ){
              x@names <- value[1]
            } else {
              x@names <- character()
            }
            return(x)
          }
)



#' @rdname chromR-method
#' @export
#'
setMethod( f="length",
  signature = "chromR",
  definition=function(x){
    return(x@len)
  }
)


# EOF.