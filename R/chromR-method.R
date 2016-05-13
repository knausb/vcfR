#'
#' @rdname chromR-method
#' @title chromR-method
#'
#' @aliases chromR-method
#'
#' @description Methods that act on objects of class chromR
#'
#'
#' @param x an object of class chromR
#' @param y some sort of object???
#' @param object an object of class chromR
#' @param value a character containing a name
#' @param ... Arguments to be passed to methods
#'
#'
#' @details
#' Methods that act on objects of class chromR.
#'
#'



##### Generic methods. #####

setMethod( f="show",
  signature = "chromR",
  definition=function(object){
            #1234567890123456789012345678901234567890
    cat( "*****   Class chromR, method Show   *****\n" )
    cat( paste("Name: ", object@name, "\n") )
    cat( paste("Length: ", object@len, "\n") )
    cat( paste("Object size:", print(object.size(object), units="MB"), "MB.\n") )
    cat( "Use head(object) for more details.\n" )
    cat( "*****      End Show (chromR)        *****\n" )
  }
)


#' @rdname chromR-method
#'
setMethod( f="plot",
  signature= "chromR",
  definition=function (x,y,...){
    DP <- x@var.info$DP[x@var.info$mask]
    MQ <- x@var.info$MQ[x@var.info$mask]
    QUAL <- as.numeric(x@vcf@fix[x@var.info$mask, 'QUAL'])

#    if( length(DP) < 0 )

    if( nrow(x@win.info ) > 0 ){
#    if( na.omit(x@win.info$variants) > 0 ){
      SNPS <- x@win.info$variants/x@win.info$length
    } else {
      SNPS <- NULL
    }


    graphics::par(mfrow=c(2,2))
    if( length(DP) > 0 ){
      graphics::hist(DP, col=3, main="Read depth (DP)", xlab="")
      graphics::rug(DP)
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No depths found")
    }
    if( length(MQ) > 0 ){
      graphics::hist(MQ, col=4, main="Mapping quality (MQ)", xlab="")
      graphics::rug(MQ)
    } else {
      plot(1:2,1:2, type='n', xlab="", ylab="")
      graphics::title(main="No mapping qualities found")
    }
    if( length(QUAL) > 0 ){
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
  }
)



##### ##### ##### ##### #####


#' @rdname chromR-method
#'
setMethod( f="print",
  signature="chromR",
  definition=function (x,y,...){
    cat("***** Object of class 'chromR' *****\n")
    cat(paste("Name: ", x@name, "\n"))
    cat(paste("Length: ", x@len, "\n"))
    cat("\nVCF fixed data:\n")
    cat("Last column (info) omitted.\n")
    cat("\nVCF variable data:\n")
    cat(paste("Columns: ", ncol(x@vcf@gt), "\n"))
    cat(paste("Rows: ", nrow(x@vcf@gt), "\n"))
    cat("(First column is format.)\n")
    cat("\nAnnotation data:\n")
    if(length(x@ann)>0){
      cat(head(x@ann[,1:8], n=4))
      cat("Last column (attributes) omitted.\n")
    } else {
      cat("Empty slot.\n")
    }
    cat("***** End print (chromR) ***** \n")
  }
)


setMethod( f="head",
  signature = "chromR",
  definition=function(x){
            #1234567890123456789012345678901234567890
    cat("*****   Class chromR, method head   *****")
    cat(paste("Name: ", x@name))
    cat(paste("Length: ", x@len))
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****     Sample names (chromR)     *****")
    cat(colnames(x@vcf@gt)[-1])
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****    Vcf fixed data (chromR)    *****")
    cat(x@vcf@fix[1:6,1:7])
    cat('', quote=FALSE)
    cat("INFO column has been suppressed, first INFO record:")
    cat(unlist(strsplit(as.character(x@vcf@fix[1, 'INFO']), split=";")))
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****   Vcf genotype data (chromR)  *****")
    if(ncol(x@vcf@gt)>=6){
              #1234567890123456789012345678901234567890
      cat("*****     First 6 columns      *********")
      cat(x@vcf@gt[1:6,1:6])
    } else {
      cat(x@vcf@gt[1:6,])
    }
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****      Var info (chromR)        *****")
    if(ncol(x@var.info)>=6){
              #1234567890123456789012345678901234567890
      cat("*****       First 6 columns        *****")
      cat(x@var.info[1:6,1:6])
    } else {
      cat(x@var.info[1:6,])
    }
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****      Vcf mask (chromR)        *****")
    cat( paste("Percent unmasked:", round(100*(sum(x@var.info$mask)/length(x@var.info$mask)), digits=2 ) ) )
    cat('', quote=FALSE)
            #1234567890123456789012345678901234567890
    cat("*****      End head (chromR)        *****")
    cat('', quote=FALSE)
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


# EOF.