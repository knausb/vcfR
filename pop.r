# Pop.
##### ##### ##### ##### #####
# Class definition.

setClass(
  Class="Pop",
  representation=representation(
    stats = "data.frame"
  ),
  prototype=prototype(
    stats = data.frame(matrix(ncol=8, nrow=0, 
                              dimnames=list(c(),
c('R_num','A_num','Ho','Ne','theta_pi','theta_w','theta_h','tajimas_d'))),
                       stringsAsFactors=FALSE)
  )
)

##### ##### ##### ##### #####
# Generic methods.

setMethod(
  f="show",
  signature = "Pop",
  definition=function(object){
    cat("*** Class pop, method Show *** \n")
    cat("\nStats:\n")
    print(head(object@stats))
    cat("\nSFS:\n")
    print(head(object@sfs[,1:min(ncol(object@sfs), 8)]))
    cat("Use print(object) for more details.\n")
    cat("******* End Show 'Pop' ******* \n")
  }
)

setMethod(
  f="print",
  signature="Pop",
  definition=function (x,y,...){
    cat("***** Object of class 'Pop' *****\n")
    cat("\nStats:\n")
    print(head(object@stats))
    cat("\nSFS:\n")
    print(head(object@sfs[,1:min(ncol(object@sfs), 8)]))
    cat("***** End print 'Pop' ***** \n")
  }
)

setMethod(
  f= "plot",
  signature= "Pop",
  definition=function (x,y,...){
    cat("***** Object of class 'Pop' *****\n")
    cat("***** Plot not yet implemented *****\n")
  }
)

##### ##### ##### ##### #####
# Functions.

create.pop <- function(x1, x2){
  if(class(x2) != "Chrom"){
    cat("Error: expecting object of class 'Chrom'")
    break
  }
  #
  vcf2gt <- function(x, cell = 1) {
    get.gt <- function(y) {unlist(lapply(strsplit(as.character(y), split = ":"), function(z) {z[cell]}))}
    gt <- apply(x[,2:ncol(x)], 2, get.gt)
    if(class(gt) != "matrix"){
      gt <- matrix(gt, ncol=length(2:ncol(x)))
    }
    gt
  }
#  gt <- 
  gt
}




##### ##### ##### ##### #####
# EOF.
