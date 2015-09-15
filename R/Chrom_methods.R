
#' @title Chrom methods
#' 
#' @param x an object of class Chrom
#' @param object an object of class Chrom
#' @param y some sort of object???
#' @param ... Arguments to be passed to methods
#' 
#' @rdname Chrom_methods
#

##### ##### Generic methods. #####

setMethod(
  f="show",
  signature = "Chrom",
#  definition=function(x){
  definition=function(object){
            #1234567890123456789012345678901234567890
    print("*****   Class Chrom, method Show   *****")
    print(paste("Name: ", object@name))
    print(paste("Length: ", object@len))
    print("Use head(object) for more details.")
#    print(paste("Name: ", x@name, "\n"))
#    print(paste("Length: ", x@len, "\n"))
    #    print("Use head(x) for more details.\n")    
            #1234567890123456789012345678901234567890
    print("*****      End Show (Chrom)        *****")
  }
)

setMethod(
  f="print",
  signature="Chrom",
  definition=function (x,y,...){
    print("***** Object of class 'Chrom' *****\n")
    print(paste("Name: ", x@name, "\n"))
    print(paste("Length: ", x@len, "\n"))
    print("\nVCF fixed data:\n")
    print("Last column (info) omitted.\n")
    print("\nVCF variable data:\n")
    print(paste("Columns: ", ncol(x@vcf@gt), "\n"))
    print(paste("Rows: ", nrow(x@vcf@gt), "\n"))
    print("(First column is format.)\n")
    print("\nAnnotation data:\n")
    if(length(x@ann)>0){
      print(head(x@ann[,1:8], n=4))
      print("Last column (attributes) omitted.\n")
    } else {
      print("Empty slot.\n")
    }
    print("***** End print (Chrom) ***** \n")
  }
)



##### Basic methods (definitions) #####

#' @rdname Chrom_methods
#' @export
#' @aliases names.chrom
#' 
setMethod(
  f="names",
  signature = "Chrom",
  definition=function(x){
    #    cat("**** Class Chrom, method names **** \n")
    #    cat("Sequence name: ", as.character(names(x@seq)), "\n")
    #    cat("First annotation name: ")
    #    print(as.character(x@ann[1,1]))
    #    cat("First variant name: ")
    #    print(as.character(x@vcf.fix[1,1]))
    #    cat("\n")
    #    cat("Sample names: \n")
    temp <- names(x@vcf@gt)[-1]
    temp
  }
)


#' @rdname Chrom_methods
#' @export
#' @aliases head.chrom
#' 
setMethod(
  f="head",
  signature = "Chrom",
  definition=function(x){
            #1234567890123456789012345678901234567890
    print("*****   Class Chrom, method head   *****")
    print(paste("Name: ", x@name))
    print(paste("Length: ", x@len))
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****     Sample names (Chrom)     *****")
    print(colnames(x@vcf@gt)[-1])
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****    Vcf fixed data (Chrom)    *****")
    print(x@vcf@fix[1:6,1:7])
    print('', quote=FALSE)
    print("INFO column has been suppressed, first INFO record:")
    print(unlist(strsplit(as.character(x@vcf@fix[1, 'INFO']), split=";")))
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****   Vcf genotype data (Chrom)  *****")
    if(ncol(x@vcf@gt)>=6){
              #1234567890123456789012345678901234567890
      print("*****     First 6 columns      *********")
      print(x@vcf@gt[1:6,1:6])
    } else {
      print(x@vcf@gt[1:6,])
    }
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****      Var info (Chrom)        *****")
    if(ncol(x@var.info)>=6){
              #1234567890123456789012345678901234567890
      print("*****       First 6 columns        *****")
      print(x@var.info[1:6,1:6])
    } else {
      print(x@var.info[1:6,])
    }
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****      Vcf mask (Chrom)        *****")
    print(paste("Percent unmasked:", 100*(sum(x@var.info$mask)/length(x@var.info$mask))))
    print('', quote=FALSE)
            #1234567890123456789012345678901234567890
    print("*****      End head (Chrom)        *****")
    print('', quote=FALSE)
  }
)



#' @rdname Chrom_methods
#' @export
#' @aliases plot.chrom
#' 
setMethod(
  f= "plot",
  signature= "Chrom",
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

    
    par(mfrow=c(2,2))
    if( length(DP) > 0 ){
      hist(DP, col=3, main="Read depth (DP)", xlab="")
      rug(DP)
    } else {
      plot(1:2,1:2, type='n')
      title(main="No depths found")
    }
    if( length(MQ) > 0 ){
      hist(MQ, col=4, main="Mapping quality (MQ)", xlab="")
      rug(MQ)
    } else {
      plot(1:2,1:2, type='n')
      title(main="No mapping qualities found")
    }
    if( length(QUAL) > 0 ){
      hist(QUAL, col=5, main="Quality (QUAL)", xlab="")
      rug(QUAL)
    } else {
      plot(1:2,1:2, type='n')
      title(main="No qualities found")
    }
    if( length(SNPS) > 0 ){
      hist( SNPS, col=6, main="Variant count (per window)", xlab="")
      rug( SNPS )
    } else {
      plot(1:2,1:2, type='n')
      title(main="No SNP densities found")
    }
    par(mfrow=c(1,1))
  }
)


##### Accessors.  #####

#### Getter for "names" ####
setGeneric("getName",function(object){standardGeneric ("getName")})

setMethod("getName","Chrom",
          function(object){
            return(object@name)
          }
)

#### Setter for name. ####

setGeneric("setName<-",function(object,value){standardGeneric("setName<-")})

setReplaceMethod(
  f="setName",
  signature="Chrom",
  definition=function(object,value){
    object@name <-value
    return (object)
  }
)





# Setter for seq. ####

#setGeneric("seq2chrom<-",function(object,value){standardGeneric("seq2chrom<-")})

#setReplaceMethod(
#  f="seq2chrom",
#  signature="Chrom",
#  definition=function(object,value){
#    # A DNAbin will store in a list when the fasta contains
#    # multiple sequences, but as a matrix when the fasta
#    # only contains one sequence.
#    if(!is.list(class(as.character(value)))){
#      object@seq <- as.list(value)
#    } else {
#      object@seq <-value      
#    }
#    object@len <-length(object@seq[[1]])
#    return (object)
#  }
#)

##### ##### ##### ##### #####









##### ##### win.info functions #####

#' @rdname Chrom_methods
#' @export
#' @aliases windowize
#'
# @description
# Creates windows
#'
#' @param win.size window size, in base pairs
#' @param max.win maximum window size
#'
#' @details
#' Reads in a vcf file and stores it in a Chrom class.
#' 
#'
windowize <- function(x, win.size=1000, max.win=10000){
  #  acgt.w <- x@acgt.w
  acgt.w <- x@seq.info$nuc.win
  windows <- matrix(NA, ncol=2, nrow=max.win)
  i <- 1
  for(j in 1:nrow(acgt.w)){
    while(acgt.w[j,2]-acgt.w[j,1] > win.size){
      windows[i,1] <- acgt.w[j,1]
      windows[i,2] <- acgt.w[j,1] + win.size - 1
      acgt.w[j,1] <- acgt.w[j,1] + win.size + 0
      i <- i+1
      if(i > max.win){
        print(paste("max i equals", max.win))
        print(paste("i equals", i))
        print(paste("j equals", j))
        print("chrom.r error: max.win is too small.\n")
        break
      }
    }
    windows[i,1] <- acgt.w[j,1]
    windows[i,2] <- acgt.w[j,2]
    i <- i+1
  }
  x@windows <- windows[1:i-1,]
  return(x)
}

gc.win <- function(x){
  win <- matrix(ncol=7,
                nrow=nrow(x@windows),
                dimnames=list(c(),
                              c('index','start','stop','gc','at','gcf','atf'))
  )
  win[,1] <- 1:nrow(win)
  win[,2] <- x@windows[,1]
  win[,3] <- x@windows[,2]
  chrom <- as.character(x@seq)[[1]]
  #
  count.nucs <- function(x){
    chrom <- chrom[x[2]:x[3]]
    win[x[1],4] <<- length(grep("[GgCc]", chrom, perl=TRUE))
    win[x[1],5] <<- length(grep("[AaTt]", chrom, perl=TRUE))
  }
  #
  apply(win, MARGIN=1, count.nucs)
  win[,6] <- win[,4]/(win[,3]-win[,2])
  win[,7] <- win[,5]/(win[,3]-win[,2])
  x@nuccomp.w <- as.data.frame(win)
  return(x)
}

snp.win <- function(x){
  snp <- matrix(ncol=5,
                nrow=nrow(x@windows),
                dimnames=list(c(),
                              c('index','start','stop','count','density'))
  )
  snp[,1] <- 1:nrow(snp)
  snp[,2] <- x@windows[,1]
  snp[,3] <- x@windows[,2]
  vcf <- x@var.info$POS[x@var.info$mask]
  #
  count.snps <- function(x){
    vcf2 <- vcf[vcf >= x[2] & vcf <= x[3]]
    snp[x[1],4] <<- length(vcf2)
  }
  apply(snp, MARGIN=1, count.snps)
  snp[,5] <- snp[,4]/(snp[,3]-snp[,2]+1)
  x@snpden.w <- as.data.frame(snp)
  return(x)
}



##### ##### vcf functions #####

vcf.fix2gt.m <- function(x){
  snames <- names(x@vcf@gt)[-1]
  pos <- paste(x@vcf@fix[,1], x@vcf@fix[,2], sep="_")
  #  pos <- x@vcf.fix[,2]
  x1 <- as.matrix(x@vcf@gt)
  nsamp <- ncol(x1) - 1
  #
  x1 <- cbind(unlist(lapply(strsplit(x1[,1], ":"), function(x){grep("GT", x)})),x1)
  #
  get.gt <- function(x){
    cell <- as.numeric(x[1])
    x <- lapply(strsplit(x[3:length(x)], ":"), function(x){x[cell]})
    unlist(x)
  }
  x1 <- apply(x1, MARGIN=1, get.gt)
  x1[x1=="0/0"] <- 0
  x1[x1=="0/1"] <- 1
  x1[x1=="1/0"] <- 1
  x1[x1=="1/1"] <- 2
  x1 <- as.numeric(x1)
  x1 <- matrix(data=x1, ncol=nsamp, byrow=TRUE,
               dimnames=list(pos, snames))
  x@gt.m <- x1
  return(x)
}





