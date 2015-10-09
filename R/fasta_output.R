#'
#'
#' @title Create fasta format output
#' @rdname fasta_output
#' @aliases write.fasta
#' 
#' @description Generate fasta format output
#' 
#' @param x object of class chromR
#' @param file name for output file
#' @param gt_split character which delimits alleles in genotype
#' @param rowlength number of characters each row should not exceed
#' @param tolower convert all characters to lowercase (T/F)
#' @param verbose should verbose output be generated (T/F)
#' @param APPEND should data be appended to an existing file (T/F)
#' 
#' 
#' @details 
#' The function \strong{write_fasta} takes an object of class chromR and writes it to a fasta.gz (gzipped text) format file.
#' The sequence in the seq slot of the chromR object is used to fill in the invariant sites.
#' The parameter 'tolower', when set to TRUE, converts all the characters in teh sequence to lower case.
#' This is important because some software, such as ape::DNAbin, requires sequences to be in lower case.
#' 
#' 
#' @export
#' 
write.fasta <- function(x, file = "", gt_split = "|", rowlength=80, tolower=TRUE, verbose=TRUE, APPEND = FALSE){
  if(class(x) != "chromR"){
    stop("Expected object of class chromR")
  }
  if(APPEND == FALSE){
    if(file.exists(file)){
      file.remove(file)
    }
  }
  haps <- extract_haps(x, gt_split = gt_split)
  if(tolower == TRUE){
    haps <- apply(haps, MARGIN=2, tolower)
  }
  
  for(i in 1:ncol(haps)){
    seq <- as.character(x@seq)[1,]
#    seq[x@vcf.fix$POS] <- haps[,i]
    seq[x@var.info$POS] <- haps[,i]
    invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, colnames(haps)[i], file, rowlength, as.integer(verbose)))
  }
  
  #invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, seqname, filename, rowlength, verbose))  
}
