#' @title Convert a vcfR object to MigrateN input file
#' @description The function converts a vcfR object to a text format that can be used as an infile for MigrateN.
#'
#' @param vcf a vcfR object.
#' @param pop factor indicating population membership for each sample.
#' @param in_pop vector of population names indicating which population to include in migrate output file.
#' @param out_file name of output file.
#' @param method should 'N' or 'H' format data be generated?
#' 
#' @return a text file that can be used as an input for MigrateN software (SNP format).
#'
#' @details
#' This function converts a vcfR object to a text file which can be used as input for MigrateN.
#' The function will remove loci with missing data, indels, and loci that are not bialleleic (loci with more than two alleles). 
#' Thus, only SNP data analysed where the length of each locus (inmutational steps) is 1 (as opposed to microsatellites or indels).
#' 
#' The output file should contain Unix line endings ("\\n").
#' Note that opening the output file in a Windows text editor (just to validate number of markers, individuals or populations) might change the end of line character (eol) to a Windows line ending ("\\r\\n"). 
#' This may produce an error running migrate-n.
#' Because these are typically non-printing characters, this may be a difficult problem to troubleshoot.
#' The easiest way to circumvent the problem is to transfer the output file to Unix machine and view it there.
#' If you do introduce Windows line endings you can convert them back to Unix with a program such as `dos2unix` or `fromdos` to change the line endings.
#' 
#' @author Shankar Shakya and Brian J. Knaus
#' 
#' @seealso \href{http://popgen.sc.fsu.edu/Migrate/Migrate-n.html}{Migrate-N} website.
#'
#' @examples
#' \dontrun{
# ' pkg <- "pinfsc50"
# ' my_vcf <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
# ' my_vcf <- read.vcfR( my_vcf, verbose = FALSE )
# ' 
#' data(vcfR_example)
#' my_pop <- as.factor(paste("pop_", rep(c("A", "B", "C"), each = 6), sep = ""))
#' vcfR2migrate(vcf = vcf , pop = my_pop , in_pop = c("pop_A","pop_C"),
#'              out_file = "my2pop.txt", method = 'H')
#' }
#'
#'
#' @export
vcfR2migrate <- function(vcf, pop, in_pop, out_file = "MigrateN_infile.txt", method = c('N','H') ) {

  method <- match.arg(method, c('N','H'), several.ok = FALSE)
  
  # Validate the input.
  if( class(vcf) != "vcfR"){
    stop(paste("Expecting an object of class vcfR, received a", class(vcf), "instead")) 
  }
  if( class(pop) != "factor"){
    stop(paste("Expecting population vector, received a", class(pop), "instead")) 
  }

  # Remove indels and non-biallelic loci
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf),]

  # Remove loci containing missing genotypes.
  gt <- extract.gt(vcf, convertNA = T)
  vcf <- vcf[!rowSums((is.na(gt))),]

#  FORMAT <- vcf@gt[1:nrow(gt),1]
#  vcf@gt <- cbind(FORMAT, gt)

  # Subset VCF data to populations.
#  vcf_list <- lapply(levels(my_pop), function(x){ vcf[,c(TRUE, x == my_pop)] })
#  names(vcf_list) <- levels(pop)
  vcf_list <- lapply(in_pop, function(x){ vcf[,c(TRUE, x == pop)] })
  names(vcf_list) <- in_pop
  
#  for (i in (1:length(vcf_list))) {
#    temp_pop <- names(vcf_list[i])
#    temp_vcf <- vcf
#    FORMAT <- vcf@gt[,1]
#    gt <- temp_vcf@gt[, -1]
#    cols <- gt[ , which(names(vcf_list[i]) == pop)]
#    temp_vcf@gt <- cbind(FORMAT, cols)
#    vcf_list[[i]] <- temp_vcf
#  }

  if(method == 'N'){
    
    myHeader <- c('N', length(vcf_list), nrow(vcf_list[[1]]))
        
    pop_list <- vector(mode = 'list', length=length(vcf_list))
    names(pop_list) <- names(vcf_list)
    
    # Extract alleles
    for(i in 1:length(vcf_list)){    
      gt <- extract.gt(vcf_list[[i]], return.alleles = T)

      allele1 <- apply(gt, MARGIN = 2, function(x){ substr(x, 1, 1) })
      rownames(allele1) <- NULL
      allele1 <- t(allele1)
      rownames(allele1) <- paste(rownames(allele1), "_1", sep = "")
      allele2 <- apply(gt, MARGIN = 2, function(x){ substr(x, 3, 3) })
      rownames(allele2) <- NULL
      allele2 <- t(allele2)
      rownames(allele2) <- paste(rownames(allele2), "_2", sep = "")
      pop_list[[i]][[1]] <- allele1
      pop_list[[i]][[2]] <- allele2
    }
    
    # Write to file
    write(myHeader, file = out_file, ncolumns = length(myHeader), sep = "\t")
    write(rep(1, times = ncol(pop_list[[1]][[1]])), file = out_file, ncolumns = ncol(pop_list[[1]][[1]]), append = TRUE, sep = "\t")
    for(i in 1:length(pop_list)){
      popName <- c(2*nrow(pop_list[[i]][[1]]), names(pop_list)[i])
      write(popName, file = out_file, ncolumns = length(popName), append = TRUE, sep = "\t")
      for(j in 1:ncol(pop_list[[i]][[1]])){
        utils::write.table(pop_list[[i]][[1]][,j], file = out_file, append = TRUE, quote = FALSE, 
                           sep = "\t", row.names = TRUE, col.names = FALSE)
        utils::write.table(pop_list[[i]][[2]][,j], file = out_file, append = TRUE, quote = FALSE, 
                           sep = "\t", row.names = TRUE, col.names = FALSE)
      }
    }
    
  } else if(method == 'H'){
    
    myHeader <- c('H', length(vcf_list), nrow(vcf_list[[1]]))
    
    # Summarize populations
    pop_list <- vector(mode = 'list', length=length(vcf_list))
    names(pop_list) <- names(vcf_list)
    for(i in 1:length(vcf_list)){
      # Matrix to hold the summary
      myMat <- matrix(nrow = nrow(vcf_list[[i]]), ncol = 6)
      
      # Population summary
      var_info <- as.data.frame(vcf_list[[i]]@fix[,1:2, drop = FALSE])
      var_info$mask <- TRUE
      gt <- extract.gt(vcf_list[[i]])
      popSum <- .gt_to_popsum(var_info = var_info, gt = gt)
#      popSum <- matrix(unlist(strsplit(as.character(popSum$Allele_counts), split = ",", fixed = TRUE)), ncol = 2, byrow = TRUE)
      
      # Populate matrix
      myMat[,1] <- paste(vcf_list[[i]]@fix[,'CHROM'], vcf_list[[i]]@fix[,'POS'], sep = "_")
      myMat[,2] <- vcf_list[[i]]@fix[,'REF']
      myMat[,4] <- vcf_list[[i]]@fix[,'ALT']
      myMat[,3] <- unlist(lapply(strsplit(as.character(popSum$Allele_counts), split = ",", fixed = TRUE), function(x){x[1]}))
      myMat[,3][is.na(myMat[,3])] <- 0
      myMat[,5] <- unlist(lapply(strsplit(as.character(popSum$Allele_counts), split = ",", fixed = TRUE), function(x){x[2]}))
      myMat[,5][is.na(myMat[,5])] <- 0
      myMat[,6] <- as.numeric(myMat[,3]) + as.numeric(myMat[,5])
      
      pop_list[[i]] <- myMat
    }
    
    # Write to file
    write(myHeader, file = out_file, ncolumns = length(myHeader), sep = "\t")
    #write(rep(1, times = nrow(pop_list[[1]])), file = out_file, ncolumns = nrow(pop_list[[1]]), append = TRUE, sep = "\t")
    for(i in 1:length(pop_list)){
      popName <- c(pop_list[[i]][1,6], names(pop_list[i]))
      write(popName, file = out_file, ncolumns = length(popName), append = TRUE, sep = "\t")
      utils::write.table(pop_list[[i]], file = out_file, append = TRUE, quote = FALSE, 
                         sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
  } else {
    stop("You should never get here!")
  }
  
  return( invisible(NULL) )
}

