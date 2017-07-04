#' @title Convert vcfR object to MigrateN format
#' @description The function converts vcfR format data to text format that can be used as infile for MigrateN.
#'
#' @param vcf a vcfR object.
#' @param pop factor of population.
#' @param in_pop population of interest to report in output file.
#' @param outfile name of output file.
#' @return a text file that can be used as an input for MigrateN software (SNP format).
#'
#' @details
#' This function converts the vcfR object to a text file which can be used downstream for MigrateN.
#' The function will remove loci with missing data, indels, and multiallelic loci. Thus, only SNPs data are analysed
#' where length of each locus is 1.
#'
#'
#' @examples
#' library(vcfR)
#' pkg <- "pinfsc50"
#' my_vcf <- system.file("extdata", "pinf_sc50.vcf.gz", package = pkg)
#' my_vcf <- read.vcfR( my_vcf, verbose = FALSE )
#' my_pop <- as.factor(rep(c("A", "B", "C"), each = 6))
#' vcfR2migrate(vcf = my_vcf , pop = my_pop , in_pop = c("A","C"), out_file = "my2pop.txt")
#'
#'
#' @export
vcfR2migrate <- function(vcf, pop, in_pop, out_file) {
  library(vcfR)

  #Removing indels and multiallelic loci
  vcf <- extract.indels(vcf, return.indels = F)
  vcf <- vcf[is.biallelic(vcf)]

  gt <- extract.gt(vcf, return.alleles = T, convertNA = T)
  gt <- gt[!rowSums((is.na(gt))),]

  FORMAT <- vcf@gt[1:nrow(gt),1]
  vcf@gt <- cbind(FORMAT, gt)


  vcf_list <- vector("list", length(unique(pop)))
  names(vcf_list) <- unique(pop)

  for (i in (1:length(vcf_list))) {
    temp_pop <- names(vcf_list[i])
    temp_vcf <- vcf
    FORMAT <- vcf@gt[,1]
    gt <- temp_vcf@gt[, -1]

    cols <- gt[ , which(names(vcf_list[i]) == pop)]

    temp_vcf@gt <- cbind(FORMAT, cols)


    vcf_list[[i]] <- temp_vcf

  }

  n.pop <- vector("list", length(vcf_list))
  names(n.pop) <- names(vcf_list)
  for ( i in 1:length(vcf_list)) {

    pop1 <- vcf_list[[i]]

    gt <- pop1@gt

    all <- vector("list", nrow(gt))

    for (j in 1:nrow(gt)){
      ind <- gt[j, -1]
      all[[j]] <- as.data.frame(unlist(strsplit(ind, "\\||/")))
    }

    n.pop[[i]] <- do.call("rbind", all)

  }


  migrate_pop <- n.pop

  for (i in 1:length(n.pop)) {
    temp.pop <- n.pop[[i]]
    names <- rownames(temp.pop)
    rownames(temp.pop) <- NULL
    temp.pop <- cbind(names,temp.pop)
    indnum <- nrow(temp.pop)/ (nrow(vcf@gt))
    colnames(temp.pop) <- c(indnum, names(n.pop)[[i]])

    migrate_pop[[i]] <- temp.pop



  }

  loci_length <- rep(1, nrow(vcf@gt))

  x <- list(loci_length = loci_length, migrate_pop = migrate_pop)

  in_pop <- in_pop
  output_path <- file.path(out_file)
  data_to_write <- x$migrate_pop[sapply(x$migrate_pop, function(x) any(colnames(x) %in% in_pop))]

  # write header
  header <- paste0("N ", length(data_to_write), " ", length(x$loci_length), "\n", paste(x$loci_length, collapse = " "))
  writeLines(header, output_path)

  # write data
  for (frame in data_to_write) {
    write.table(frame, file = output_path, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  }

}
