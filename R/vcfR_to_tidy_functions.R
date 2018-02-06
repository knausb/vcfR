



#### Document all the vcf2tidy related functions together ####
#' @title Convert vcfR objects to tidy data frames
#' @name Convert to tidy data frames
#' @rdname vcfR_to_tidy_conversion
#' @description
#' Convert the information in a vcfR object to a long-format data frame
#' suitable for analysis or use with Hadley Wickham's packages, 
#' \href{https://cran.r-project.org/package=dplyr}{dplyr},
#' \href{https://cran.r-project.org/package=tidyr}{tidyr}, and
#' \href{https://cran.r-project.org/package=ggplot2}{ggplot2}.
#' These packages have been
#' optimized for operation on large data frames, and, though they can bog down
#' with very large data sets, they provide a good framework for handling and filtering
#' large variant data sets.  For some background
#' on the benefits of such "tidy" data frames, see 
#' \href{https://www.jstatsoft.org/article/view/v059i10}{this article}.
#' 
#' For some filtering operations, such as those where one wants to filter genotypes
#' upon GT fields in combination with INFO fields, or more complex 
#' operations in which one wants to filter
#' loci based upon the number of individuals having greater than a certain quality score,
#' it will be advantageous to put all the information into a long format data frame 
#' and use \code{dplyr} to perform the operations.  Additionally, a long data format is
#' required for using \code{ggplot2}.  These functions convert vcfR objects to long format
#' data frames.  
#'  
#' @param x an object of class vcfR
#' 
#' @details 
#' The function \strong{vcfR2tidy} is the main function in this series.  It takes a vcfR
#' object and converts the information to a list of long-format data frames.  The user can
#' specify whether only the INFO or both the INFO and the FORMAT columns should be extracted, and also
#' which INFO and FORMAT fields to extract.  If no specific INFO or FORMAT fields are asked
#' for, then they will all be returned.  If \code{single_frame == FALSE} and 
#' \code{info_only == FALSE} (the default), 
#' the function returns a list with three components: \code{fix}, \code{gt}, and \code{meta} as follows:
#' \enumerate{
#' \item \code{fix} A data frame of the fixed information columns and the parsed INFO columns, and 
#' an additional column, \code{ChromKey}---an integer identifier
#' for each locus, ordered by their appearance in the original data frame---that serves
#' together with POS as a key back to rows in \code{gt}.  
#' \item \code{gt} A data frame of the genotype-related fields. Column names are the names of the 
#' FORMAT fields with \code{gt_column_prepend} (by default, "gt_") prepended to them.  Additionally
#' there are columns \code{ChromKey}, and \code{POS} that can be used to associate
#' each row in \code{gt} with a row in \code{fix}.
#' \item\code{meta} The meta-data associated with the columns that were extracted from the INFO and FORMAT
#' columns in a tbl_df-ed data frame.  
#' }
#' This is the default return object because it might be space-inefficient to
#' return a single tidy data frame if there are many individuals and the CHROM names are
#' long and/or there are many INFO fields.  However, if
#' \code{single_frame = TRUE}, then the results are returned as a list with component \code{meta}
#' as before, but rather than having \code{fix} and \code{gt} as before, both those data frames
#' have been joined into component \code{dat} and a ChromKey column is not returned, because
#' the CHROM column is available.
#' 
#' If \code{info_only == FALSE}, then just the fixed columns and the parsed INFO columns are 
#' returned, and the FORMAT fields are not parsed at all.  The return value is a list with
#' components \code{fix} and \code{meta}.  No column ChromKey appears.
#' 
#' The following functions are called by \strong{vcfR2tidy} but are documented below because
#' they may be useful individually.
#' 
#' The function \strong{extract_info_tidy} let's you pass in a vector of the INFO fields that
#' you want extracted to a long format data frame. If you don't tell it which fields to 
#' extract it will extract all the INFO columns detailed in the VCF meta section.
#' The function returns a tbl_df data frame of the INFO fields along with with an additional
#' integer column \code{Key} that associates
#' each row in the output data frame with each row (i.e. each CHROM-POS combination) 
#' in the original vcfR object \code{x}.  
#' 
#' The function \strong{extract_gt_tidy} let's you pass in a vector of the FORMAT fields that
#' you want extracted to a long format data frame. If you don't tell it which fields to 
#' extract it will extract all the FORMAT columns detailed in the VCF meta section.
#' The function returns a tbl_df data frame of the FORMAT fields with an additional
#' integer column \code{Key} that associates
#' each row in the output data frame with each row (i.e. each CHROM-POS combination),
#' in the original vcfR object \code{x}, and an additional column \code{Indiv} that gives
#' the name of the individual.  
#' 
#' The function \strong{vcf_field_names} is a helper function that
#' parses information from the metadata section of the
#' VCF file to return a data frame with the \emph{metadata} information about either the INFO 
#' or FORMAT tags.  It
#' returns a \code{tbl_df}-ed data frame with column names: "Tag", "ID", "Number","Type",
#' "Description", "Source", and "Version".
#' 
#' @return An object of class tidy::data_frame or a list where every element is of class tidy::data_frame.
#' 
#' @note  To run all the examples, you can issue this:
#' \code{example("vcfR2tidy")}
#' 
#' @author Eric C. Anderson <eric.anderson@@noaa.gov>
#' @seealso
#' \href{https://cran.r-project.org/package=dplyr}{dplyr},
#' \href{https://cran.r-project.org/package=tidyr}{tidyr}.
#' 
#' @examples 
#' # load the data
# data(vcfR_example)
#' data("vcfR_test")
#' vcf <- vcfR_test
#' 
#' 
#' # extract all the INFO and FORMAT fields into a list of tidy
#' # data frames: fix, gt, and meta. Here we don't coerce columns
#' # to integer or numeric types...
#' Z <- vcfR2tidy(vcf)
#' names(Z)
#' 
#' 
#' # here is the meta data in a table
#' Z$meta
#' 
#' 
#' # here is the fixed info
#' Z$fix
#' 
#' 
#' # here are the GT fields.  Note that ChromKey and POS are keys
#' # back to Z$fix
#' Z$gt
#' 
#' 
#' # Note that if you wanted to tidy this data set even further
#' # you could break up the comma-delimited columns easily
#' # using tidyr::separate
#' 
#' 
#' 
#' 
#' # here we put the data into a single, joined data frame (list component
#' # dat in the returned list) and the meta data.  Let's just pick out a 
#' # few fields:
#' vcfR2tidy(vcf, 
#'           single_frame = TRUE, 
#'           info_fields = c("AC", "AN", "MQ"), 
#'           format_fields = c("GT", "PL"))
#'
#' 
#' # note that the "gt_GT_alleles" column is always returned when any
#' # FORMAT fields are extracted.
#' 
#' 
#' 
#' 
#' # Here we extract a single frame with all fields but we automatically change
#' # types of the columns according to the entries in the metadata.
#' vcfR2tidy(vcf, single_frame = TRUE, info_types = TRUE, format_types = TRUE)
#' 
#' 
#' 
#' 
#' # for comparison, here note that all the INFO and FORMAT fields that were
#' # extracted are left as character ("chr" in the dplyr summary)
#' vcfR2tidy(vcf, single_frame = TRUE)
#' 
#' 
#' 
#' 
#' 
#' # Below are some examples with the vcfR2tidy "subfunctions"
#' 
#' 
#' # extract the AC, AN, and MQ fields from the INFO column into
#' # a data frame and convert the AN values integers and the MQ
#' # values into numerics.
#' extract_info_tidy(vcf, info_fields = c("AC", "AN", "MQ"), info_types = c(AN = "i", MQ = "n"))
#' 
#' # extract all fields from the INFO column but leave 
#' # them as character vectors
#' extract_info_tidy(vcf)
#' 
#' # extract all fields from the INFO column and coerce 
#' # types according to metadata info
#' extract_info_tidy(vcf, info_types = TRUE)
#' 
#' # get the INFO field metadata in a data frame
#' vcf_field_names(vcf, tag = "INFO")
#' 
#' # get the FORMAT field metadata in a data frame
#' vcf_field_names(vcf, tag = "FORMAT")
#' 
#' 
#' 




#### vcfR2tidy ####
#' @rdname vcfR_to_tidy_conversion
#' @aliases vcfR2tidy
#' 
#' @param info_only if TRUE return a list with only a \code{fix} component
#' (a single data frame that has the parsed INFO information) and 
#' a \code{meta} component. Don't extract any of the FORMAT fields. 
#' @param single_frame return a single tidy data frame in list component
#' \code{dat} rather returning it in components
#' \code{fix} and/or \code{gt}. 
#' @param toss_INFO_column if TRUE (the default) the INFO column will be removed from output as
#' its consituent parts will have been parsed into separate columns.
#' @param ... more options to pass to \code{\link{extract_info_tidy}} and 
#' \code{\link{extract_gt_tidy}}.  See parameters listed below.
#' 
# @importFrom dplyr everything
# @import dplyr
#' 
#' @export
vcfR2tidy <- function(x, 
                      info_only = FALSE, 
                      single_frame = FALSE, 
                      toss_INFO_column = TRUE,
                      ...) {
  
#### Some Error Checking and Preliminaries ####
  if(single_frame == TRUE && info_only == TRUE) 
    stop("You can pass both single_frame and info_only as TRUE")
  
  #check to make sure that the user didn't pass in unacceptable params in ...
  dotslist <- list(...)
  unk_parm <- setdiff(
    names(dotslist), 
    c("info_fields", "info_types", "info_sep", "format_fields", "format_types", "dot_is_NA",
      "alleles", "allele.sep", "gt_column_prepend", "verbose")
  )
  
  if(length(unk_parm) > 0){
    stop("Unknown \"...\" parameters ", 
        paste(unk_parm, collapse = " "), 
        " to function vcfR2tidy"
        )
  }
  
  info_dots <- dotslist[names(dotslist) %in% c("info_fields", "info_types", "info_sep")]
  info_dots$x = x
  format_dots <- dotslist[names(dotslist) %in% c("format_fields", "format_types", "dot_is_NA",
                                                 "alleles", "allele.sep", "gt_column_prepend", "verbose")]
  format_dots$x = x
  
  # klugie hack for dealing with the gt_column_prepend
  if(!is.null(format_dots[["gt_column_prepend"]])) {
    gt_prep <- format_dots[["gt_column_prepend"]]
  } else {
    gt_prep = "gt_"
  }
  
#### extract the INFO data. and return if that is all the is requested ####
  # get the base fix data as a data frame
  base <- as.data.frame(x@fix, stringsAsFactors = FALSE) %>% tibble::as.tibble()
  base$POS <- as.integer(base$POS)
  base$QUAL <- as.numeric(base$QUAL)
  if(toss_INFO_column == TRUE) {
    base <- base %>% dplyr::select_(~ -INFO)
  }
  
  # also get the  full meta data for all the INFO fields
  info_meta_full <- vcf_field_names(x, tag = "INFO") 
  
  fix <- do.call(what = extract_info_tidy, args = info_dots)
  if(info_only == TRUE) {
#    ret <- cbind(base, fix) %>% 
    ret <- dplyr::bind_cols(base, fix) %>% 
      tibble::as.tibble() %>%
      dplyr::select_(~ -Key)
    
    # only retain meta info for the fields that we are returning
    info_meta <- info_meta_full %>%
      dplyr::filter_(~ID %in% names(ret))
    
    return(list(fix = ret, meta = info_meta))
  }
  
#### Extract the GT data, and return what is appropriate ####
  # if you got here then we need to extract some gt fields, too
  gt <- do.call(what = extract_gt_tidy, args = format_dots)
  
  # get the full FORMAT meta data and add the gt_column_prepend to them
  gt_meta_full <- vcf_field_names(x, tag = "FORMAT") %>%
    dplyr::mutate_(ID = ~paste(gt_prep, ID, sep = ""))
  
  # if the user is asking for a single data frame we give it to them here:
  if(single_frame == TRUE) {
#    ret <- cbind(base, fix) %>%
    ret <- dplyr::bind_cols(base, fix) %>%
      dplyr::left_join(gt, by = "Key") %>%
      tibble::as.tibble() %>%
      dplyr::select_(~ -Key)  # no point in keeping Key around at this point
    
    info_meta <- info_meta_full %>%
      dplyr::filter_(~ID %in% names(ret))
    
    gt_meta <-  gt_meta_full %>%
      dplyr::filter_(~ID %in% names(ret))
    
      return(list(dat = ret, meta = dplyr::bind_rows(info_meta, gt_meta)))
  }
  
  # if the user is not asking for a single data frame then we return a list 
  # which has appropriate keys for getting the fix and the gt associated
  # appropriately.
#  retfix <- cbind(base, fix) %>%
  retfix <- dplyr::bind_cols(base, fix) %>%
    tibble::as.tibble() %>%
    dplyr::mutate_(ChromKey = ~as.integer(factor(CHROM), levels = unique(CHROM))) %>%
    dplyr::select_(~ChromKey, ~dplyr::everything())  # note that we will drop Key from this after we have used it
  
  retgt <- gt %>%
    dplyr::left_join(dplyr::select_(retfix, ~ChromKey, ~Key, ~POS), by = "Key") %>%
    dplyr::select_(~ChromKey, ~POS, ~dplyr::everything()) %>%
    dplyr::select_(~ -Key)
  
  info_meta <- info_meta_full %>%
    dplyr::filter_(~ID %in% names(retfix))
  
  gt_meta <-  gt_meta_full %>%
    dplyr::filter_(~ID %in% names(retgt))
    
  
  # return the list
  list(
    fix = retfix %>% 
      dplyr::select_(~ -Key),
    gt = retgt,
    meta = dplyr::bind_rows(info_meta, gt_meta)
  )
}


#### extract_info_tidy ####
#' @rdname vcfR_to_tidy_conversion
#' @aliases extract_info_tidy
#' 
#' @param info_fields names of the fields to be extracted from the INFO column
#' into a long format data frame.  If this is left as NULL (the default) then
#' the function returns a column for every INFO field listed in the metadata.
#' @param info_types named vector of "i" or "n" if you want the fields extracted from the INFO column to be converted to integer or numeric types, respectively.
#' When set to NULL they will be characters.  
#' The names have to be the exact names of the fields.  
#' For example \code{info_types = c(AF = "n", DP = "i")} will convert column AF to numeric and DP to integer.
#' If you would like the function to try to figure out the conversion from the metadata information, then set \code{info_types = TRUE}.  
#' Anything with Number == 1 and (Type == Integer or Type == Numeric) will then be converted accordingly.
#' @param info_sep the delimiter used in the data portion of the INFO fields to 
#' separate different entries.  By default it is ";", but earlier versions of the VCF
#' standard apparently used ":" as a delimiter.
#' @export
extract_info_tidy <- function(x, info_fields = NULL, info_types = TRUE, info_sep = ";") {
  
  if(!is.null(info_fields) && any(duplicated(info_fields))) stop("Requesting extraction of duplicate info_field names")
  if(class(x) != "vcfR") stop("Expecting x to be a vcfR object, not a ", class(x))
  
  vcf <- x
  x <- as.data.frame(x@fix, stringsAsFactors = FALSE) %>% 
    tibble::as.tibble()
  
  # if info_fields is NULL then we try to do all of them
  if(is.null(info_fields)) {
    info_df <- vcfR::vcf_field_names(vcf, tag = "INFO")
    info_fields <- info_df$ID
  }
  # if info_types == TRUE
  # then we try to discern the fields amongst info_fields that should be coerced to integer and
  # numeric
  if(!is.null(info_types) && length(info_types) == 1 && info_types[1] == TRUE) {
    info_df <- vcfR::vcf_field_names(vcf, tag = "INFO") %>%
      dplyr::filter_(~ID %in% info_fields)
    info_types <- guess_types(info_df)
  }
  
  # here is where the action is
  # first split into a list of vectors and then make them named vectors of values and 
  # pick them out in order using info_fields
  ret <- stringr::str_split(string = x$INFO, pattern = info_sep) %>%
    lapply(function(x) {
      y <- stringr::str_split(x, pattern = "=", n = 2)
      vals <- unlist(lapply(y, function(z) z[2]))
      names(vals) <- unlist(lapply(y, function(z) z[1]))
      unname(vals[info_fields])
    }) %>% 
    unlist %>%
    matrix(ncol = length(info_fields), byrow = TRUE) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(info_fields) %>%
    tibble::as.tibble()
  
  if(!is.null(info_types)) {
    ns <- info_types[!is.na(info_types) & info_types == "n"]
    is <- info_types[!is.na(info_types) & info_types == "i"]
    fs <- info_types[!is.na(info_types) & info_types == "f"]
    
    if(length(ns) > 0) {
      ret[names(ns)] <- lapply(ret[names(ns)], as.numeric)
    }
    if(length(is) > 0) {
      ret[names(is)] <- lapply(ret[names(is)], as.integer)
    }
    
    proc_flag <- function(x, INFO){
      x2 <- rep(FALSE, times = length(INFO))
      x2[ grep(x, INFO) ] <- TRUE
      x2
    }
    if(length(fs) > 0) {
      ret[names(fs)] <- lapply(names(fs), proc_flag, x$INFO)
    }
    
  }
  cbind(Key = 1:nrow(ret), ret) %>% tibble::as.tibble()
}

#### extract_gt_tidy ####
#' @rdname vcfR_to_tidy_conversion
#' @aliases extract_gt_tidy
#' 
#' @param format_fields names of the fields in the FORMAT column to be extracted from 
#' each individual in the vcfR object into 
#' a long format data frame.  If left as NULL, the function will extract all the FORMAT
#' columns that were documented in the meta section of the VCF file.
#' @param format_types named vector of "i" or "n" if you want the fields extracted according to the FORMAT column to be converted to integer or numeric types, respectively.
#' When set to TRUE an attempt to determine their type will be made from the meta information.
#' When set to NULL they will be characters.  
#' The names have to be the exact names of the format_fields.  
#' Works equivalently to the \code{info_types} argument in 
#' \code{\link{extract_info_tidy}}, i.e., if you set it to TRUE then it uses the information in the
#' meta section of the VCF to coerce to types as indicated.
#' @param dot_is_NA if TRUE then a single "." in a character field will be set to NA.  If FALSE
#' no conversion is done.  Note that "." in a numeric or integer field 
#' (according to format_types) with Number == 1 is always
#' going to be set to NA.
#' @param alleles if TRUE (the default) then this will return a column, \code{gt_GT_alleles} that
#' has the genotype of the individual expressed as the alleles rather than as 0/1.
#' @param allele.sep	character which delimits the alleles in a genotype (/ or |) to be passed to
#' \code{\link{extract.gt}}. Here this is not used for a regex (as it is in other functions), but merely
#' for output formatting.
#' @param gt_column_prepend string to prepend to the names of the FORMAT columns
#' @param verbose logical to specify if verbose output should be produced
#' in the output so that they
#' do not conflict with any INFO columns in the output.  Default is "gt_". Should be a 
#' valid R name. (i.e. don't start with a number, have a space in it, etc.)
#' @export
extract_gt_tidy <- function(x, 
                            format_fields = NULL, 
                            format_types = TRUE, 
                            dot_is_NA = TRUE,
                            alleles = TRUE,
                            allele.sep = "/",
                            gt_column_prepend = "gt_",
                            verbose = TRUE) {
  
  if(!is.null(format_fields) && any(duplicated(format_fields))){
    stop("Requesting extraction of duplicate format_field names")
  }
  if(class(x) != "vcfR"){
    stop("Expecting x to be a vcfR object, not a ", class(x))
  }
  
  vcf <- x  # Rename it.
  
  # Get this, because we may need it.
  # Extracts FORMAT acronyms from the meta region.
  format_df <- vcfR::vcf_field_names(vcf, tag = "FORMAT")
  
  # If format_fields is NULL then we try to do all of them
  if(is.null(format_fields)) {
    format_fields <- format_df$ID
  }
  # If info_types == TRUE
  # then we try to discern the fields amongst info_fields that should be coerced to integer and
  # numeric.
  if(!is.null(format_types) && length(format_types) == 1 && format_types[1] == TRUE) {
    format_types <- guess_types(format_df %>% dplyr::filter_(~ID %in% format_fields))
  }
  
  # Make a parallel vector that indicates which fields should be numeric or not
  # so we can tell extract.gt to take care of it.
  coerce_numeric <- rep(FALSE, length(format_fields))
  coerce_numeric[format_fields %in% names(format_types)] <- TRUE
  
  # Now get all the gt fields
  ex <- 1:length(format_fields)
  names(ex) <- format_fields
  
  get_gt <- function(i, ...){
    if(verbose == TRUE){
      message("Extracting gt element ", names(ex)[i])
    }
    ret <- extract.gt(x = vcf, element = format_fields[i], as.numeric = coerce_numeric[i])
    ret <- as.vector(ret)

    ret
  }
  geno_info <- lapply(ex, get_gt)
  
  geno_info <- dplyr::as_data_frame(geno_info)
  geno_info <- dplyr::mutate_(Key = ~rep(1:nrow(vcf@fix), times = ncol(vcf@gt) - 1),
                              Indiv = ~rep(colnames(vcf@gt)[-1], each = nrow(vcf@fix)),
                              geno_info
                              )
  geno_info <- dplyr::select_(geno_info, ~Key, ~Indiv, ~dplyr::everything())
#  geno_info <- dplyr::select_(geno_info, ~Key, ~Indiv, ~everything())
#  geno_info <- dplyr::select_(geno_info, ~Key, ~Indiv, grep(c('Key|Indiv'), names(Z), invert=TRUE))
  
#  geno_info <- lapply(ex, function(i) {
#    message("Extracting gt element ", names(ex)[i])
#    ret <- extract.gt(x = vcf, element = format_fields[i], as.numeric = coerce_numeric[i])
#    if(dot_is_NA == TRUE) ret[ret == "."] <- NA 
#    as.vector(ret)
#  }) %>%
#    dplyr::as_data_frame() %>%
    #setNames(paste("gt_", names(.), sep = "")) %>%
#    dplyr::mutate_(Key = ~rep(1:nrow(vcf@fix), times = ncol(vcf@gt) - 1),
 #          ChromKey = rep(fix$ChromKey, times = ncol(V@gt) - 1),
#           Indiv = ~rep(colnames(vcf@gt)[-1], each = nrow(vcf@fix))) %>%
#    dplyr::select_(~Key, ~Indiv, ~everything())

  
    
  # now coerce numerics that should be integers to ints:
#  if( length(format_types) > 0 ){
  if( sum( format_types == "i" ) > 0 ){
  geno_info[names(format_types)[format_types == "i"]] <- 
    lapply(geno_info[names(format_types)[format_types == "i"]], as.integer)
  }
  
  # and now, if alleles == TRUE, get the GT column expressed as alleles
  if(alleles == TRUE) {
    geno_info$GT_alleles <- as.vector(extract.gt(x = vcf, element = "GT", return.alleles = TRUE))
  }
  # now prepend gt_ to every column name except the Key and Indiv columns:
  names(geno_info)[-c(1,2)] <- paste(gt_column_prepend, names(geno_info)[-c(1,2)], sep = "")
  
  geno_info
}



# given a data frame of INFO or FORMAT info, this sets everything
# that has Number = 1 to Integer or Numeric as appropriate.  Returns
# a named vector suitable for passing to, for example, info_types.
# this is not exported
guess_types <- function(D) {
  tmp <- D %>%
    dplyr::filter_(~Number == 1) %>%
    dplyr::mutate_(tt = ~ifelse(Type == "Integer", "i", ifelse(Type == "Numeric" | Type == "Float", "n", ""))) %>%
    dplyr::filter_(~tt %in% c("n", "i")) %>%
    dplyr::select_(~ID, ~Number, ~Type, ~tt)
  
  tmp <- D %>%  dplyr::filter_(~Number == 0 & Type == 'Flag')  %>%
    dplyr::mutate_(tt = ~ifelse(Type == "Flag", "f")) %>%
    dplyr::filter_(~tt %in% c("f")) %>%
    dplyr::select_(~ID, ~Number, ~Type, ~tt)  %>%
    dplyr::bind_rows(tmp)   
  
  ret <- tmp$tt
  names(ret) <- tmp$ID
  ret
}

#### vcf_field_names ####
#' @rdname vcfR_to_tidy_conversion
#' @aliases vcf_field_names
#' 
#' @param tag name of the lines in the metadata section of the VCF file to parse out.
#' Default is "INFO".  The only other one tested and supported, currently is, "FORMAT".
#' 
#' @export
vcf_field_names <- function(x, tag = "INFO") {
  if(class(x) != "vcfR") stop("Expecting x to be a vcfR object, not a ", class(x))
  if( tag != 'INFO' & tag != 'FORMAT') stop("Expecting tag to either be INFO or FORMAT")

  # Subset to tag.
  x <- x@meta
  left_regx <- paste("^##", tag, "=<", sep = "")  # regex to match and replace 
  x <- x[grep(left_regx, x)]
  # Clean up the string ends.
  x <- sub(left_regx, "", x)
  x <- sub(">$", "", x)

  # Delimit on quote protected commas.
  x <- lapply(x, function(x){scan(text=x, what="character", sep=",", quiet = TRUE)})

  # Get unique keys.
  myKeys <- unique(unlist(lapply(strsplit(unlist(x), split = "="), function(x){x[1]})))
  # Omit default keys so we can make them first.
  myKeys <- grep("^ID$|^Number$|^Type$|^Description$", myKeys, invert = TRUE, value = TRUE)
  myKeys <- c("ID", "Number", "Type", "Description", myKeys)

  myReturn <- data.frame(matrix(ncol=length(myKeys) + 1, nrow=length(x)))
  colnames(myReturn) <- c("Tag", myKeys)
  myReturn[,'Tag'] <- tag
  getValue <- function(x){
    myValue <- grep(paste("^", myKeys[i], "=", sep=""), x, value = TRUE)
    if(length(myValue) == 0){
      is.na(myValue) <- TRUE
    } else {
      myValue <- sub(".*=", "", myValue)
    }
    myValue
  }
  
  for(i in 1:length(myKeys)){
    myReturn[,i+1] <- unlist(lapply(x, function(x){ getValue(x) }))
  }
  tibble::as_tibble(myReturn)
}


