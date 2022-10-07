library(testthat)
library(vcfR)
context("vcfRtidy functions")

data("vcfR_example")
#data("vcfR_test")
data("vcfR_test_snpEff")

##### ##### ##### ##### #####
# vcf_field_names


test_that("vcf_field_names works",{
  data("vcfR_example")
  Z <- vcf_field_names(vcf, tag = "INFO")

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")

  Z <- vcf_field_names(vcf, tag = "FORMAT")

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")
})


test_that("vcf_field_names works, comma in quotes not parsed",{
   data("vcfR_test")
   myMeta <- vcfR_test@meta
   vcfR_test@meta <- c(myMeta[1:12], '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">', myMeta[13:18])

   Z <- vcf_field_names(vcfR_test, tag = "INFO")
   expect_is(Z, "tbl_df")
   expect_is(Z, "tbl")
   expect_is(Z, "data.frame")

   vcfR_test@meta[13] <- '##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">'
   Z <- vcf_field_names(vcfR_test)

   expect_is(Z, "tbl_df")
   expect_is(Z, "tbl")
   expect_is(Z, "data.frame")
})


test_that("vcf_field_names works, zero INFO records in meta",{
  data("vcfR_test")
  ## cause error ##
  # remove all INFO lines in @meta object
  INFO.meta.lines <- grepl("^##INFO", vcfR_test@meta);
  vcfR_test@meta <- vcfR_test@meta[!INFO.meta.lines];

  #debug(vcfR2tidy)
  #debug(vcf_field_names)
  tidyVCF <- vcfR2tidy(vcfR_test)

  expect_is(tidyVCF, "list")
  expect_identical(names(tidyVCF), c("fix", "gt", "meta"))
  expect_identical(
    names(tidyVCF$fix),
    c("ChromKey", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  )
})



##### ##### ##### ##### #####
# extract_gt_tidy

test_that("extract_gt_tidy works for GT element",{
#  Z <- extract_gt_tidy(vcf)
  suppressMessages(
#
    Z <- extract_gt_tidy(vcf, format_fields = c('GT'))
#    Z <- extract_gt_tidy(vcf, format_fields = c('GT'), format_types = TRUE )
#    Z <- extract_gt_tidy( vcf, format_fields = c('GT'), format_types = character(0) )
  )
  expect_is(Z, 'tbl_df')
  expect_equal(names(Z)[1], 'Key')
  expect_equal(names(Z)[2], 'Indiv')
  expect_equal(names(Z)[3], 'gt_GT')
  expect_equal(names(Z)[4], 'gt_GT_alleles')
})


test_that("extract_gt_tidy works for all elements",{
  suppressMessages( Z <- extract_gt_tidy(vcf) )
  expect_is(Z, 'tbl_df')

})


##### ##### ##### ##### #####
# vcfR2tidy

test_that("vcfR2tidy works",{
  data("vcfR_test")
  Z <- vcfR2tidy(vcfR_test, info_only = FALSE)

  expect_is(Z, 'list')
  expect_equal( length(Z), 3 )

  expect_is(Z[['fix']], "tbl_df")
  expect_is(Z[['fix']], "tbl")
  expect_is(Z[['fix']], "data.frame")

})

# test_that("vcfR2tidy works, ID=REF",{
#    data("vcfR_test")
#    myMeta <- vcfR_test@meta
#    vcfR_test@meta <- c(myMeta[1:12], '##INFO=<ID=REF,Number=0,Type=Flag,Description="Has reference A coding region variation where one allele in the set is identical to the reference sequence. FxnCode = 8">', myMeta[13:18])
#    Z <- vcfR2tidy(vcf, info_only = FALSE, verbose = FALSE)
#    Z$meta$ID
# })

##### ##### ##### ##### #####
# extract_info_tidy


test_that("extract_info_tidy works",{
  Z <- extract_info_tidy(vcf, info_fields = c("AC", "AN", "MQ"), info_types = c(AN = "i", MQ = "n"))

  expect_is(Z, "tbl_df")
  expect_is(Z, "tbl")
  expect_is(Z, "data.frame")

#  expect_equal( length(Z), 3 )

})

test_that("extract_info_tidy works with Flags",{
  data(vcfR_test)
  Z <- extract_info_tidy(vcfR_test, info_types = TRUE)
  expect_is(Z$DB, "logical")
  expect_is(Z$H2, "logical")
})

##### ##### ##### ##### #####
# tidy field separation (column-wise) for snpEff and VEP annotated VCFs
data("vcfR_test_snpEff")
vcfR_test_snpEff_tidy <- vcfR2tidy(vcfR_test_snpEff)
data("vcfR_test_VEP")
vcfR_test_VEP_tidy <- vcfR2tidy(vcfR_test_VEP)

test_that("extraction of annotation fields works for snpEff", {
  ann_cols_test <- c(
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
    "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos / cDNA.length", "CDS.pos / CDS.length", "AA.pos / AA.length",
    "Distance", "ERRORS / WARNINGS / INFO"
  )
  lof_cols_test <- c("lof_Gene_Name", "lof_Gene_ID", "lof_Number_of_transcripts_in_gene", "lof_Percent_of_transcripts_affected")
  nmd_cols_test <- c("nmd_Gene_Name", "nmd_Gene_ID", "nmd_Number_of_transcripts_in_gene", "nmd_Percent_of_transcripts_affected")

  snpEff_ann_cols <- extract_ann_cols_snpeff(vcfR_test_snpEff_tidy$meta, clean_names = FALSE)
  expect_identical(snpEff_ann_cols$ANN, ann_cols_test)
  expect_identical(snpEff_ann_cols$LOF, lof_cols_test)
  expect_identical(snpEff_ann_cols$NMD, nmd_cols_test)

  ann_cols_test <- janitor::make_clean_names(ann_cols_test)
  snpEff_ann_cols_clean <- extract_ann_cols_snpeff(vcfR_test_snpEff_tidy$meta, clean_names = TRUE)
  expect_identical(snpEff_ann_cols_clean$ANN, ann_cols_test)

  # Test the custom definition of column names.
  ann_cols_test[ann_cols_test == "gene_id"] <- "ensembl_id"
  snpEff_ann_cols_clean_custom <- extract_ann_cols_snpeff(
    vcfR_test_snpEff_tidy$meta, clean_names = TRUE, custom_names = list("gene_id" = "ensembl_id")
  )
  expect_identical(snpEff_ann_cols_clean_custom$ANN, ann_cols_test)

  # Test the case an annotation column is missing.
  vcfR_test_snpEff_tidy_filtered <- vcfR_test_snpEff_tidy
  vcfR_test_snpEff_tidy_filtered$fix <- dplyr::select(vcfR_test_snpEff_tidy_filtered$fix, -LOF)
  vcfR_test_snpEff_tidy_filtered$meta <- dplyr::filter(vcfR_test_snpEff_tidy_filtered$meta, ID != "LOF")
  snpEff_ann_cols <- extract_ann_cols_snpeff(vcfR_test_snpEff_tidy_filtered$meta)
  expect_true(all(c("ANN", "NMD") %in% names(snpEff_ann_cols)))
})

test_that("extraction of annotation fields works for VEP", {
  csq_cols_test <- c(
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene",
    "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc",
    "HGVSp", "cDNA_position", "CDS_position", "Protein_position",
    "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND",
    "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "CCDS",
    "ENSP", "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM",
    "GENE_PHENO", "SIFT", "PolyPhen", "DOMAINS", "miRNA", "AF", "AFR_AF",
    "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "AA_AF", "EA_AF", "gnomAD_AF",
    "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF",
    "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF",
    "MAX_AF", "MAX_AF_POPS", "FREQS", "CLIN_SIG", "SOMATIC", "PHENO",
    "PUBMED", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS"
  )

  VEP_csq_cols <- extract_ann_cols_vep(vcfR_test_VEP_tidy$meta, clean_names = FALSE)
  expect_identical(VEP_csq_cols$CSQ, csq_cols_test)
})

split_ann <- function(ann, cols) {
  ann %>%
    # Split multiple annotations.
    stringr::str_split(",") %>%
    unlist() %>%
    # Split annotation fields.
    stringr::str_split("\\|") %>%
    data.frame() %>%
    t() %>%
    magrittr::set_colnames(cols) %>%
    tibble::as_tibble()
}

test_that("separation of annotation columns in tidy VCF from snpEff works", {
  snpEff_ann_cols <- extract_ann_cols_snpeff(vcfR_test_snpEff_tidy$meta)
  # We will just test the first multi-annotation, which will be separated into two rows.
  snpEff_sep_test <- split_ann(vcfR_test_snpEff_tidy$fix$ANN[[1]], snpEff_ann_cols$ANN)
  snpEff_sep <- separate_ann_snpeff(vcfR_test_snpEff_tidy$fix, vcfR_test_snpEff_tidy$meta)
  expect_identical(
    snpEff_sep_test,
    snpEff_sep %>%
      dplyr::slice(1:2) %>%
      dplyr::select(dplyr::all_of(snpEff_ann_cols$ANN))
  )

  for (col in c("LOF", "NMD")) {
    snpEff_sep_ <- snpEff_sep %>%
      dplyr::filter(!is.na(!!rlang::sym(col))) %>%
      dplyr::slice(1)
    snpEff_sep_test <- vcfR_test_snpEff_tidy$fix %>%
      dplyr::filter(!is.na(!!rlang::sym(col))) %>%
      {.[[col]][[1]]} %>%
      stringr::str_remove("^\\(") %>%
      stringr::str_remove("\\)$") %>%
      stringr::str_split("\\|") %>%
      unlist() %>%
      data.frame() %>%
      t() %>%
      magrittr::set_colnames(snpEff_ann_cols[[col]]) %>%
      tibble::as_tibble()
    expect_identical(snpEff_sep_test, dplyr::select(snpEff_sep_, dplyr::all_of(snpEff_ann_cols[[col]])))
  }
})

test_that("separation of annotation columns in tidy VCF from VEP works", {
  VEP_ann_cols <- extract_ann_cols_vep(vcfR_test_VEP_tidy$meta)
  vcfR_test_VEP_tidy_multi <- vcfR_test_VEP_tidy$fix %>%
    dplyr::filter(stringr::str_detect(CSQ, stringr::fixed(","))) %>%
    dplyr::slice(1)
  VEP_sep_test <- split_ann(vcfR_test_VEP_tidy_multi$CSQ[[1]], VEP_ann_cols$CSQ)
  VEP_sep <- separate_ann_vep(vcfR_test_VEP_tidy_multi, vcfR_test_VEP_tidy$meta)
  expect_identical(VEP_sep_test, dplyr::select(VEP_sep, dplyr::all_of(VEP_ann_cols$CSQ)))
})
