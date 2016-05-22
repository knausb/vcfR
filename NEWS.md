

# vcfR 1.1.0
Released on CRAN 2016-05-XX.

This release includes the incorporation of suggestions made by reviewers of the manuscript submitted to Molecular Ecology Resources.

* added `is.het()` to identify heterozygotes in a matrix of genotypes.
* Fixed one-off error in `vcfR2DNAbin` where a variant one position beyond the locus would attempt to be included but threw an error.
* Added examples to VCF input and output.
* Added `vcfR_test` as lightweight test VCF data.
* Changed chromR@name to chromR@names for consistency with other R objects.
* Added `AD_frequency` calculates allele frequencies from matrices of AD data.
* `read.vcfR()` handles VCF data with no GT region (ala LoFreq).
* `gt2alleles` handles missing data ('.').
* `read.vcfR()` checks for and removes carriage returns (Windows).
* `vcfR2DNAbin` converts 'NA' to 'n' prior to conversion to DNAbin.
* `chromR2vcfR` implements use.mask.
* `extract.gt()` converts "." to NA.
* Added tidyr compatibility - thank you Eric Anderson!
* `write.vcf()` now uses mask = TRUE.
* `maf()` provides counts and frequency for the minor (or other) allele.
* `create.chromR()` now handles instances with no seq and the annotation position exceeds the greatest VCF POS.
* `read.vcfR()` now handles tilde expansion.
* `addID()` populates the non-missing values in the ID column of VCF data by concatenating the chromosome and position. 

# vcfR 1.0.0
Released on CRAN 2016-02-22.
This release was used to prepare the manuscript that was submitted to Molecular Ecology Resources.

