

# vcfR 1.0.0.9000

* changed chromR@name to chromR@names for consistency with other R objects.
* `added AD_frequency` calculates allele frequencies from matrices of AD data.
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

