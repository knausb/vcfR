
# vcfR 2.0.0.
There is currently no plan to release vcfR 2.0.0.
If and when this 'major' release occurs it will include changes that will break backward compatibility.
At the present, this is simply a to-do list for ideas to include in the next major release.

* Move 'FORMAT' column to its own slot. We can then cbind FORMAT and gt when passing to compiled code.
This may have been addressed at 64a308ba50b9119108e8946737460de5997b805b by adding `samples` to vcfR method `[`.


# vcfR 1.7.0
Released on CRAN 2018-02-07.
* `vcf_field_names()` now delimts on KEY= of key/value pairs, allows commas to be used within value.
* `read.vcfR()` will download files when provided with a link.
* Added example data from the Variant Effect Predictor (vep) `data(vep)`.


# vcfR 1.6.0
Released on CRAN 2017-12-08.
* `vcf_field_names()` now handles keys that are out of order and multiple optional keys.
* `vcfR2DNAbin()` can include indels and maintains alignment.
* `write.vcf()` now handles tilde expansion.
* `rePOS()` attempts to create a non-overlapping coordinate system from POS and CHROM.
* `vcfR2DNAbin()` manages the asterisk allele.
* `extract.indels()` ignores GATK's <NON_REF>.
* Added support for chromR objects with no gt slot to `proc.chromR()`.
* Created `peak_to_ploid()` to call peaks and calculate dfe from `freq_peak()` output.
* Created `freq_peak_plot()` to help visualize the output of `freq_peak()`.
* `.vcf_stats_gz` now has nrows and skip parameters.
* removed `.Call()` statements to standardize style.
* Created `vcfR2migrate()` to output MigrateN format data.
* Addressed clang-UBSAN memory leak in `freq_peak()`.
* Created `pairwise_genetic_diff()` to calculate pairwise differentiation.

# vcfR 1.5.0
Released on CRAN 2017-05-18.

* Created `genetic_diff()` to calculate fixation indicies.
* Addressed symbol recognition NOTE: https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661.
* Moved `pinfsc50.png` to tools.
* Added `samples` parameter to vcfR method `[`.
* Deprecated the parameters 'chrom.s' and 'chrom.e' of 'chromo()', please use 'xlim' instead.
* Added `length()` method for chromR objects.
* `[` method throws warning if FORMAT is omitted.
* `plot()` for signature 'chromR' handles INFO column when its all NA.
* `create.chrom()` subsets to first chromosome when more than one is provided.
* adegenet::nLoc(NULL) appears to generate an error when converting data types.


# vcfR 1.4.0
Released on CRAN 2017-01-07.

* `masplit()` converts '.' to NA.
* `extract.indels()` does not recognize NA as a deletion.
* Added parameter `getINFO` to `getFIX()` to suppress INFO column.
* Prof Brian Ripley brought to my attention that I have new memory access issues:
The memory-access errors are new this version, and there is also undefined behavour (trying to coerce NaN to integer).


# vcfR 1.3.0
Released on CRAN 2016-12-08.

* `extract.gt()` no longer uses parameter `allele.sep()`. 
* Added more info to chromR show method.
* When annotation data include more than one chromosome in `create.chromR()` the data are subset to the first chromosome. Thank you Christian!
* added `convertNA` parameter to `extract.gt()` to allow preservation of VCF encoding of missing data. Thank you Thierry!
* added `convertNA` parameter to `read.vcfR()` to allow preservation of VCF encoding of missing data. Thank you Thierry!
* extract.haps omits gt.split and implements unphased_as_NA
* gtsplit handles a mixture of phased and unphased data
* Added 'getters' for vcfR and chromR slots. Thanks Zhian!
* Created `freq_peak()` to find peaks in allele balance frequency data.
* Created `masplit()` to parse matrices contains delimited strings.
* Created `ordisample()` to ordinate sample information.
* `extract.gt()` can now use the ID column from the fix region for rownames.
* Created `INFO2df()` and `metaINFO2df()`.
* Prof Brian Ripley made me aware of memory leaks reported by valgrind.
  Conditional jump or move depends on uninitialised value(s) - write_vcf_body file initialization issue resolved.

# vcfR 1.2.0
Released on CRAN 2016-07-25.

* `vcfR2genind()` greps genotypes containing a missing allele ('.') and sets to NA.
* dplyr v0.5.0 broke some vcfR2tidy functionality. This functionality should be fixed in this release.
* `is_het()` rapidly identifies heterozygotes.
* `extract.info()` scores missing elements as NA.


# vcfR 1.1.0
Released on CRAN 2016-05-26.

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

