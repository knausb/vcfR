
# vcfR 2.0.0.
There is currently no plan to release vcfR 2.0.0.
If and when this 'major' release occurs it will include changes that will break backward compatibility.
At the present, this is simply a to-do list for ideas to include in the next major release.

* `gzread` needs us to tell it a maximum number of bytes to read in, so we need to make an arbitrary guess.
I think I encountered a situation where 4-96 was not enough so I've bumped it to 16,384 B.
* Move 'FORMAT' column to its own slot. We can then cbind FORMAT and gt when passing to compiled code.
This may have been addressed at 64a308ba50b9119108e8946737460de5997b805b by adding `samples` to vcfR method `[`.
* In issue #92 (vcfR2genlight big data #92), JimWhiting91 has documented that `extract.gt()` could be greatly improved with multithreading. While he used `mclapply()` I do not feel this is the best solution because it does not work on Windows. I think a better solution would be [RCppParallel](https://rcppcore.github.io/RcppParallel/) because this should work on all CRAN platforms.

# vcfR 1.13.0
Released on CRAN 202X-XX-XX

# vcfR 1.12.0
Released on CRAN 2020-09-01
* Added ```PKGTYPE: both``` to appveyor.yml so Windows packages can be built from source
* Omitted ```configure``` file that unnecessarily tried to invoke checkbashisms
* Incorporated help from https://stackoverflow.com/a/62721142 to use checkbashisms when checking on Debian flavors of Linux (ended up omitting this change but left this here to document it and the link)

# vcfR 1.11.0
Released on CRAN 2020-06-05
* Now compatible with R 4.0.0 and dplyr 1.0.0

# vcfR 1.10.0
Released on CRAN 2020-02-06
* Handled deprecated "dplyr::verb_" function in vcfR2tidy
* Omitted unused elipses from proc.chromR()

# vcfR 1.9.0
Released on CRAN 2020-01-10
* Changed class(x) == "matrix" to inherited(x, "matrix")
* Changed license from `GPL` to `GPL-3` (#144).
* `extract.haps()` reports the correct number of variants processed when verbose.
* The square brackets ([]) handle @gt slots with no samples.
* `vcfR2loci()` now has the option `return.alleles = FALSE`.
* `vcfR2genind()` now has the option `retrun.alleles = FALSE`.
* Error handling code moved into the C++ functions called by read.vcfR so that errors are thrown earlier when reading a VCF. read.vcfR no longer checks that a file is readable first, which solves issues sometimes seen with shared files. (Issue #109, reported and fixed by @NikNakk).
* `extract.haps()` did not include the parameter `return.alleles = TRUE` in it's call to `extract.gt()` in the haploid branch of the function. This parameter has now been added. This also affects `vcfR2DNAbin()` which calls this function.
* `vcfR2genlight()` includes the parameter `...` to pass parameters to `adegenet::df2genind()`.
* `is.indel()` returns logical vector to identify indels.
* gt.to.popsum now handles genotypes that include some, but not all, missing alleles.

# vcfR 1.8.0
Released on CRAN 2018-04-17
* Attempted to address CRAN's 'Note: break used in wrong context: no loop is visible' issue.
* `.vcf_stats_gz()` reports number of elements in header as well as the file's last line. This is used by `read.vcfR()` to check for poorly formed files.
* `show` method for vcfR now queries @fix instead of @gt.
* `check_keys()` checks key definitions in the meta section to make sure they are unique.
* `freq_peak_plot()` has parameter `posUnits` to adjust units of scatterplot.
* `vcfR2migrate()` manual discusses Unix and Windows line endings.


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
* Created `pairwise_genetic_diff()` to calculate pairwise differentiation. Thanks Javier!

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

