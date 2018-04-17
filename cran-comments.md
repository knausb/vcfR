
## Resubmission
This is a resubmission.
CRAN reported the following after submission.

Found the following (possibly) invalid URLs:
  URL: https://software.broadinstitute.org/gatk/
    From: inst/doc/converting_data.html
    Status: 500
    Message: Internal Server Error

Also fails manually for me. A temporary issue?

This link does work from my computer.
I have also asked several lab mates to check it and they have all validated that it works.
I have also rerun tests on winbuilder R version 3.4.4.
https://win-builder.r-project.org/N5kC8mk7aG29/
and winbuilder R version 3.5.0 RC (2018-04-16 r74611)
https://win-builder.r-project.org/Jh34i9ic2m07/
and have failed to reproduce the behaviour.
I think this indicates that the reported behaviour was temporary.

Incidentally, this is an updated link in this release (1.8.0) versus the previous (1.7.0).
The reason for this update is that R CMD check --as-cran told me I had an old link and it even gave me the updated link that is in this current version.
So I think its doing its job.

## Test environments
* local: ubuntu 16.04 LTS, R 3.4.4
* local: OS X install, R 3.4.4
* ubuntu 14.04.5 LTS (on travis-ci), R 3.4.4
* rocker/r-devel: R Under development (unstable) (2018-03-16 r74418) with valgrind
* winbuilder: R version 3.4.4 (2018-03-15)
* winbuilder: R version 3.5.0 RC (2018-04-15 r74605)

* Windows Server 2012 R2 x64 (build 9600; on AppVeyor), R version 3.5.0 RC (2018-04-15 r74605)
Error : package 'Rcpp' was installed by an R version with different internals; it needs to be reinstalled for use with this R version
I do not believe that this ERROR is due to vcfR.

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Brian J. Knaus <briank.lists@gmail.com>'

Version contains large components (1.7.0.9000)


Possibly mis-spelled words in DESCRIPTION:
  DNAbin (9:46)
  VCF (2:33, 3:68, 4:62, 5:5, 8:30, 10:5)
  VcfR (9:55)
  genlight (9:36)

I have reviewed these words and feel they are spelled correctly.
'DNAbin' refers to an object of class ape::DNAbin.
'VCF' refers to the variant call format specification, a format of file handled by this package.
'VcfR' refers to this package.
'genlight' refers to an object of class adegenet::genlight.


https://cran.r-project.org/web/checks/check_results_vcfR.html
Additional issues:

I believe I have addressed the WARNings and valgrind issue detected at CRAN.


## Downstream dependencies

I have also run R CMD check on downstream dependencies of vcfR
All packages that I could install passed:

devtools::revdep_check()

With one exception:
annovar
https://CRAN.R-project.org/package=annovarR 

Checked annovarR: 1 error  | 0 warnings | 0 notes

I feel that this is not due to vcfR.


## Thank you CRAN Core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

v1.7.0 Thank you Uwe Ligges for processing my submission!

v1.6.0 Thank you Uwe Ligges for processing my submission!

v1.5.0 Thank you Kurt Hornik for helping me!

v1.4.0 Thank you Uwe Ligges and Kurt Hornik for helping me!

v1.3.0 Thank you Uwe Ligges and Kurt Hornik for helping me!
Prof Brian Ripley brought to my attention that I had new memory access issues.
Thank you again for bringing these to my attention!

v1.2.0 Thank you Uwe Ligges for helping me!
This version was accepted on the first try.
See Uwe, I'm learning!

Prof. Brian Ripley also brought to my attention that I had overlooked memory-access errors and that valgrind was reporting use of uninitialized memory and many small memory leaks.
Thank you for bringing this to my attention, its a big help for those of us who are still learning valgrind!

v1.1.0 Thank you Uwe Ligges for helping me get my title in title case, my Description in order and handling my submission!

