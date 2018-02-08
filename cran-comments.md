

## Test environments
* local: ubuntu 16.04 LTS, R 3.4.3
* local: OS X install, R 3.4.3
* ubuntu 14.04.5 LTS (on travis-ci), R 3.4.2
* Windows Server 2012 R2 x64 (build 9600; on AppVeyor), R version 3.4.3 Patched (2018-02-03 r74215)
* winbuilder: R version 3.4.3 (2017-11-30)
* winbuilder: R Under development (unstable) (2018-02-01 r74194)


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... Note_to_CRAN_maintainers
Maintainer: ‘Brian J. Knaus <briank.lists@gmail.com>’

* checking installed package size ... NOTE
  installed size is  9.5Mb
  sub-directories of 1Mb or more:
    libs   7.8Mb


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


## Downstream dependencies

I have also run R CMD check on downstream dependencies of vcfR
All packages that I could install passed:

devtools::revdep_check()

Checked annovarR: 1 error  | 0 warnings | 0 notes
Checked pcadapt : 0 errors | 0 warnings | 3 notes

The ERROR from annovarR is below, I do not feel this is related to vcfR.

ERROR: dependencies ‘RMySQL’, ‘AnnotationDbi’ are not available for package ‘annovarR’
* removing ‘/tmp/Rtmp4bz1PA/R-lib/annovarR’

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

