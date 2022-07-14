
## Submission

This package, vcfR, is generating WARNings on debian/fedora with clang. I believe these are addressed here. This submission also includes updates to teh package.

## Test environments

* local:
ubuntu 22.04 LTS and R 4.2.1
ubuntu 22.04 LTS and R Under development (unstable) (2022-07-07 r82559)

* local:
OS X Catalina 10.15.6 and R 4.0.2 and clang

win-builder:
* using R version 4.0.2 (2020-06-22)
* using R Under development (unstable) (2020-08-23 r79071)

travis-ci:
* Ubuntu 16.04.6 LTS, R version 4.0.0 (2020-04-24)
* Ubuntu 16.04.6 LTS, R Under development (unstable) (2020-08-26 r79084)

Currently failing:
Error in dyn.load(file, DLLpath = DLLpath, ...) : 
  unable to load shared object '/home/travis/R/Library/ape/libs/ape.so':
  libRlapack.so: cannot open shared object file: No such file or directory

My interpretation is that this is not due to vcfR.

AppVeyor:
* Windows Server 2012 R2 x64 (build 9600), R Under development (unstable) (2020-08-24 r79088)
* Windows Server 2012 R2 x64 (build 9600), R version 4.0.2 (2020-06-22)

rhub:
* None for this submission


## R CMD check results


* checking installed package size ... NOTE
  installed size is  6.5Mb
  sub-directories of 1Mb or more:
    libs   4.1Mb








# Delete checks below.


Package was archived on CRAN

Possibly mis-spelled words in DESCRIPTION:
  DNAbin (9:46)
  VCF (2:33, 3:68, 4:62, 5:5, 8:30, 10:5)
  VcfR (9:55)
  genlight (9:36)

These words and acronyms are esoteric to working with genomic data and are all correctly spelled.

Found the following (possibly) invalid URLs:
  URL: https://uswest.ensembl.org/info/docs/tools/vep/index.html
    From: man/vep.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)
  URL: https://www.internationalgenome.org/node/101
    From: inst/doc/intro_to_vcfR.html
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)

These URLs all work when pasted into firefox.


* checking installed package size ... NOTE
  installed size is 10.4Mb
  sub-directories of 1Mb or more:
    libs   8.5Mb

This has not been an issue in the past.


* checking for future file timestamps ... NOTE
unable to verify current time

I interpret this as not an issue with vcfR.


## Thank you CRAN Core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

v1.12.0 Thank you Uwe Ligges for processing my submission!

v1.11.0 Thank you Uwe Ligges for processing my submission!

v1.10.0 Thank you Uwe Ligges for processing my submission!

v1.9.0 Thank you Uwe Ligges for processing my submission!

v1.8.0 Thank you Uwe Ligges for processing my submission!

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

