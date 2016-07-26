

Release of dplyr 0.5.0 affected functionality of some vcfR code.
This release should address this functionality.

## Test environments
* local ubuntu 12.04, R 3.2.5
* ubuntu 12.04 (on travis-ci), R 3.3.1
* local OS X install, R 3.3.0
* win-builder (devel and release)


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Brian J. Knaus <briank.lists@gmail.com>’

* checking installed package size ... NOTE
  installed size is  7.6Mb
  sub-directories of 1Mb or more:
    doc    2.2Mb
    libs   4.8Mb

There are 8 vignettes (HTML) which contribute to the 2.8Mb in doc.
These could be migrated to a separate documentation package.
This could reduce this package's size.
Because this is a NOTE and not an ERROR or WARNING I have left it as one package for now.


Possibly mis-spelled words in DESCRIPTION:
  DNAbin (9:76)
  VCF (2:33, 3:68, 4:62, 5:10, 8:58, 10:34)
  VcfR (10:5)
  genlight (9:66)
  genomic (7:51)

I have reviewed these word and feel they are spelled correctly.
'DNAbin' referr to an object of class ape::DNAbin.
'VCF' refers to the variant call format specification, a format of file handled by this package.
'VcfR' refers to this package.
'genlight' refers to an object of class adegenet::genlight.
'genomic' refers to properties of a genome.


## Downstream dependencies



## Thank you CRAN core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

v1.2.0 Thank you Uwe Ligges for helping me!
This version was accepted on the first try.
See Uwe, I'm learning!

v1.1.0 Thank you Uwe Ligges for helping me get my title in title case, my Description in order and handling my submission!

