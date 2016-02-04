## Test environments
* local ubuntu 12.04, R 3.2.3
* ubuntu 12.04 (on travis-ci), R 3.2.3
* local OS X install, R 3.2.3
* win-builder (devel and release)
* ubuntu 12.04, R Under development (unstable) (2016-01-11 r69918) (rocker r-devel container)


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Brian Knaus <briank.lists@gmail.com>’
  New submission

* checking package dependencies ... NOTE
  No repository set, so cyclic dependency check skipped
  
Only occurs on travis-ci.
Link of relevance = https://github.com/travis-ci/travis-ci/issues/4125.

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    doc    2.8Mb
    libs   4.3Mb

There are 8 vignettes (HTML) which contribute to the 2.8Mb in doc.
These could be migrated to a separate documentation package.
This could reduce this package's size to 4.9Mb.
Because this is a NOTE and not a ERROR or WARNING I have left it as one package for now.


Possibly mis-spelled words in DESCRIPTION:
  fasta (3:92)
  gff (3:84)
  vcf (2:52, 3:58, 3:79, 3:167)
  
The words fasta, gff and vcf are data file formats and are correctly spelled.


## Downstream dependencies


