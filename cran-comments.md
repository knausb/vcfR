## Test environments
* local ubuntu 12.04, R 3.2.3
* ubuntu 12.04 (on travis-ci), R 3.2.3
* local OS X install, R 3.2.3
* win-builder (devel and release)


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Brian Knaus <briank.lists@gmail.com>’
  New submission

* checking package dependencies ... NOTE
  No repository set, so cyclic dependency check skipped
  
Occurs on travis-ci because no CRAN repository is set (e.g., in ~/.Rprofile)
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
  Vcf (2:52)
  fasta (3:95)
  gff (3:86)
  vcf (3:58, 3:80, 3:171)
  
The words fasta, gff and vcf are data file formats and are correctly spelled.


## Downstream dependencies


