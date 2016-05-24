## Test environments
* local ubuntu 12.04, R 3.2.5
* ubuntu 12.04 (on travis-ci), R 3.3.0
* local OS X install, R 3.3.0
* win-builder (devel and release)


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Brian J. Knaus <briank.lists@gmail.com>’

* checking installed package size ... NOTE
  installed size is  7.7Mb
  sub-directories of 1Mb or more:
    doc    2.8Mb
    libs   4.3Mb

There are 8 vignettes (HTML) which contribute to the 2.8Mb in doc.
These could be migrated to a separate documentation package.
This could reduce this package's size.
Because this is a NOTE and not an ERROR or WARNING I have left it as one package for now.


## Downstream dependencies


