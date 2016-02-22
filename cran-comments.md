## Resubmission
This is a second resubmission. In this version I have addressed the below comments:

Thanks, processing the vignettes takes a while:
* checking re-building of vignette outputs ... [334s] OK
Can this be reduced?

This package is intended to process genomic data.
These datasets can be rather large.
There are 8 vignettes in this package (this also contributeds to its size in NOTE below).
I am of the opinion that one of the goals of these vignettes is that they should demonstrate the code performs reasonably on modest size datasets.
(This opinion appears to be causing me be a bit of trouble.)
I've created a data only package that is just under 5 MB that is used in these vignettes:

https://cran.r-project.org/web/packages/pinfsc50/index.html

This data set includes 18 samples and 22,031 variants.
This dataset is processed in the vignettes and is why they take long to build.
I have tried to reduce the build time of the vignettes by:

- I have set eval=FALSE on code chuncks that demonstrate translation to data structures supported by other CRAN packages, some of which do not perform well on this size dataset (i.e., the slow down is due to code from these other packages).
- I have tried to remove redundancy among the vignettes which has reduced their size.

This has reduced the vignette build time to 128s (at win-builder).
If this is not an adequate reduction perhaps I should omit the vignettes from the package and host them elsewhere?
Feed back from CRAN would be appreciated on this.

Thank you and apologies for a wordy response!



## Resubmission
This is a resubmission. In this version I have addressed the below comments:

> Description: Tools for working with variant call format (vcf) files. Reads
>   in "vcf", gff and "fasta" data for visualization.  Includes

* Can you pls use single quotes?
Double quotes in the DESCRIPTION have been changed to single quotes.

* Can you pls provide a refernce/URL for "vcf format"?
http://samtools.github.io/hts-specs/

* Also, is it vcf?  Vcf (as in the title)?  VCF?
Review of the definition indicates that it should be VCF (my error). This has been updated.


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
  VCF (2:53, 3:59)
  
The word VCF is a data file format and is correctly spelled.


## Downstream dependencies


