
## Resubmission
This is a resubmission. In this version I have:

Changed:
http://cran.r-project.org/package=pinfsc50
to:
https://cran.r-project.org/package=pinfsc50
in vcfR.R to respect the CRAN URL cannonical form.

The files in 'inst/doc' have been deleted and rebuilt.


## Test environments
* local ubuntu 16.04 LTS, R 3.3.2
* ubuntu 12.04 (on travis-ci), R 3.3.1
* local OS X install, R 3.3.2
* rhub::check( platform = "ubuntu-gcc-devel" )
* rhub::check( platform = "windows-x86_64-release" )
* rhub::check( platform = "windows-x86_64-devel" )


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Brian J. Knaus <briank.lists@gmail.com>’

* checking installed package size ... NOTE
  installed size is  8.6Mb
  sub-directories of 1Mb or more:
    doc    2.2Mb
    libs   5.7Mb

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

I have reviewed these words and feel they are spelled correctly.
'DNAbin' referr to an object of class ape::DNAbin.
'VCF' refers to the variant call format specification, a format of file handled by this package.
'VcfR' refers to this package.
'genlight' refers to an object of class adegenet::genlight.
'genomic' refers to properties of a genome.


## Downstream dependencies

I have also run R CMD check on downstream dependencies of vcfR
All packages that I could install passed except:

* pcadapt: 
* checking installed package size ... NOTE
  installed size is  5.1Mb
Does not appear to be related to vcfR


## Memory-access errors

During the last submission Prof. Brian Ripley brought to my attention that vcfR contained memory access errors.
I believe I have addressed this by running commands such as:

R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_conversion.R

On various suspect functions.


## Thank you CRAN core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

v1.3.0 Thank you Uwe Ligges and Kurt Hornik for helping me!

v1.2.0 Thank you Uwe Ligges for helping me!
This version was accepted on the first try.
See Uwe, I'm learning!

Prof. Brian Ripley also brought to my attention that I had overlooked memory-access errors and that valgrind was reporting use of uninitialized memory and many small memory leaks.
Thank you for bringing this to my attention, its a big help for those of us who are still learning valgrind!

v1.1.0 Thank you Uwe Ligges for helping me get my title in title case, my Description in order and handling my submission!

