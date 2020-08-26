
This release includes updates to address requests from CRAN to incorporate checkbashisms on Debian flavors or Linux.

## Test environments

* local:
ubuntu 18.04 LTS and R 4.0.2

* local:
OS X Catalina 10.15.6 and R 4.0.2 and clang

win-builder:
* using R version 4.0.2 (2020-06-22)
* using R Under development (unstable) (2020-07-19 r78884)

travis-ci:
* ubuntu 16.04 LTS, R version 4.0.0 (2020-04-24)
* ubuntu 16.04 LTS, R Under development (unstable) (2020-07-23 r78905)

AppVeyor:
* Windows Server 2012 R2 x64 (build 9600) R version 3.6.2 Patched (2020-01-25 r77764)

rhub:



Issues:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Brian J. Knaus <briank.lists@gmail.com>’

New submission

Package was archived on CRAN

Found the following (possibly) invalid URLs:
  URL: https://gatkforums.broadinstitute.org/gatk/discussion/6926/spanning-or-overlapping-deletions-allele
    From: man/vcfR2DNAbin.Rd
    Status: Error
    Message: libcurl error code 6:
      	Could not resolve host: gatkforums.broadinstitute.org
        
        
        


## R CMD check results
There were no ERRORs.


There were 2 NOTEs:

Found the following (possibly) invalid URLs:
  URL: http://www.1000genomes.org/node/101
    From: inst/doc/intro_to_vcfR.html
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)
  URL: https://uswest.ensembl.org/info/docs/tools/vep/index.html
    From: man/vep.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)
      	
This url works when I copy and paste it into firefox.


* checking installed package size ... NOTE
  installed size is 10.4Mb
  sub-directories of 1Mb or more:
    libs   8.5Mb




## Thank you CRAN Core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

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

