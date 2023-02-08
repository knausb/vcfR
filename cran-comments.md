
## Submission

This package, vcfR, is generating WARNings on debian/fedora with clang. I believe these are addressed here. This submission also includes updates to the package.

## Test environments

* local:
R 4.2.2; Ubuntu 22.04.1 LTS; gcc (Ubuntu 11.3.0-1ubuntu1~22.04) 11.3.0
ubuntu 22.04.1 LTS and R Under development (unstable) (2022-07-07 r82559) clang++-14 -Wall

* local:
OS X Monterey version 12.4 and R 4.2.1 and clang 13.1.6

* GitHub Actions
MacOS-latest (release) - macOS 12.6.2 21G320; R-4.2.2
Ubuntu-latest (devel) - Ubuntu 22.04.1 LTS; R Under development (unstable) (2023-02-01 r83747)
Ubuntu-latest (oldrel-1) - Ubuntu 22.04.1 LTS; R version 4.1.3 (2022-03-10)
Ubuntu-latest (release) - Ubuntu 22.04.1 LTS; R version 4.2.2 (2022-10-31)
Windows-latest (release) - Microsoft Windows Server 2022 10.0.20348 Datacenter; R version 4.2.2 (2022-10-31 ucrt)

* AppVeyor:
Windows Server 2012 R2 x64 (build 9600); R Under development (unstable) (2023-02-01 r83747 ucrt)
Windows Server 2012 R2 x64 (build 9600); R version 4.2.2 (2022-10-31 ucrt)

* win-builder:
None for this submission

* rhub:
None for this submission


## R CMD check results

* checking installed package size ... NOTE
  installed size is  6.5Mb
  sub-directories of 1Mb or more:
    libs   4.1Mb


## revdepcheck results

We checked 18 reverse dependencies (4 from CRAN + 14 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


## Thank you CRAN Core Team!

[CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) states that all correspondence should be with CRAN and not members of the team.
However, I think its polite to thank those who have helped this project.
So I've decided to start a list of thanks with the hope that these individuals may see this in the future.

v1.13.0 Thank you Uwe Ligges for processing my submission!

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

