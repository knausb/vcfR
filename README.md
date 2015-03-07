
## Explore and manipulate variant call format files

VcfR, tools created to work with vcf files.

[![Travis-CI Build Status](https://travis-ci.org/knausb/vcfR.png?branch=master)](https://travis-ci.org/knausb/vcfR)

*vcfR* - S4 class to contain a vcf file as well as functions to load this object.

*Chrom* - S4 class to cantain a vcf file as well as sequence (fatsta) and annotation (gff) information.

Provides functions to load data, subset, filter and visualize these data.


While this project is in development it can be installed through github:

    library(devtools)
    install_github(repo="knausb/vcfR")
    library(vcfR)


The development version (which may not be stable) can also be installed:

    devtools::install_github(repo="knausb/vcfR@devel")
    library(vcfR)


If you would like the vignettes use:

    devtools::install_github(repo="knausb/vcfR@devel", build_vignettes=TRUE)

Enjoy!
