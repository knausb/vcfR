
## Explore and manipulate variant call format (vcf) files

VcfR, tools created to work with [vcf](https://github.com/samtools/hts-specs) files.

[![Travis-CI Build Status](https://travis-ci.org/knausb/vcfR.png?branch=master)](https://travis-ci.org/knausb/vcfR)

*vcfR* - S4 class to contain a vcf file as well as functions to load this object.

*Chrom* - S4 class to cantain variant information (vcf) as well as sequence (fasta) and annotation (gff) information.

Provides functions to load data, subset, filter and visualize these data.


While this project is in development it can be installed through github:

    library(devtools)
    install_github(repo="knausb/vcfR")
    library(vcfR)


The development version (which may not be stable) can also be installed:

    devtools::install_github(repo="knausb/vcfR@devel")
    library(vcfR)


If you would like the vignettes use:

    devtools::install_github(repo="knausb/vcfR", build_vignettes=TRUE)
    
or:

    devtools::install_github(repo="knausb/vcfR@devel", build_vignettes=TRUE)

If you've built the vignettes you can browse them with:

    browseVignettes(package="vcfR")
    
which should open a page in your web browser.

If you've installed this package with devtools you will probably need to run:

    devtools::install(build_vignettes = TRUE)

In order to build the vignettes (which does not happen by default).

Enjoy!
