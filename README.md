
## Manipulate and visualize variant call format (vcf) data

VcfR: tools to read, write, parse and visualize [vcf](https://github.com/samtools/hts-specs) data.

[![Travis-CI Build Status](https://travis-ci.org/knausb/vcfR.png?branch=master)](https://travis-ci.org/knausb/vcfR)
[![Coverage Status](https://coveralls.io/repos/github/knausb/vcfR/badge.svg?branch=master)](https://coveralls.io/github/knausb/vcfR?branch=master)


vcfR is an R package intended to allow easy manipulation and visualization of variant call format (vcf) data.
Functions are provided to rapidly read from and write to vcf files.
Once data are read into memory they can be stored in either of two data structures.

*vcfR* - S4 class to contain a vcf file as well as functions to load this object.

*chromR* - S4 class to contain variant information (vcf) as well as sequence (fasta) and annotation (gff) information.

Additional functions provide tha ability to subset vcf data as well as to extract and parse subsets of the data.
For example, individual genotypes, sequence depths or genotype likelihoods (when provided in the file) can easily be accessed.
These tools are provided to aid researchers in rapidly surveying the quality and other characteristics of data provided as vcf data.
With this information in hand, researchers should be able to determine criteria for hard filtering in order to attempt to maximize biological variation and minimize technical variation.




## Download

While this project is in development it can be installed through github:

    devtools::install_github(repo="knausb/vcfR")
    library(vcfR)


If you would like the vignettes use:

    devtools::install_github(repo="knausb/vcfR", build_vignettes=TRUE)


If you've built the vignettes, you can browse them with:

    browseVignettes(package="vcfR")


If you've installed this package with devtools you will probably need to run:

    devtools::install(build_vignettes = TRUE)
    

------

## Development version

The development version (which may not be stable) can also be installed:

    devtools::install_github(repo="knausb/vcfR@devel")
    library(vcfR)


And to build the vignettes:

    devtools::install_github(repo="knausb/vcfR@devel", build_vignettes=TRUE)


------


Enjoy!
