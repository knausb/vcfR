library("testthat")
#library("vcfR")
#
#test_check("vcfR")


#
test_file("testthat/test_10_write_vcf.R")
#
test_file("testthat/test_1_vcf.R")
#
test_file("testthat/test_2_chromR.R")
#
test_file("testthat/test_3_extract_gt.R")
#
test_file("testthat/test_4_vcfR2DNAbin.R")
#
test_file("testthat/test_io.R")
#
test_file("testthat/test_conversion.R")
#
test_file("testthat/test_vcfR_methods.R")
#
test_file("testthat/test_addID.R")
#
test_file("testthat/test_maf.R")
#
test_file("testthat/test_genotype_matrix_functions.R")
#
test_file("testthat/test_vcfRtidy.R")
#
test_file("testthat/test_heatmapbp.R")
#
test_file("testthat/test_ad_frequency.R")
#
test_file("testthat/test_chromR_method.R")
#
test_file("testthat/test_get.R")
#
test_file("testthat/test_freq_peak.R")


# This test will write the plots to a file on the filesystem.
# This should cause CRAN to complain.
# Comment out this test prior to CRAN submission.
#
#
test_file("testthat/test_chromo_plot.R")
#
test_file("testthat/test_drplot.R")


##### ##### ##### ##### #####


#test_file("testthat/test_rank_variants.R")
#test_file("testthat/test_var_window.R")
#test_file("testthat/test_windowing.R")



##### ##### ##### ##### #####
# Notes on debugging compiled code.

# lldg gdb
# http://kevinushey.github.io/blog/2015/04/13/debugging-with-lldb/

# R -d gdb
# run
# q()
# quit

# Also control-c returns to gdb.
# continue returns from gdb to R (plus a second return)

# in gdb
# b rcppMean
# sets a breakpoint at the function rcppMean

# b read_body_gz
# c

# breakpoint set --line 46



# valgrind
# http://kevinushey.github.io/blog/2015/04/05/debugging-with-valgrind/

# R -d "valgrind --db-attach=yes" -f test_1_vcf.R
# R -d "valgrind --track-origins=yes" -f test_1_vcf.R

# R -d "valgrind --track-origins=yes" -f testthat/test_conversion.R
# R -d "valgrind --track-origins=yes" -f testthat/test_windowing.R


# 2016-08-05

# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_1_vcf.R
# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_10_write_vcf.R
# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_freq_peak.R

# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_3_extract_gt.R


# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_conversion.R 
# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat/test_conversion.R 

# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat.R

# R -d "valgrind --leak-check=full --vgdb-error=1" -f tests/testthat.R --restore --save --no-readline --vanilla > log.txt 2>&1


##### ##### ##### ##### #####
# g++

# g++ -std=c++0x -I"C:/PROGRA~1/R/R-33~1.2/include" -DNDEBUG -I"C:/Users/USERtqITJhLcoQ/R/Rcpp/include" -I"d:/Compiler/gcc-4.9.3/local330/include" -O2 -Wall -mtune=core2 -c gt_to_popsum.cpp -o gt_to_popsum.o
# g++ -std=c++11 -I/usr/share/R/include -DNDEBUG   -I"/home/knausb/R/x86_64-pc-linux-gnu-library/3.3/Rcpp/include"   -fpic  -g -O2 -fstack-protector-strong -Wsign-compare -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g -c gt_to_popsum.cpp -o gt_to_popsum.o


##### ##### ##### ##### #####
# EOF.