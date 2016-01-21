library("testthat")
#library("vcfR")
#
#test_check("vcfR")


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
#test_file("testthat/test_chromo_plot.R")

test_file("testthat/test_vcfR_methods.R")


#test_file("testthat/test_rank_variants.R")
#test_file("testthat/test_var_window.R")
#test_file("testthat/test_windowing.R")



##### ##### ##### ##### #####
# Debugging compiled code.

# lldg gdb
# http://kevinushey.github.io/blog/2015/04/13/debugging-with-lldb/

# R -d gdb
# run
# q()
# quit

# Also control-c returns to gdb.
# continue retruns from gdb to R (plus a second return)

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


##### ##### ##### ##### #####
# EOF.