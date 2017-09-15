#
#
#library("testthat")
context("chromo_plot  functions")

library(vcfR)



##### ##### ##### ##### #####
# chromo, vcf only

test_that("chromo works, variant data only",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
  chrom <- proc.chromR(chrom, verbose = FALSE)

  expect_is(chrom, "chromR")
  
  plot1 <- chromo( chrom, boxp = FALSE )
  expect_true( is.null(plot1) )
  
  plot2 <- chromo( chrom, boxp = TRUE )
  expect_true( is.null(plot2) )
})


test_that("chromR plot when DP == NA",{
#  library(vcfR)
#  trace('plot', browser, exit=browser, signature='chromR')
  data(vcfR_test)
  is.na(vcfR_test@fix[,'INFO']) <- TRUE
  chrom <- create.chromR(vcfR_test, verbose = FALSE)
  myPlot <- plot(chrom)
  expect_true( is.null(myPlot) )
})


test_that("chromR plot when MQ == NA",{
#  trace('plot', browser, exit=browser, signature='chromR')
  data(vcfR_test)
  is.na(vcfR_test@fix[,'INFO']) <- TRUE
  chrom <- create.chromR(vcfR_test, verbose = FALSE)
  myPlot <- plot(chrom)
  expect_true( is.null(myPlot) )
})



##### ##### ##### ##### #####
# chromo, vcf and seq


test_that("chromo works, variant and seq data",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
#  chromo( chrom ) # Should error!
  chrom <- proc.chromR(chrom, verbose = FALSE)

  expect_is(chrom, "chromR")
    
  plot1 <- chromo( chrom, boxp = FALSE )
  expect_true( is.null(plot1) )
  plot2 <- chromo( chrom, boxp = TRUE )
  expect_true( is.null(plot2) )
})


##### ##### ##### ##### #####
# chromo, vcf and annotation

test_that("chromo works, variant and annotation data",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig", vcf=vcf, ann=gff, verbose=FALSE)
#  chromo( chrom )
  chrom <- proc.chromR(chrom, verbose = FALSE)
  
  expect_is(chrom, "chromR")
  
  plot1 <- chromo( chrom, boxp = FALSE )
  expect_true( is.null(plot1) )
  plot2 <- chromo( chrom, boxp = TRUE )
  expect_true( is.null(plot2) )
  
})


##### ##### ##### ##### #####
# chromo, vcf, seq and annotation

test_that("chromo works, variant and annotation data",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
#  chromo( chrom )
  chrom <- proc.chromR(chrom, verbose = FALSE)

  expect_is(chrom, "chromR")
  
  plot1 <- chromo( chrom, boxp = FALSE )
  expect_true( is.null(plot1) )
  plot2 <- chromo( chrom, boxp = TRUE )
  expect_true( is.null(plot2) )
})


test_that("chromo works, vcf with no variants, seq or ann",{
  data("vcfR_example")
  vcf <- vcf[0,]
  
#  chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
#  chrom@len
  
  chrom <- proc.chromR(chrom, verbose = FALSE)
  
  expect_is(chrom, 'chromR')
  expect_equal( nrow(chrom@win.info), 0 )
  
})



test_that("chromo works, vcf with no variants",{
  data("vcfR_example")
  vcf <- vcf[0,]
  
  chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
# chrom@len
  
  chrom <- proc.chromR(chrom, verbose = FALSE)
  
  expect_is(chrom, 'chromR')
  expect_gt( nrow(chrom@win.info), 0 )
  
})



##### ##### ##### ##### #####
# chromo with custom tracks.


test_that("chromo works, custom tracks",{
  data("vcfR_example")
  chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
  
  ##### ##### ##### ##### #####
  # Create lists of dots and rectangles
  rlst1 <- cbind( gff[ seq(1,23, by=2) , 4 ],
                  0, 
                  gff[ seq(1,23, by=2) , 5 ],
                  500
                )
  rlist <- list(rlst1)
  myList1 <- list(title = "Track1",
                  dmat  = chrom@var.info[,2:4],
                  rlst = rlst1,
                  rcol=4,
                  bwcol=1:2
                  )

  ##### ##### ##### ##### #####
  # Custom plots.
  plot1 <- chromo( chrom, boxp = FALSE, drlist1 = myList1 )
  expect_true( is.null(plot1) )

  chrom <- proc.chromR(chrom, verbose = FALSE)
  plot2 <- chromo( chrom, boxp = FALSE , drlist1 = myList1 )
  expect_true( is.null(plot2) )
  
  plot3 <- chromo( chrom, boxp = TRUE , drlist1 = myList1  )
  expect_true( is.null(plot3) )
  
  plot4 <- chromo( chrom, boxp = TRUE , drlist1 = myList1, drlist2 = myList1  )
  expect_true( is.null(plot4) )

  plot5 <- chromo( chrom, boxp = TRUE , drlist1 = myList1, drlist2 = myList1, drlist3 = myList1  )
  expect_true( is.null(plot5) )
})


# rlist, dcol, rcol, rbcol and bwcol.)
#names(myList1)

##### ##### ##### ##### #####
# Clean up when we're done.

unlink('Rplots.pdf')

# EOF.