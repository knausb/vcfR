.onAttach <- function(libname, pkgname){
#  pkg.version <- packageDescription("vcfR", fields = "Version")
  pkg.version <- utils::packageVersion("vcfR")
  
  startup.txt <- paste("\n",
                       "   *****       ***   vcfR   ***       *****\n",
#                       "   *****       *****      *****       *****\n",
                       "   This is vcfR ", pkg.version, " \n",
                       "     browseVignettes('vcfR') # Documentation\n",
                       "     citation('vcfR') # Citation\n",
#                       "   > To cite: citation('vcfR')\n",
#                       "   > Documentation: browseVignettes('vcfR')\n",
                       "   *****       *****      *****       *****",
                       "\n",
                       sep="")

  packageStartupMessage(startup.txt)
}

