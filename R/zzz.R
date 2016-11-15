.onAttach <- function(libname, pkgname){
#  pkg.version <- packageDescription("vcfR", fields = "Version")
  pkg.version <- utils::packageVersion("vcfR")
  
  startup.txt <- paste("\n",
                       "   >>> This is vcfR ", pkg.version, " <<<\n",
                       "   > To cite: citation('vcfR')\n",
                       "   > Documentation: browseVignettes('vcfR')",
                       "\n",
                       sep="")

  packageStartupMessage(startup.txt)
}

