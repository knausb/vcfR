.onUnload <- function (libpath) {
  library.dynam.unload("vcfR", libpath)
}
