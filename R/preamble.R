.onAttach <- function(...) {
  # Begin Exclude Linting
  greet <- paste("# This research was partially supported under NSF grants 1463642 and 1915277")
  packageStartupMessage(greet)
  # End Exclude Linting
}

.onUnload <- function(libpath) {
  library.dynam.unload("smerc", libpath)
}
