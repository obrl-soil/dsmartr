.onAttach <- function(libname, pkgname)  {
  packageStartupMessage("dsmartr ", utils::packageDescription("dsmartr", field="Version"), " attached", appendLF = TRUE)
}
