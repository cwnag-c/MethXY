.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  bioc_pkgs <- c(
    "SummarizedExperiment",
    "S4Vectors",
    "minfi",
    "limma",
    "genefilter",
    "FlowSorted.Blood.EPIC"
  )
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}
