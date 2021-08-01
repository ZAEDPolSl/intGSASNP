validate_pathways <- function(pathways) {
  if (is.null(pathways) ||
      is.na(pathways) ||
      !is.list(pathways$MODULES) ||
      !is.list(pathways$GENES) ||
      !is.list(pathways$MODULES2GENES) ||
      !is.list(pathways$GENES2MODULES))
    stop("Pathways object should contain a list of GENES and a list of MODULES.")
}

validate_input <- function(input) {
  if (is.null(input) ||
      is.na(input) ||
      colnames(input) != c("snp", "entrez", "p.value"))
    stop("Input object should be a dataframe with column names: snp, entrez, p.value.")
}

validate_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dep = TRUE)
    if (!require(pkg, character.only = TRUE))
      stop(paste("Package", pkg, "not found.", sep = " "))
  }
}
