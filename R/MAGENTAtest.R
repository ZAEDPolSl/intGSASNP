#' MAGENTA test implementation
#'
#' Performs a MAGENTA test on SNP p-values
#'
#' @param data input dataframe with column names: "snp", "entrez", "p.value"
#' @param pathways object containing a list of GENES and a list of MODULES
#' @param int.method integration method (either "min", "fisher', or "stouffer")
#' @param adjust.int flag for multiple testing correction in integration
#' @param int.cor.LD flag for LD correction in integration
#' @param cutoff enrichment p-value cutoff 
#' @param adjust.method multiple testing correction method
#' @param permutations number of gene set permutations
#'
#' @return A data frame with module names, calculated p-value, and additional statistics. 
#'
#' @examples
#' MAGENTAtest(data = example_dataset, pathways = pathway_library, int.method = "min", 
#' adjust.int = FALSE, int.cor.LD = FALSE, cutoff = 0.05, adjust.method = "BH", permutations = 1000)
#'
#' @export
MAGENTAtest <-
  function(data,
           pathways = NA,
           int.method = "min",
           adjust.int = FALSE,
           int.cor.LD = FALSE,
           cutoff = 0.05,
           adjust.method = "BH",
           permutations = 1000) {
    
    validate_input(data)
    validate_pathways(pathways)
    validate_package("fastmatch")
    
    M2G = pathways$MODULES2GENES
    no_of_modules <- length(M2G)
    MOD = pathways$MODULES
    
    cat("Performing p-value integration... \n")
    data <-
      integrate(data,
                method = int.method,
                adjust = adjust.int,
                cor.LD = int.cor.LD)
    N <- nrow(data)
    
    cat("MAGENTA test... \n")
    result <- data.frame(matrix(NA, nrow = 0, ncol = 8))
    colnames(result) <-
      c("ID",
        "Title",
        "GS.size",
        "LEF",
        "LEFa",
        "p",
        "p.adj",
        "p.cutoff")
    p_cutoff <- quantile(data[, 2], cutoff)
    percs <- floor(c(1:10) * 0.1 * no_of_modules) # logging output
    for (i in 1:no_of_modules) {
      module = names(M2G)[i]
      result[i, 1] <- module
      result[i, 2] <- MOD$Title[MOD$ID == module]
      genes <- M2G[[i]]
      GS <- data$entrez[data$entrez %fin% genes]
      GS.size <- length(GS)
      result[i, 3] <- GS.size
      GS.p <- data[data$entrez %fin% GS,]
      LEF <- sum(GS.p[, 2] < p_cutoff)
      result[i, 4] <- LEF
      LEFas <- c()
      for (j in 1:permutations) {
        GS.art <- sample(data[, 1], GS.size)
        GS.art.p <- data[data$entrez %fin% GS.art,]
        LEFas <- c(LEFas, sum(GS.art.p[, 2] < p_cutoff))
      }
      LEFa <- sum(LEFas >= LEF)
      p <- LEFa / permutations
      result[i, 5] <- LEFa
      result[i, 6] <- p
      
      if (i %fin% percs) {
        perc <- which(percs == i) * 10
        cat(paste(perc, "% ... \n", sep = ""))
      }
    }
    result$p.cutoff <- p_cutoff
    result$p.adj <- p.adjust(p = result$p, method = adjust.method)
    cat("MAGENTA test completed. \n")
    return(result)
  }
