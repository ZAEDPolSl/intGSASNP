#' ORA test implementation
#'
#' Performs a ORA test on SNP p-values
#'
#' @param data input dataframe with column names: "snp", "entrez", "p.value"
#' @param pathways object containing a list of GENES and a list of MODULES
#' @param int.method integration method (either "min", "fisher', or "stouffer")
#' @param adjust.int flag for multiple testing correction in integration
#' @param int.cor.LD flag for LD correction in integration
#' @param adjust.method multiple testing correction method
#'
#' @return A data frame with module names, calculated p-value, and additional statistics. 
#'
#' @examples
#' ORAtest(data = example_dataset, pathways = pathway_library, int.method = "min", 
#' adjust.int = FALSE, int.cor.LD = TRUE, adjust.method = "BH")
#'
#' @export
ORAtest <-
  function(data,
           pathways = NA,
           int.method = "min",
           adjust.int = FALSE,
           int.cor.LD = FALSE,
           adjust.method = "BH") {
    
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
    
    cat("ORA test... \n")
    result <- data.frame(matrix(NA, nrow = 0, ncol = 5))
    colnames(result) <- c("ID", "Title", "p", "p.adj", "OR")
    percs <- floor(c(1:10) * 0.1 * no_of_modules) # logging output
    for (i in 1:no_of_modules) {
      module = names(M2G)[i]
      result[i, 1] <- module
      result[i, 2] <- MOD$Title[MOD$ID == module]
      genes <- M2G[[i]]
      significant <- data[, 2] < 0.05
      K <- sum(significant)
      M <- length(data$entrez[data$entrez %fin% genes])
      significant_genes <- data$entrez[significant]
      x <- length(significant_genes[significant_genes %fin% genes])
      tab <- cbind(c(x, K - x), c(M - x, N - M - K + x))
      p <- min(1 - cumsum(dhyper(0:(x - 1), M, N - M, K)))
      OR <- (x * (N - M - K + x)) / ((M - x) * (K - x))
      result[i, 3] <- p
      result[i, 5] <- OR
      
      if (i %fin% percs) {
        perc <- which(percs == i) * 10
        cat(paste(perc, "% ... \n", sep = ""))
      }
    }
    result$p.adj <- p.adjust(p = result$p, method = adjust.method)
    cat("ORA test completed. \n")
    return(result)
  }
