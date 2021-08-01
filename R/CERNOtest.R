#' CERNO test implementation
#'
#' Performs a CERNO test on SNP p-values
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
#' CERNOtest(data = example_dataset, pathways = pathway_library, int.method = "min", 
#' adjust.int = FALSE, int.cor.LD = FALSE, adjust.method = "BH")
#'
#' @export
CERNOtest <-
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
    
    cat("CERNO test... \n")
    result <- data.frame(matrix(NA, nrow = 0, ncol = 8))
    colnames(result) <- c("ID", "Title", "N", "F-test", "ES", "p", "p.adj", "AUC")
    x <- data[order(data[, 2]),]
    ranks <- c(1:N) / N
    entrez = x[, 1]
    
    cP <- function(cF) {
      len <- length(cF)
      df <- 2 * len
      p_val <- 1 - pchisq(cF, df = df)
      return(p_val)
    }
    
    percs <- floor(c(1:10) * 0.1 * no_of_modules) # logging output
    for (i in 1:no_of_modules) {
      module = names(M2G)[i]
      genes <- M2G[[i]]
      result[i, 1] <- module
      result[i, 2] <- MOD$Title[MOD$ID == module]
      present <- entrez %fin% genes
      N1 <- length(present[present == TRUE])
      cF <- -2 * sum(log(ranks[present]))
      result[i, 3] <- N1
      result[i, 4] <- cF
      if (N1 == 0) {
        result[i, 5] <- 0
        result[i, 8] <- 0
      } else {
        result[i, 5] <- cF / (2 * N1)
        N2 <- N - N1
        result[i, 8] <- ((N1 * N2) + ((N1 * (N1 + 1)) / 2) - sum(ranks[present])) / (N1 * N2)
      }
      
      if (i %fin% percs) {
        perc <- which(percs == i) * 10
        cat(paste(perc, "% ... \n", sep = ""))
      }
    }
    
    result$p <- pchisq(result[, 4], 2 * result[, 3], lower.tail = FALSE)
    result$p.adj <- p.adjust(p = result$p, method = adjust.method)
    cat("CERNO test completed. \n")
    return(result)
  }
