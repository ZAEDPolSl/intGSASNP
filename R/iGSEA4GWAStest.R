#' iGSEA4GWAS test implementation
#'
#' Performs an iGSEA4GWAS test on SNP p-values
#'
#' @param data input dataframe with column names: "snp", "entrez", "p.value"
#' @param pathways object containing a list of GENES and a list of MODULES
#' @param int.method integration method (either "min", "fisher', or "stouffer")
#' @param adjust.int flag for multiple testing correction in integration
#' @param int.cor.LD flag for LD correction in integration
#' @param adjust.method multiple testing correction method
#' @param permutations number of gene set permutations
#' @param cores number of cores; if cores > 2 use parallel computing
#'
#' @return A data frame with module names, calculated p-value, and additional statistics. 
#'
#' @examples
#' iGSEA4GWAStest(data = example_dataset, pathways = pathway_library, int.method = "min", 
#' adjust.int = FALSE, int.cor.LD = FALSE, adjust.method = "BH", permutations = 10000)
#'
#' @export
iGSEA4GWAStest <- function(data,
                           pathways = NA,
                           int.method = "min",
                           adjust.int = FALSE,
                           int.cor.LD = FALSE,
                           adjust.method = "BH",
                           permutations = 1000,
                           cores = 2) {
  
  
  iGSEA4GWAStest_parallel <-
    function(data,
             pathways = NA,
             int.method = "min",
             adjust.int = FALSE,
             int.cor.LD = FALSE,
             adjust.method = "BH",
             permutations = 1000,
             cores = 2) {
      
      validate_input(data)
      validate_pathways(pathways)
      validate_package("fastmatch")
      validate_package("foreach")
      validate_package("doSNOW")
      validate_package("doParallel")
      
      M2G = pathways$MODULES2GENES
      no_of_modules <- length(M2G)
      MOD = pathways$MODULES
      
      data_cpy <- data
      data_cpy <- data_cpy[order(data_cpy$p.value), ]
      
      cat("Performing p-value integration... \n")
      data <-
        integrate(data,
                  method = int.method,
                  adjust = adjust.int,
                  cor.LD = int.cor.LD)
      N <- nrow(data)
      ranked <- data[order(data$p.value),]
      ranked$rank <- -log10(ranked$p.value)
      
      cat("Generating permutation matrix... \n")
      samples <-
        permutation_matrix(
          by = "snp",
          ranked,
          permutations,
          data_cpy,
          int.method,
          adjust.int,
          adjust.method,
          int.cor.LD
        )
      get_new_ranks <- snp_ranks
      
      cat("iGSEA4GWAS test... \n")
      registerDoParallel(cores = cores)
      pb <- txtProgressBar(max = 341, style = 3)
      progress <- function(n)
        setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      result_para <- foreach (
        i = 1:no_of_modules,
        .combine = 'rbind',
        .options.snow = opts,
        .packages = c("fastmatch"),
        .export = c(
          "gsea",
          "integrate",
          "permutation_matrix",
          "snp_ranks",
          "phenotype_ranks"
        )
      ) %dopar% {
        result <- data.frame(matrix(NA, nrow = 1, ncol = 8))
        
        module = names(M2G)[i]
        genes = M2G[[i]]
        result[1, 1] <- module
        result[1, 2] <- MOD$Title[MOD$ID == module]
        S <- ranked
        c = c("p.miss", "p.hit", "ES")
        S[c] <- NA
        S$present = ranked$entrez %fin% genes
        S <- gsea(S, N)
        ES = S$ES[which.max(abs(S$ES))]
        result[1, 3] <- ES
        
        ESs <- matrix(NA, ncol = permutations, nrow = 2)
        for (k in 1:permutations) {
          S.art = S
          S.art[c] <- NA
          S.art <- get_new_ranks(S.art, samples, k)
          S.art <- S.art[order(S.art$rank, decreasing = TRUE),]
          S.art$present = S.art$entrez %fin% genes
          S.art <- gsea(S.art, N)
          ESs[1, k] = S.art$ES[which.max(abs(S.art$ES))]
          
          data_c <- samples[[k]]$dataSPES
          data_c$present <- data_c$entrez %fin% genes
          cutoff <- round(0.05 * nrow(data_c))
          K <- length(unique(data_c$entrez[1:cutoff]))
          data_c <- data_c[data_c$present,]
          cutoff <- round(0.05 * nrow(data_c))
          kk <- length(unique(data_c$entrez[1:cutoff]))
          ESs[2, k] <- ESs[1, k] * kk / K #SPES
        }
        
        SPESs <- ESs[2, ]
        ESs <- ESs[1, ]
        
        if (ES >= 0) {
          ESs <- ESs[which(ESs > 0)]
          count <- sum(ESs >= ES)
        } else {
          ESs <- ESs[which(ESs < 0)]
          count <- sum(ESs <= ES)
        }
        
        result[1, 4] <- count / length(ESs)  # p-value ES
        
        # SPES
        cutoff <- round(0.05 * nrow(data_cpy))
        K <- length(unique(data_cpy$entrez[1:cutoff]))
        S.c <- data_cpy[data_cpy$entrez %fin% S$entrez[S$present == TRUE],]
        cutoff <- round(0.05 * nrow(S.c))
        k <- length(unique(S.c$entrez[1:cutoff]))
        SPES <- ES * k / K
        result[1, 6] <- SPES
        
        if (ES >= 0) {
          SPESs <- SPESs[which(SPESs > 0)]
          count <- sum(SPESs >= SPES)
        } else {
          SPESs <- SPESs[which(SPESs < 0)]
          count <- sum(SPESs <= SPES)
        }
        
        result[1, 7] <- count / length(SPESs)
        
        result
      }
      
      result <- result_para
      colnames(result) <-
        c(
          "ID",
          "Title",
          "ES",
          "p.value",
          "p.value.adj",
          "SPES",
          "SPES.p.val",
          "SPES.p.val.adj"
        )
      
      
      result$p.value.adj <- p.adjust(p = result$p.value, method = adjust.method)
      result$SPES.p.val.adj <- p.adjust(p = result$SPES.p.val, method = adjust.method)
      cat("\n")
      cat("iGSEA4GWAS test completed. \n")
      return(result)
    }
  
  iGSEA4GWAStest_serial <-
    function(data,
             pathways = NA,
             int.method = "min",
             adjust.int = FALSE,
             int.cor.LD = FALSE,
             adjust.method = "BH",
             permutations = 1000) {
      
      validate_input(data)
      validate_pathways(pathways)
      validate_package("fastmatch")
      
      M2G = pathways$MODULES2GENES
      no_of_modules <- length(M2G)
      MOD = pathways$MODULES
      
      data_cpy <- data
      data_cpy <- data_cpy[order(data_cpy$p.value), ]
      
      cat("Performing p-value integration... \n")
      data <-
        integrate(data,
                  method = int.method,
                  adjust = adjust.int,
                  cor.LD = int.cor.LD)
      N <- nrow(data)
      ranked <- data[order(data$p.value),]
      ranked$rank <- -log10(ranked$p.value)
      
      cat("Generating permutation matrix... \n")
      samples <-
        permutation_matrix(
          by = "snp",
          ranked,
          permutations,
          data_cpy,
          int.method,
          adjust.int,
          adjust.method,
          int.cor.LD
        )
      get_new_ranks <- snp_ranks
      
      cat("iGSEA4GWAS test... \n")
      
      result <- data.frame(matrix(NA, nrow = 0, ncol = 8))
      colnames(result) <- c("ID",
                            "Title",
                            "ES",
                            "p.value",
                            "p.value.adj",
                            "SPES",
                            "SPES.p.val",
                            "SPES.p.val.adj")
      percs <- floor(c(1:10) * 0.1 * no_of_modules) # logging output
      for (i in 1:no_of_modules) {
        module = names(M2G)[i]
        genes = M2G[[i]]
        result[i, 1] <- module
        result[i, 2] <- MOD$Title[MOD$ID == module]
        S <- ranked
        c = c("p.miss", "p.hit", "ES")
        S[c] <- NA
        S$present = ranked$entrez %fin% genes
        S <- gsea(S, N)
        ES = S$ES[which.max(abs(S$ES))]
        result[i, 3] <- ES
        
        ESs <- matrix(NA, ncol = permutations, nrow = 2)
        for (k in 1:permutations) {
          S.art = S
          S.art[c] <- NA
          S.art <- get_new_ranks(S.art, samples, k)
          S.art <- S.art[order(S.art$rank, decreasing = TRUE),]
          S.art$present = S.art$entrez %fin% genes
          S.art <- gsea(S.art, N)
          ESs[1, k] = S.art$ES[which.max(abs(S.art$ES))]
          
          data_c <- samples[[k]]$dataSPES
          data_c$present <- data_c$entrez %fin% genes
          cutoff <- round(0.05 * nrow(data_c))
          K <- length(unique(data_c$entrez[1:cutoff]))
          data_c <- data_c[data_c$present,]
          cutoff <- round(0.05 * nrow(data_c))
          kk <- length(unique(data_c$entrez[1:cutoff]))
          ESs[2, k] <- ESs[1, k] * kk / K #SPES
        }
        
        SPESs <- ESs[2, ]
        ESs <- ESs[1, ]
        
        if (ES >= 0) {
          ESs <- ESs[which(ESs > 0)]
          count <- sum(ESs >= ES)
        } else {
          ESs <- ESs[which(ESs < 0)]
          count <- sum(ESs <= ES)
        }
        
        result[i, 4] <- count / length(ESs)  # p-value ES
        
        # SPES
        cutoff <- round(0.05 * nrow(data_cpy))
        K <- length(unique(data_cpy$entrez[1:cutoff]))
        S.c <- data_cpy[data_cpy$entrez %fin% S$entrez[S$present == TRUE], ]
        cutoff <- round(0.05 * nrow(S.c))
        k <- length(unique(S.c$entrez[1:cutoff]))
        SPES <- ES * k / K
        result[i, 6] <- SPES
        
        if (ES >= 0) {
          SPESs <- SPESs[which(SPESs > 0)]
          count <- sum(SPESs >= SPES)
        } else {
          SPESs <- SPESs[which(SPESs < 0)]
          count <- sum(SPESs <= SPES)
        }
        
        result[i, 7] <- count / length(SPESs)
        
        if (i %fin% percs) {
          perc <- which(percs == i) * 10
          cat(paste(perc, "% ... \n", sep = ""))
        }
      }
      
      result$p.value.adj <- p.adjust(p = result$p.value, method = adjust.method)
      result$SPES.p.val.adj <- p.adjust(p = result$SPES.p.val, method = adjust.method)
      cat("iGSEA4GWAS test completed. \n")
      return(result)
    }
  
  if (cores > 1) {
    result <- iGSEA4GWAStest_parallel(
      data,
      pathways,
      int.method,
      adjust.int,
      int.cor.LD,
      adjust.method,
      permutations,
      cores
    )
  } else if (cores == 1) {
    result <- iGSEA4GWAStest_serial(data,
                                    pathways,
                                    int.method,
                                    adjust.int,
                                    int.cor.LD,
                                    adjust.method,
                                    permutations)
  } else {
    stop("Incorrect number of cores.")
  }
  
  return(result)
  
}