#' GSEAtest test implementation
#'
#' Performs an GSEAtest test on SNP p-values
#'
#' @param data input dataframe with column names: "snp", "entrez", "p.value"
#' @param pathways object containing a list of GENES and a list of MODULES
#' @param int.method integration method (either "min", "fisher', or "stouffer")
#' @param adjust.int flag for multiple testing correction in integration
#' @param int.cor.LD flag for LD correction in integration
#' @param adjust.method multiple testing correction method
#' @param permutations number of gene set permutations
#' @param permutation.method permutation method (either "entrez", "snp", or "snp_raw")
#' @param cores number of cores; if cores > 2 use parallel computing
#'
#' @return A data frame with module names, calculated p-value, and additional statistics. 
#'
#' @examples
#' GSEAtest(data = example_dataset, pathways = pathway_library, int.method = "min", 
#' adjust.int = FALSE, int.cor.LD = FALSE, adjust.method = "BH", permutations = 10000, 
#' permutation.method = "entrez", cores = 2)
#'
#' @export
GSEAtest <- function(data,
                     pathways = NA,
                     int.method = "min",
                     adjust.int = FALSE,
                     int.cor.LD = FALSE,
                     adjust.method = "BH",
                     permutations = 1000,
                     permutation.method = "entrez",
                     cores = 2) {
    
    # data input: snp | entrez | p.value
    GSEAtest_parallel <-
      function(data,
               pathways = NA,
               int.method = "min",
               adjust.int = FALSE,
               int.cor.LD = FALSE,
               adjust.method = "BH",
               permutations = 1000,
               permutation.method = "entrez",
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
        
        data.int <- data
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
            by = permutation.method,
            ranked,
            permutations,
            data.int,
            int.method,
            adjust.int,
            adjust.method,
            int.cor.LD
          )
        get_new_ranks <-
          if (permutation.method == "entrez")
            entrez_ranks
        else
          snp_ranks
        
        cat("GSEA test... \n")
        registerDoParallel(cores = cores)
        pb <- txtProgressBar(max = 341, style = 3)
        progress <- function(n)
          setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        result_para <- foreach(
          i = 1:no_of_modules,
          .combine = 'rbind',
          .options.snow = opts,
          .packages = c("fastmatch"),
          .export = "gsea"
        ) %dopar% {
          result <- data.frame(matrix(NA, nrow = 1, ncol = 8))
          ES_NES_perm <- matrix(0, nrow = 1, ncol = permutations)
          
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
          
          ESs <- c()
          for (k in 1:permutations) {
            S.art = S
            S.art[c] <- NA
            S.art <- get_new_ranks(S.art, samples, k)
            S.art <- S.art[order(S.art$rank, decreasing = TRUE),]
            S.art$present = S.art$entrez %fin% genes
            S.art <- gsea(S.art, N)
            ES.art = S.art$ES[which.max(abs(S.art$ES))]
            ESs <- c(ES.art, ESs)
          }
          
          if (ES >= 0) {
            pos <- which(ESs > 0)
            ES_NES_perm[1, pos] <- ESs[pos] / abs(mean(ESs[pos]))
            ESs <- ESs[pos]
            count <- sum(ESs >= ES)
          } else {
            pos <- which(ESs < 0)
            ES_NES_perm[1, pos] <- ESs[pos] / abs(mean(ESs[pos]))
            ESs <- ESs[pos]
            count <- sum(ESs <= ES)
          }
          
          result[1, 4] <- count / length(ESs)  # p-value ES
          result[1, 6] <- ES / abs(mean(ESs)) # NES
          
          out <- cbind(result, ES_NES_perm)
          out
        }
        
        result <- result_para[, 1:8]
        colnames(result) <-
          c("ID",
            "Title",
            "ES",
            "p.value",
            "p.value.adj",
            "NES",
            "p.NES",
            "q.NES")
        
        ES_NES_perm <- result_para[, 9:ncol(result_para)]
        
        cat("\n")
        cat("Calculate NES... \n")
        result <- calculate_nes(no_of_modules, ES_NES_perm, result)
        result$p.value.adj <- p.adjust(p = result$p.value, method = adjust.method)
        cat("GSEA test completed. \n")
        return(result)
      }
    
    GSEAtest_serial <-
      function(data,
               pathways = NA,
               int.method = "min",
               adjust.int = FALSE,
               int.cor.LD = FALSE,
               adjust.method = "BH",
               permutations = 1000,
               permutation.method = "entrez") {
        
        validate_input(data)
        validate_pathways(pathways)
        validate_package("fastmatch")
        
        M2G = pathways$MODULES2GENES
        no_of_modules <- length(M2G)
        MOD = pathways$MODULES
        
        data.int <- data
        data <-
          integrate(data,
                    method = int.method,
                    adjust = adjust.int,
                    cor.LD = int.cor.LD)
        N <- nrow(data)
        ranked <- data[order(data$p.value),]
        ranked$rank <- -log10(ranked$p.value)
        
        result <- data.frame(matrix(NA, nrow = 0, ncol = 8))
        colnames(result) <-
          c("ID",
            "Title",
            "ES",
            "p.value",
            "p.value.adj",
            "NES",
            "p.NES",
            "q.NES")
        
        cat("Generating permutation matrix... \n")
        samples <-
          permutation_matrix(
            by = permutation.method,
            ranked,
            permutations,
            data.int,
            int.method,
            adjust.int,
            adjust.method,
            int.cor.LD
          )
        get_new_ranks <-
          if (permutation.method == "entrez")
            entrez_ranks
        else
          snp_ranks
        
        cat("GSEA test... \n")
        ES_NES_perm <- matrix(0, nrow = no_of_modules, ncol = permutations)
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
          
          ESs <- c()
          for (k in 1:permutations) {
            S.art = S
            S.art[c] <- NA
            S.art <- get_new_ranks(S.art, samples, k)
            S.art <- S.art[order(S.art$rank, decreasing = TRUE),]
            S.art$present = S.art$entrez %fin% genes
            S.art <- gsea(S.art, N)
            ES.art = S.art$ES[which.max(abs(S.art$ES))]
            ESs <- c(ES.art, ESs)
          }
          
          if (ES >= 0) {
            pos <- which(ESs > 0)
            ES_NES_perm[i, pos] <- ESs[pos] / abs(mean(ESs[pos]))
            ESs <- ESs[pos]
            count <- sum(ESs >= ES)
          } else {
            pos <- which(ESs < 0)
            ES_NES_perm[i, pos] <- ESs[pos] / abs(mean(ESs[pos]))
            ESs <- ESs[pos]
            count <- sum(ESs <= ES)
          }
          
          result[i, 4] <- count / length(ESs)  # p-value ES
          result[i, 6] <- ES / abs(mean(ESs)) # NES
          
          if (i %fin% percs) {
            perc <- which(percs == i) * 10
            cat(paste(perc, "% ... \n", sep = ""))
          }
        }
        
        cat("Calculate NES... \n")
        result <- calculate_nes(no_of_modules, ES_NES_perm, result)
        result$p.value.adj <- p.adjust(p = result$p.value, method = adjust.method)
        cat("GSEA test completed. \n")
        return(result)
      }
    
    if (cores > 1) {
      result <- GSEAtest_parallel(data,
                                  pathways,
                                  int.method,
                                  adjust.int,
                                  int.cor.LD,
                                  adjust.method,
                                  permutations,
                                  permutation.method,
                                  cores)
    } else if (cores == 1) {
      result <- GSEAtest_serial(data,
                                pathways,
                                int.method,
                                adjust.int,
                                int.cor.LD,
                                adjust.method,
                                permutations,
                                permutation.method)
    } else {
      stop("Incorrect number of cores.")
    }
    
    return(result)
  }