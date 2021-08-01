#' GSEASNP test implementation
#'
#' Performs an GSEASNP test on SNP p-values
#'
#' @param data input dataframe with column names: "snp", "entrez", "p.value"
#' @param pathways object containing a list of GENES and a list of MODULES
#' @param adjust.method multiple testing correction method
#' @param permutations number of gene set permutations
#' @param cores number of cores; default is 2
#'
#' @return A data frame with module names, calculated p-value, and additional statistics. 
#'
#' @examples
#' GSEASNPtest(data = example_dataset, pathways = pathway_library, adjust.method = "BH", 
#' permutations = 10000)
#'
#' @export
GSEASNPtest <- function(data,
                        pathways = NA,
                        adjust.method = "BH",
                        permutations = 10000,
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
  
  cat("Initialization... \n")
  N <- nrow(data)
  ranked <- data[order(data$p.value),]
  ranked$rank <- -log10(ranked$p.value)
  
  cat("Generating permutation matrix... \n")
  samples <-
    permutation_matrix(by = "snp_raw",
                       ranked = ranked,
                       permutations = permutations)
  
  cat("GSEA-SNP test... \n")
  registerDoParallel(cores = cores)
  pb <- txtProgressBar(max = 341, style = 3)
  progress <- function(n)
    setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  result_para <-
    foreach (
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
      S$present = S$entrez %fin% genes
      S <- gsea(S, N)
      ES = S$ES[which.max(abs(S$ES))]
      result[1, 3] <- ES
      
      ESs <- c()
      for (k in 1:permutations) {
        S[c] <- NA
        S$entrez <- samples[, k + 1]
        S$present = S$entrez %fin% genes
        S <- gsea(S, N)
        ESs <- c(ESs, S$ES[which.max(abs(S$ES))])
      }
      
      if (ES >= 0) {
        poz <- which(ESs > 0)
        ES_NES_perm[1, poz] <- ESs[poz] / abs(mean(ESs[poz]))
        ESs <- ESs[poz]
        count <- sum(ESs >= ES)
      } else {
        poz <- which(ESs < 0)
        ES_NES_perm[1, poz] <- ESs[poz] / abs(mean(ESs[poz]))
        ESs <- ESs[poz]
        count <- sum(ESs <= ES)
      }
      
      result[1, 4] <- count / length(ESs)
      result[1, 6] <- ES / abs(mean(ESs))
      
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
  stopCluster(cl)
  
  cat("\n")
  cat("Calculate NES... \n")
  result <- calculate_nes(no_of_modules, ES_NES_perm, result)
  result$p.value.adj <- p.adjust(p = result$p.value, method = adjust.method)
  cat("GSEASNP test completed. \n")
  return(result)
}
