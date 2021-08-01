gsea <- function(S, N) {
  is_hit <- S$present
  NH <- sum(is_hit)
  NR <- sum(abs(S$rank[is_hit]))
  P_hit <- rep(0, each = nrow(S))
  P_hit[is_hit] <- (abs(S[is_hit, "rank"]) / NR)
  miss_inc <- 1 / (N - NH)
  
  P_miss <- rep(miss_inc, each = nrow(S))
  P_miss[is_hit] <- 0
  S$p.miss <- cumsum(P_miss)
  S$p.hit <- cumsum(P_hit)
  S$ES <- cumsum(P_hit - P_miss)
  return(S)
}

permutation_matrix <-
  function(by,
           ranked,
           permutations,
           data.int = NA,
           int.method = "min",
           adjust.int = FALSE,
           adjust.method = "BH",
           int.cor.LD = FALSE) {
    if (by == "entrez") {
      permutation_by_entrez <- function(data, permutations) {
        set.seed(0)
        result <- data.frame(matrix(NA, nrow = nrow(data), ncol = permutations + 1))
        colnames(result) <- c("entrez", c(1:permutations))
        result[, 1] <- data$entrez
        for (i in 2:(permutations + 1)) {
          result[, i] <- sample(data$rank)
        }
        return(result)
      }
      
      permutation_by_entrez(ranked, permutations)
    } else if (by == "snp") {
      permutation_by_snp <-
        function(data,
                 permutations,
                 int.method,
                 adjust.int,
                 adjust.method,
                 int.cor.LD) {
          set.seed(0)
          result <- list()
          for (i in 1:permutations) {
            data$entrez <- sample(data$entrez)
            integrated <-
              integrate(
                data,
                method = int.method,
                adjust = adjust.int,
                adjust.method = adjust.method,
                cor.LD = int.cor.LD
              )
            l <-
              list(
                entrez = integrated$entrez,
                p.value = integrated$p.value,
                dataSPES = data
              )
            result <- append(result, list(l))
          }
          return(result)
        }
      
      permutation_by_snp(data.int,
                         permutations,
                         int.method,
                         adjust.int,
                         adjust.method,
                         int.cor.LD)
    } else if (by == "snp_raw") {
      permutation_by_snp_raw <- function(data, permutations) {
        set.seed(0)
        result <-
          data.frame(matrix(
            NA,
            nrow = nrow(data),
            ncol = (permutations + 1)
          ))
        colnames(result) <- c("snp", c(1:permutations))
        result[, 1] <- data[, 1]
        for (i in 2:(permutations + 1)) {
          result[, i] <- sample(data$entrez)
        }
        return(result)
      }
      
      permutation_by_snp_raw(ranked, permutations)
    } else {
      stop("Invalid 'by' parameter: correct values are 'snp', 'snp_raw', or 'entrez'.")
    }
  }

entrez_ranks <- function(data, samples, idx) {
  data$rank <- samples[, idx + 1]
  return(data)
}

snp_ranks <- function(data, samples, idx) {
  data$entrez <- samples[[idx]]$entrez
  data$p.value <- samples[[idx]]$p.value
  data$rank <- -log10(data$p.value)
  return(data)
}

calculate_nes <- function(no_of_modules, ES_NES_perm, result) {
  NES_perm_pos <- ES_NES_perm[ES_NES_perm > 0]
  NES_perm_neg <- ES_NES_perm[ES_NES_perm < 0]
  ind_pos = result[, 3] >= 0
  for (i in 1:no_of_modules) {
    if (result[i, 3] >= 0) {
      result[i, 7] = max(sum(NES_perm_pos >= result[i, 6]) / length(NES_perm_pos),
                         1 / length(NES_perm_pos))
      result[i, 8] = min(result[i, 7] / (sum(result[ind_pos, 6] >= result[i, 6]) /
                                           sum(ind_pos)), 1)
    } else {
      result[i, 7] = max(sum(NES_perm_neg <= result[i, 6]) / length(NES_perm_neg),
                         1 / length(NES_perm_neg))
      result[i, 8] = min(result[i, 7] / (sum(result[ind_pos == FALSE, 6] <= result[i, 6]) /
                                           sum(ind_pos == FALSE)), 1)
    }
  }
  return(result)
}
