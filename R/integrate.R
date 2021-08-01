integrate <- function(data = NA,
                      method = "fisher",
                      adjust = TRUE,
                      adjust.method = "BH",
                      cor.LD = FALSE) {
  
  validate_input(data)
  validate_package("dplyr")
  
  ld <- function(p.value, len) {
    return (1 - (1 - p.value) ^ ((len + 1) * 0.5))
  }
  
  integrate.min <- function(x) {
    x <- x[!is.na(x)]
    len <- length(x)
    if (len == 0) {
      p.value <- NA
    } else {
      p.value <- min(x)
    }
    
    if (cor.LD == TRUE) {
      p.value <- ld(p.value, len)
    }
    
    return(p.value)
  }
  
  integrate.fisher <- function(x) {
    x <- x[!is.na(x)]
    len <- length(x)
    if (len == 0) {
      p.value <- NA
    } else {
      df <- 2 * len
      test <- -2 * sum(log(x))
      p.value <- 1 - pchisq(test, df = df)
      p.value <- if (is.nan(p.value))
        1.110223e-16
      else
        p.value
      p.value <- as.numeric(p.value)
    }
    
    if (cor.LD == TRUE) {
      p.value <- ld(p.value, len)
    }
    
    return(p.value)
  }
  
  integrate.stouffer <- function(x) {
    x <- x[!is.na(x)]
    len <- length(x)
    if (len == 0) {
      p.value <- NA
    } else {
      cdf <- qnorm(1 - x)
      cdf <- cdf[is.infinite(cdf) == F]
      cdf <- cdf[is.na(cdf) == F]
      len <- length(cdf)
      z <- sum(cdf) / sqrt(len)
      p.value <- 1 - pnorm(z)
      p.value <- if (is.nan(p.value))
        1.110223e-16
      else
        p.value
      p.value <- as.numeric(p.value)
    }
    
    if (cor.LD == TRUE) {
      p.value <- ld(p.value, len)
    }
    
    return(p.value)
  }
  
  if (grepl(method, "min", ignore.case = TRUE)) {
    fun = integrate.min
  } else if (grepl(method, "fisher", ignore.case = TRUE)) {
    fun = integrate.fisher
  } else if (grepl(method, "stouffer", ignore.case = TRUE)) {
    fun = integrate.stouffer
  } else {
    stop("Invalid method parameter.")
  }
  
  data = data[,-1]
  
  result <- data %>%
    group_by(entrez) %>%
    summarise(p.value = fun(p.value)) %>%
    as.data.frame()
  
  result$p.value[result$p.value == 0] <- 1.110223e-16
  
  if (adjust == TRUE) {
    result$p.value <- p.adjust(result$p.value, adjust.method)
  }
  
  return(result)
}
