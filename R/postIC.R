postIC <- function(x, posterior = NULL) {
  ## Basic quantities
  S <- nrow(x)
  n <- ncol(x)
  beta <- 1 / log(n)
  
  ## WBIC
  loglik_total <- rowSums(x)
  log_weights <- (beta - 1) * loglik_total
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  weights <- weights / sum(weights)
  wbic <- -sum(2 * weights * loglik_total)
  
  ## WAIC
  lppd <- sum(apply(x, 2, function(v) {
    m <- max(v)
    m + log(mean(exp(v - m)))
  }))
  p_waic <- sum(apply(x, 2, var))
  waic <- -2 * (lppd - p_waic)
  
  ## Bayesian MDL
  mdl <- NA
  if (!is.null(posterior)) {
    d <- ncol(posterior)
    E_loglik <- mean(loglik_total)
    sigma <- cov(posterior)
    logdet_sigma <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
    mdl <- E_loglik - (d / 2) * log(n) + (1/2) * logdet_sigma + (d / 2) * log(2 * pi * exp(1))
    mdl <- -2 * mdl
  }
  
  ## Output
  indices <- list(WBIC = wbic, WAIC = waic, p_WAIC = p_waic,
                  lppd = lppd, MDL = mdl)
  return(indices)
}