postIC <- function(x, posterior = NULL) {
  ## Basic quantities
  S <- nrow(x)
  n <- ncol(x)
  beta <- 1 / log(n)
  deviance_draws <- -2 * rowSums(x)
  
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
  
  ## TIC
  dev_tic <- -2 * sum(colMeans(x))
  pD_tic  <-      sum(apply(x, 2, var))
  tic     <- dev_tic + 2 * pD_tic
  
  ## EIC
  ll_orig       <- colMeans(x)
  mean_ll_orig  <- mean(ll_orig)
  bias_resample <- numeric(nrow(x))
  for (i in 1:nrow(x)) {
    bias_resample[i] <- mean(x[i,]) - mean_ll_orig
  }
  bias_estimate <- mean(bias_resample)
  dev_eic       <- -2 * sum(ll_orig)
  eic           <- dev_eic - 2 * n * bias_estimate
  
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
  
  ## DIC family
  mean_dev        <- mean(deviance_draws)
  deviance_plugin <- which(min(deviance_draws))
  var_dev         <- var(deviance_draws)
  p_DIC           <- mean_dev - deviance_plugin
  DIC_classic     <- deviance_plugin + 2 * p_DIC
  p_V             <- 0.5 * var_dev               # variance‑based penalty
  DIC_p           <- deviance_plugin + 2 * p_V
  DIC_i           <- (DIC_classic + DIC_p)/2     # parameterization‑invariant DIC
  
  ## Output
  indices <- list(WBIC = wbic, WAIC = waic, p_WAIC = p_waic,
                  lppd = lppd, TIC = tic, EIC = eic, MDL = mdl,
                  DIC = DIC_classic, p_DIC = p_DIC,
                  DIC_p = DIC_p, p_V = p_V,
                  DIC_i = DIC_i)
  return(indices)
}