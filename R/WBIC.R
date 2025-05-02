WBIC <- function(x) {
  # x is an S x N matrix
  # S = number of posterior draws
  # N = number of observations
  n <- ncol(x)            # sample size N
  beta <- 1 / log(n)
  
  # Step 1: Compute total log-likelihood per posterior draw
  loglik_total <- rowSums(x)  # length S
  
  # Step 2: Compute importance weights
  log_weights <- (beta - 1) * loglik_total
  log_weights <- log_weights - max(log_weights)  # for numerical stability
  weights <- exp(log_weights)
  weights <- weights / sum(weights)             # normalize
  
  # Step 3: Estimate WBIC as twice the negative of the weighted average of total log-likelihood
  wbic <- -sum(2 * weights * loglik_total)
  return(wbic)
}
