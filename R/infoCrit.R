infoCrit <- function(posterior, deviance, sampleSize, infoC=T) {
  ## Initial settings
  posterior <- as.matrix(posterior)
  n         <- sampleSize
  d         <- ncol(posterior)
  C         <- {n^{2/3}}
  if(infoC) {
    sigma <- cov(posterior)
    FC <- function(sigma) {
      {{1/d} * tr(sigma %*% t(sigma))} - {{tr(sigma)/d}^2}
    }
    tr <- function(sigma) sum(diag(sigma))
    CC <- function(sigma) {
      {.5 * sum(log(diag(sigma)))} - {.5 * log(det(sigma))}
    }
    C1 <- function(sigma) {
      {d*log(tr(sigma)/d)/2} - {.5 * log(det(sigma))}
    }
  }
 
  ## Fit measures
  # Akaike
  aic  <- deviance + 2 * d
  # Finite-sample corrected Akaike 
  aicc <- aic + {2*{d + 1}*{d + 2}}/{n - d - 2}
  # Bridge
  bc   <- deviance + C * sum(1/{1:d})
  # Bayesian (Schwarz)
  bic  <- deviance + d * log(n)
  # Hannan and Quinn
  hqc  <- deviance + 2 * d * log(log(n))
  if(infoC) {
    # Frobenius complexity
    fcomp <- deviance + 2 * FC(sigma)
    # Informational complexity
    icomp <- deviance + 2 * CC(sigma)
    # Cramer-Rao lower bound
    ifim  <- deviance + 2 * C1(sigma)
    # Posterior expected utility
    peu   <- deviance + d + 2 * C1(sigma)
    # Posterior expected utility Log-sampleSize
    peuln <- deviance + d + 2 * log(n) * C1(sigma)
    # Finite-sample corrected Posterior expected utility
    peums <- deviance + d + {2 * {{{n * d}/{n - d - 2}} + C1(sigma)}}
    # Minimum Description Length
    logdet_sigma <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
    mdl <- deviance + (d * log(n)) - logdet_sigma - {d * log(2 * pi * exp(1))}
    # Stochastic complexity 2
    sc2   <- deviance + {d*log(n/{2*pi})} - logdet_sigma
    # Stochastic complexity 3
    sc3   <- deviance + {d*log(n/{2*pi})}
  }
  
  ## Final results
  if(infoC) {
    results <- data.frame("aic"=aic, "aicc"=aicc, "bc"=bc, "bic"=bic, "hqc"=hqc,
                          "fcomp"=fcomp, "icomp"=icomp, "ifim"=ifim, "peu"=peu,
                          "peuln"=peuln, "peums"=peums, "mdl"=mdl, "sc2"=sc2, "sc3"=sc3)
  } else {
    results <- data.frame("aic"=aic, "aicc"=aicc, "bc"=bc, "bic"=bic, "hqc"=hqc)
  }
  return(results)
}
