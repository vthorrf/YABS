samplingImportanceResampling <- function(MAP, VarCov, Model, Data, iterations) {
  # Sample values from the posterior
  theta <- mvrnorm(iterations, MAP, Sigma=VarCov, empirical=TRUE)
  # New values
  yhat <- matrix(NA, nrow=iterations, ncol=length(Model(theta[1, ], Data)[["yhat"]]))
  Monitor <- matrix(NA, nrow=iterations, ncol=length(Model(theta[1, ], Data)[["Monitor"]]))
  LP <- vector("numeric", length=iterations)
  Dev <- vector("numeric", length=iterations)
  # Calculate log-posterior, deviance, and predictions
  for (i in 1:iterations) {
    if(i == 1) {
      temp <- Model(MAP, Data)
    } else {
      temp <- Model(theta[i, ], Data)
    }
    LP[i] <- temp[["LP"]]
    Dev[i] <- temp[["Dev"]]
    yhat[i,] <- temp[["yhat"]]
    Monitor[i,] <- temp[["Monitor"]]
  }
  # Calculate AIC
  AIC        <- Dev + 2 * length(MAP)
  # Calculate the densities for the samples
  LMVN       <- dmnorm(x=theta, mu=MAP, sigma=VarCov, log=TRUE)
  # Calculate the log weight of the samples
  MaxDens    <- max(LP - LMVN)
  LogWeights <- LP - LMVN - MaxDens
  if (any(!is.finite(LogWeights))) {
    LogWeights[!is.finite(LogWeights)] <- min(LogWeights[is.finite(LogWeights)])
  }
  # Resample the posterior samples
  probs <- exp(LogWeights - logadd(LogWeights))
  indices <- sample(1:iterations, size=iterations, replace=TRUE, prob=probs)
  posterior <- theta[indices, ]
  colnames(posterior) <- Data[["parm.names"]]
  # Final results
  return(list("indices"=indices, "yhat"=yhat, "Monitor"=Monitor,
              "LP"=LP, "Dev"=Dev, "posterior"=posterior))
}

logadd <- function (x, add = TRUE) {
  x <- as.vector(x)
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  x <- c(x[1], x[which(x != x[1])])
  if (length(x) == 1)  {
    return(x)
  }
  n <- length(x)
  if (add == TRUE) {
    z <- x[1] + log(1 + sum(exp(x[-1] - x[1])))
  } else {
    z <- x[1] + sum(log(1 - exp(x[-1] - x[1])))
  }
  return(z)
}