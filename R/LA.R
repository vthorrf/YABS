LA <- function(Model, Data, Initial.Values=NULL, lower = -Inf, upper = Inf,
               control = list(), hessian = FALSE, par.cov=NULL, SIR=TRUE, 
               method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
               iterations=NULL) {
  ################=============== Initial settings
  ## Default values for the arguments and some error handling
  if(length(method) > 1)       method          <- "BFGS"
  if(is.null(Initial.Values))  Initial.Values  <- Data$PGF(Data)
  if(is.null(iterations))      iterations      <- 1000
  if(iterations <= 0)          iterations      <- 1000
  if(!is.matrix(par.cov) & !is.null(par.cov)) stop("'par.cov' should be a matrix. Please check the documentation.")
  
  ################=============== Performing Laplace Approximation
  startTime = proc.time()
  ## Objetive function
  post <- function(par) return( suppressWarnings(-Model(par, Data)$LP) )
  ## Optimization routine
  fit  <- optim( par=Initial.Values, fn=post, hessian=hessian, lower=lower,
                 upper=upper, control=control, method=method )
  MAP <- fit$par
  names(MAP) <- Data[["parm.names"]]
  ## Parameters' covariance matrix
  if(is.null(par.cov) & is.null(fit$hessian) & hessian) {
    fit$hessian <- optimHess(MAP, fn=post, control=control)
  }
  if(is.null(par.cov) & hessian) {
    VarCov <- -solve(fit$hessian)
    diag(VarCov) <- abs(diag(VarCov))
  } else if(is.null(par.cov) & !hessian) {
    VarCov <-tryCatch( solve(-diag(gradN(Model, Data, MAP, order=2)*1e6)),
                       error=function(e) ginv(-diag(gradN(Model, Data, MAP, order=2)*1e6)) )
  }else if(!is.null(par.cov)) {
    VarCov <- par.cov
  }
  
  ################=============== Sampling Importance Resampling
  if(SIR) {
    # Sample values from the posterior
    theta <- mvrnorm(iterations, MAP, Sigma=VarCov, empirical=TRUE)
    # New values
    yhat <- matrix(NA, nrow=iterations, ncol=length(Model(theta[1, ], Data)[["yhat"]]))
    LP <- vector("numeric", length=iterations)
    Dev <- vector("numeric", length=iterations)
    # Calculate log-posterior, deviance, and predictions
    for (i in 1:iterations) {
      temp <- Model(theta[i, ], Data)
      LP[i] <- temp[["LP"]]
      Dev[i] <- temp[["Dev"]]
      yhat[i,] <- temp[["yhat"]]
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
  } else {
    temp <- Model(MAP, Data)
    LP <- temp[["LP"]]
    Dev <- temp[["Dev"]]
    AIC <- Dev + 2 * length(MAP)
    yhat <- temp[["yhat"]]
    
  }
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  
  ################=============== Final list of results
  ## Calculate DIC
  if(SIR) {
    Dbar      <- mean(Dev)
    pD        <- var(Dev)/2
    DIC       <- list(DIC=Dbar + pD, Dbar=Dbar, pD=pD)
  }
  ## Effective Sample Size
  ESS       <- apply(posterior,2,effectiveSize)
  ## List of results
  if(SIR) {
    Result    <- list(posterior=posterior,
                      MAP=MAP,
                      yhat=yhat[indices,],
                      LP=LP[indices],
                      Dev=Dev,
                      AIC=AIC,
                      DIC=DIC,
                      ESS=ESS,
                      convergence={fit$convergence == 0})
  } else {
    Result    <- list(MAP=MAP,
                      yhat=yhat,
                      LP=LP,
                      Dev=Dev,
                      AIC=AIC,
                      convergence={fit$convergence == 0})
  }
  class(Result) <- "YABS_LA"
  return(Result)
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