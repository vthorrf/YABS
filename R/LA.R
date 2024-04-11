LA <- function(Model, Data, Initial.Values=NULL, par.cov=NULL, SIR=TRUE,
               iterations=NULL, nearPD=TRUE, check.convergence=FALSE, 
               method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS", "L-BFGS-B", "SANN", "Brent", "nlm", "nlminb"),
               lower = -Inf, upper = Inf, control = list(), hessian = FALSE) {
  ################=============== Initial settings
  ## Default values for the arguments and some error handling
  if(length(method) > 1)       method          <- "nlminb"
  if(is.null(Initial.Values))  Initial.Values  <- Data$PGF(Data)
  if(is.null(iterations))      iterations      <- 1000
  if(iterations <= 0)          iterations      <- 1000
  if(!is.matrix(par.cov) & !is.null(par.cov)) stop("'par.cov' should be a matrix. Please check the documentation.")
  
  ################=============== Performing Laplace Approximation
  startTime = proc.time()
  cat("Initial optimization using the", method, "method\n")
  ## Objetive function
  post <- function(par) return( suppressWarnings(-Model(par, Data)$LP) )
  ## Optimization routine
  if(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS", "L-BFGS-B", "SANN", "Brent")) {
    fit  <- optim( par=Initial.Values, fn=post, hessian=hessian, lower=lower,
                   upper=upper, control=control, method=method )
  } else if(method == "nlm") {
    fit <- suppressWarnings( nlm( p=Initial.Values, f=post, hessian=hessian ) )
    if(fit$code %in% c(1,2)) {
      fit$convergence <- 0
    } else {
      fit$convergence <- 1
    }
  } else if(method == "nlminb") {
    fit <- suppressWarnings( nlminb( start=Initial.Values, objective=post, control=control,
                             lower=lower, upper=upper ) )
  } else {
    stop("Unknown optimization method! Please check the documentation")
  }
  ## Check if convergence was achieved
  if({fit$convergence != 0} & check.convergence) {
    continue <- readline("Convergence was not achieved. Do you want to continue anyway? (Y/N) ")
    resp <- abs(regexpr(continue, 'n', ignore.case = TRUE) - regexpr(continue, 'y', ignore.case = TRUE))
    while(resp == 0) {
      continue <- readline("Please respond Y if you want to continue or N if you want to stop: ")
      resp <- abs(regexpr(continue, 'n', ignore.case = TRUE) - regexpr(continue, 'y', ignore.case = TRUE))
    }
    if(regexpr(continue, 'y', ignore.case = TRUE) == -1) {
      stop("Convergence was not achieved. Try another optimization method or running the algorithm for longer")
    }
  }
  ## MAP estimate
  if(method == "nlm") {
    MAP <- fit$estimate
  } else {
    MAP <- fit$par
  }
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
  if({min(eigen(VarCov)$values) <= 0} & nearPD) {
    VarCov <- nearPD(VarCov, base.matrix=TRUE, ensureSymmetry=TRUE)$mat
    colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
  }
  
  ################=============== Sampling Importance Resampling
  if(SIR) {
    cat("Start Sampling-Importance Resampling with", iterations,"samples\n")
    # New values
    samples <- samplingImportanceResampling(MAP, VarCov, Model, Data, iterations)
    indices <- samples[["indices"]]
    yhat <- samples[["yhat"]]
    Monitor <- samples[["Monitor"]]
    LP <- samples[["LP"]]
    Dev <- samples[["Dev"]]
    posterior <- samples[["posterior"]]
    colnames(posterior) <- Data[["parm.names"]]
  } else {
    cat("Start sampling a posterior with", iterations,"samples\n")
    posterior <- mvrnorm(iterations, MAP, Sigma=VarCov, empirical=TRUE)
    indices <- 1:nrow(posterior)
    yhat <- matrix(NA, nrow=iterations, ncol=length(Model(posterior[1, ], Data)[["yhat"]]))
    Monitor <- matrix(NA, nrow=iterations, ncol=length(Model(posterior[1, ], Data)[["Monitor"]]))
    LP <- vector("numeric", length=iterations)
    Dev <- vector("numeric", length=iterations)
    for (i in 1:iterations) {
      if(i == 1) {
        temp <- Model(MAP, Data)
      } else {
        temp <- Model(posterior[i, ], Data)
      }
      LP[i] <- temp[["LP"]]
      Dev[i] <- temp[["Dev"]]
      yhat[i,] <- temp[["yhat"]]
      Monitor[i,] <- temp[["Monitor"]]
    }
    #temp <- Model(MAP, Data)
    #LP <- temp[["LP"]]
    #Dev <- temp[["Dev"]]
    AIC <- Dev + 2 * length(MAP)
    #yhat <- temp[["yhat"]]
    #Monitor <- temp[["Monitor"]]
  }
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  
  ################=============== Final list of results
  ## Calculate DIC
  #if(SIR) {
    Dbar      <- mean(Dev)
    pD        <- var(Dev)/2
    DIC       <- list(DIC=Dbar + pD, Dbar=Dbar, pD=pD)
    ## Effective Sample Size
    ESS       <- apply(posterior,2,effectiveSize)
  #}
  cat("Done!\n")
  ## List of results
  #if(SIR) {
    Result    <- list(posterior=posterior,
                      MAP=MAP,
                      VarCov=VarCov,
                      yhat=yhat[indices,],
                      LP=LP[indices],
                      Monitor=Monitor[indices,],
                      Dev=Dev,
                      AIC=AIC,
                      DIC=DIC,
                      ESS=ESS,
                      convergence={fit$convergence == 0})
  #} else {
  #  Result    <- list(MAP=MAP,
  #                    yhat=yhat,
  #                    LP=LP,
  #                    Monitor=Monitor,
  #                    Dev=Dev,
  #                    AIC=AIC,
  #                    convergence={fit$convergence == 0},
  #                    VarCov=VarCov)
  #}
  class(Result) <- "YABS_LA"
  return(Result)
}