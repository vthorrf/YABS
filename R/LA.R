LA <- function(Model, Data, Initial.Values = NULL, par.cov = NULL, SIR = FALSE,
               iterations = NULL, nearPD = TRUE, check.convergence = FALSE,
               method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS", "L-BFGS-B",
                          "SANN", "Brent", "nlm", "nlminb"),
               lower = -Inf, upper = Inf, control = list(), hessian = FALSE) {
  ################=============== Initial settings
  if (length(method) > 1) method <- "nlminb"
  if (is.null(Initial.Values)) Initial.Values <- Data$PGF(Data)
  if (is.null(iterations) || iterations <= 0) iterations <- 1000
  if (!is.null(par.cov) && !is.matrix(par.cov)) {
    stop("'par.cov' must be a covariance matrix.")
  }
  startTime <- proc.time()
  cat("Initial optimization using the", method, "method\n")
  
  ################=============== Objective
  post <- function(par) {
    -Model(par, Data)$LP
  }
  
  ################=============== Optimization
  if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS", "L-BFGS-B", "SANN", "Brent")) {
    fit <- optim(par = Initial.Values, fn = post, method = method,
                 lower = lower, upper = upper,
                 control = control, hessian = hessian)
    MAP <- fit$par
    conv <- fit$convergence
  } else if (method == "nlm") {
    fit <- suppressWarnings(nlm(p = Initial.Values, f = post, hessian = hessian))
    MAP <- fit$estimate
    conv <- ifelse(fit$code %in% c(1,2), 0, 1)
  } else if (method == "nlminb") {
    fit <- suppressWarnings(nlminb(start = Initial.Values, objective = post,
                                   lower = lower, upper = upper,
                                   control = control))
    MAP <- fit$par
    conv <- fit$convergence
  } else {
    stop("Unknown optimization method.")
  }
  
  ################=============== Convergence check
  if (conv != 0 && check.convergence) {
    stop("Optimization did not converge.")
  }
  names(MAP) <- Data[["parm.names"]]
  
  ################=============== Covariance Matrix
  if (is.null(par.cov)) {
    if (hessian) {
      if(is.null(fit$hessian)) {
        H <- optimHess(MAP, fn=post, control=control)
      } else {
        H <- fit$hessian
      }
    } else {
      g2 <- c(gradN(Model, Data, MAP, order = 2))
      H <- diag(g2)
    }
    VarCov <- tryCatch(solve(-H),
                       error = function(e) ginv(-H))
    diag(VarCov) <- abs(diag(VarCov))
  } else {
    VarCov <- par.cov
  }
  
  ################=============== Ensure PD
  if (min(eigen(VarCov)$values) <= 0 && nearPD) {
    VarCov <- Matrix::nearPD(VarCov, base.matrix = TRUE)$mat
    colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
  }
  
  ################=============== Sampling
  if (SIR) {
    cat("Start Sampling-Pareto-smoothed Importance Resampling with", iterations, "samples\n")
    samples <- samplingImportanceResampling(MAP, VarCov, Model, Data, iterations)
    posterior <- samples$posterior
    yhat      <- samples$yhat
    Monitor   <- samples$Monitor
    LP        <- samples$LP
    Dev       <- samples$Dev
    # Optional diagnostics
    weights   <- samples$weights_psis
    pareto_k  <- samples$pareto_k
    colnames(posterior) <- Data[["parm.names"]]
  } else {
    cat("Sampling Gaussian approximation with", iterations, "samples\n")
    theta <- MASS::mvrnorm(iterations, MAP, Sigma = VarCov)
    posterior <- matrix(NA, nrow=iterations, ncol=length(MAP))
    colnames(posterior) <- Data[["parm.names"]]
    tmp0 <- Model(MAP, Data)
    yhat_dim <- length(tmp0$yhat)
    mon_dim  <- length(tmp0$Monitor)
    yhat    <- matrix(NA, iterations, yhat_dim)
    Monitor <- matrix(NA, iterations, mon_dim)
    LP      <- numeric(iterations)
    Dev     <- numeric(iterations)
    for (i in 1:iterations) {
      tmp            <- Model(theta[i, ], Data)
      posterior[i, ] <- tmp$parm
      LP[i]          <- tmp$LP
      Dev[i]         <- tmp$Dev
      yhat[i, ]      <- tmp$yhat
      Monitor[i,]    <- tmp$Monitor
    }
    weights  <- rep(1/iterations, iterations)
    pareto_k <- NA
  }
  
  ################=============== Diagnostics
  stopTime <- proc.time()
  elapsedTime <- (stopTime - startTime)[3]
  # DIC
  Dbar <- mean(Dev)
  pD   <- var(Dev) / 2
  DIC  <- list(DIC = Dbar + pD, Dbar = Dbar, pD = pD)
  # ESS
  ESS <- apply(posterior, 2, coda::effectiveSize)
  
  ################=============== Output
  la.info <- list(algorithm = method, n.iter = iterations,
                  convergence = (conv == 0), elapsed.mins = elapsedTime / 60,
                  pareto_k = pareto_k, weights = weights)
  Result <- list(posterior = posterior,
                 MAP = colMeans(posterior),
                 VarCov = VarCov,
                 yhat = yhat,
                 LP = LP,
                 Monitor = Monitor,
                 Dev = Dev,
                 DIC = DIC,
                 ESS = ESS,
                 la.info = la.info)
  class(Result) <- "YABS_LA"
  cat("Done!\n")
  return(Result)
}