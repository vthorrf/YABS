MCMC <- function(Model, Data, Initial.Values=NULL, iterations=NULL, burnin=NULL,
                 status=NULL, thinning=NULL, adapt=NULL, nchains=1, parallel=FALSE,
                 cores=NULL, update.progress=NULL, opt.init=TRUE, par.cov=NULL,
                 algo=c("mwg","rwm","barker","ohss")) {
  ################=============== Initial settings
  ## Default values for the arguments and some error handling
  if(length(algo) != 1)        algo            <- NULL
  if(is.null(algo))            algo            <- "rwm"
  if(is.null(Initial.Values))  Initial.Values  <- Data$PGF(Data)
  if(is.null(iterations))      iterations      <- 100
  if(iterations <= 0)          iterations      <- 100
  if(is.null(burnin))          burnin          <- 0
  if(burnin < 0)               burnin          <- 0
  if(is.null(status))          status          <- round(iterations/10,0)
  if(status <= 0)              status          <- round(iterations/10,0)
  if(status == 0)              status          <- 1
  if(is.null(adapt))           adapt           <- 0
  if(adapt < 0)                adapt           <- 0
  if(is.null(cores))           cores           <- 2
  if(cores <= 0)               cores           <- 2
  if(is.null(update.progress)) update.progress <- 2
  if(update.progress <= 0)     update.progress <- 2
  if(!is.matrix(par.cov) & !is.null(par.cov)) stop("'par.cov' should be a matrix. Please check the documentation.")
  status <- round(status,0)
  h <- 1e-6; ITER <- iterations + burnin + adapt
  if(is.null(thinning))        thinning        <- 1
  if(thinning <= 0)            thinning        <- 1
  
  ## Improve initial values and estimate step sizes
  post <- function(par) return( suppressWarnings(Model(par, Data)$LP) )
  if(opt.init) {
    if(!is.null(par.cov)) {
      cat("Improving initial values with MAP estimation.\n")
    } else {
      cat("Improving initial values and finding initial step sizes with MAP estimation.\n")
    }
    if(!is.null(par.cov)) Hessian <- FALSE else Hessian <- TRUE
    MAP <- tryCatch(optim( Initial.Values, post, method="L-BFGS",
                           control=list(fnscale=-1), hessian=Hessian ),
                    error=function(e) {
                      optim( Initial.Values, post, method="BFGS",
                             control=list(fnscale=-1), hessian=Hessian )
                    })
    if(!is.null(par.cov)) epsilon <- par.cov else epsilon <- tryCatch( solve(-MAP$hessian),
                                                                       error=function(e) ginv(-MAP$hessian) )
    if(min(eigen(epsilon)$values) < 0) {
      epsilon <- tryCatch( nearPD(epsilon, base.matrix=T)$mat, 
                           error=function(e) epsilon )
    }
    Initial.Values <- mvrnorm(ncol(epsilon)+4+nchains, MAP$par, Sigma=epsilon, empirical=TRUE)
    if(MAP$convergence == 0) {
      adapt <- 0
      cat("Initial optimization was sufficient for estimating the step sizes.\n")
    }
  } else {
    if(is.null(par.cov)) {
      epsilon <- diag( length(Initial.Values) ) * .1
    } else { epsilon <- par.cov }
    if(min(eigen(epsilon)$values) < 0) {
      epsilon <- tryCatch( nearPD(epsilon, base.matrix=T)$mat, 
                           error=function(e) epsilon )
    }
    Initial.Values <- mvrnorm(ncol(epsilon)+4+nchains, Initial.Values, Sigma=epsilon, empirical=TRUE)
  }
  lt_epsilon <- t(chol(epsilon + diag(1e-4, nrow=ncol(epsilon), ncol=ncol(epsilon))))
  
  ## Prepare the parameters for the MCMC algorithms
  liv        <- length(Initial.Values[1,])
  acceptance <- 0
  MO0        <- Model(Initial.Values[1,], Data)
  thinned    <- matrix(Initial.Values[1,], floor(ITER/thinning)+1,
                       length(Initial.Values[1,]), byrow=TRUE)
  postpred   <- matrix(MO0$yhat, floor(ITER/thinning)+1,
                       length(MO0$yhat), byrow=TRUE)
  DEV  <- matrix(MO0[["Dev"]], floor(ITER/thinning)+1, 1)
  MON  <- matrix(MO0[["Monitor"]], floor(ITER/thinning)+1,
                 length(MO0[["Monitor"]]), byrow=TRUE)
  BURN <- floor({burnin + adapt}/thinning)+1
  
  ################=============== Fit model
  if(algo == "mwg") {
    ##############=============== Metropolis-within-Gibbs
    method = "MWG"
    cat("Algorithm: Metropolis-within-Gibbs\n\n")
    startTime = proc.time()
    if(parallel) {
      cat("Running ", nchains," chains in parallel\n", sep="")
      cl <- makeCluster(cores)
      pboptions(nout=update.progress)
      fits <- pblapply(X=1:nchains, function(i) {
        temp0 <- Model(Initial.Values[i,], Data)
        harmwg(Model, Data, ITER, status, thinning, acceptance,
               DEV, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      }, cl = cl)
      stopCluster(cl)
    } else {
      fits <- lapply(1:nchains, function(i) {
        cat("=========Chain number ", i,"=========\n", sep="")
        temp0 <- Model(Initial.Values[i,], Data)
        harmwg(Model, Data, ITER, status, thinning, acceptance,
               DEV, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      })
    }
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else if(algo == "rwm") {
    ##############=============== Random-walk Metropolis
    method = "RWM"
    cat("Algorithm: Random-walk Metropolis\n\n")
    startTime = proc.time()
    if(parallel) {
      cat("Running ", nchains," chains in parallel\n", sep="")
      cl <- makeCluster(cores)
      pboptions(nout=update.progress)
      fits <- pblapply(X=1:nchains, function(i) {
        temp0 <- Model(Initial.Values[i,], Data)
        harm(Model, Data, ITER, status, thinning, acceptance,
             DEV, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      }, cl = cl)
      stopCluster(cl)
    } else {
      fits <- lapply(1:nchains, function(i) {
        cat("=========Chain number ", i,"=========\n", sep="")
        temp0 <- Model(Initial.Values[i,], Data)
        harm(Model, Data, ITER, status, thinning, acceptance,
             DEV, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      })
    }
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else if(algo == "barker") {
    ##############=============== Barker Proposal Metropolis
    method = "BPM"
    cat("Algorithm: Barker Proposal Metropolis\n\n")
    startTime = proc.time()
    if(parallel) {
      cat("Running ", nchains," chains in parallel\n", sep="")
      cl <- makeCluster(cores)
      pboptions(nout=update.progress)
      fits <- pblapply(X=1:nchains, function(i) {
        temp0 <- Model(Initial.Values[i,], Data)
        gcharm(Model, Data, ITER, status, thinning, acceptance,
               DEV, h, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      }, cl = cl)
      stopCluster(cl)
    } else {
      fits <- lapply(1:nchains, function(i) {
        cat("=========Chain number ", i,"=========\n", sep="")
        temp0 <- Model(Initial.Values[i,], Data)
        gcharm(Model, Data, ITER, status, thinning, acceptance,
               DEV, h, liv, MON, temp0, thinned, postpred, adapt, lt_epsilon)
      })
    }
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else if(algo == "ohss") {
    ##############=============== Oblique Hyperrectangle Slice Sampler
    method = "OHSS"
    cat("Algorithm: Oblique Hyperrectangle Slice Sampler\n\n")
    startTime = proc.time()
    if(parallel) {
      cat("Running ", nchains," chains in parallel\n", sep="")
      cl <- makeCluster(cores)
      pboptions(nout=update.progress)
      fits <- pblapply(X=1:nchains, function(i) {
        temp0 <- Model(Initial.Values[i,], Data)
        ohss(Model, Data, ITER, status, thinning, acceptance,
             DEV, liv, MON, temp0, thinned, postpred, adapt)
      }, cl = cl)
      stopCluster(cl)
    } else {
      fits <- lapply(1:nchains, function(i) {
        cat("=========Chain number ", i,"=========\n", sep="")
        temp0 <- Model(Initial.Values[i,], Data)
        ohss(Model, Data, ITER, status, thinning, acceptance,
             DEV, liv, MON, temp0, thinned, postpred, adapt)
      })
    }
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else  stop("Unkown MCMC algorithm. Please, check documentation.")
  
  ################=============== Results
  ## Posterior list
  post_list <- lapply(fits, function(g) {
    temp0 <- as.matrix(g$thinned[{1:nrow(g$thinned)} %!in% seq_len(BURN),])
    dev <- unlist(g$Dev[{1:nrow(g$Dev)} %!in% seq_len(BURN),])
    temp1 <- cbind(temp0, dev, {{2 * liv} + dev})
    colnames(temp1) <- c(Data$parm.names, "deviance", "aic")
    return(temp1)
  })
  ## Posterior predictive
  ppred <- lapply(fits, function(g) {
    as.matrix(g$postpred[{1:nrow(g$postpred)} %!in% seq_len(BURN),])
  })
  ## Additionally monitored variables/parameters
  Monitor <- lapply(fits, function(g) {
    as.matrix(g$Mon[{1:nrow(g$Mon)} %!in% seq_len(BURN),])
  })
  ## LogPosterior
  logPost <- lapply(post_list, function(g) {
    apply(g, 1, post)
  })
  ## Deviance
  deviance <- lapply(post_list, function(g) {
    g[,"deviance"]
  })
  ## Posterior data frame
  posterior <- data.frame(do.call("rbind",post_list))
  ## Calculate DIC
  Dbar      <- mean(unlist(deviance))
  pD        <- var(unlist(deviance))/2
  DIC       <- list(DIC=Dbar + pD, Dbar=Dbar, pD=pD)
  ## Acceptance rate
  accRate   <- mean(sapply(fits, function(g) g$Acc))
  ## Effective Sample Size
  ESS       <- apply(posterior,2,effectiveSize)
  ## R_hat (Potential scale reduction factor)
  if( nchains == 1 ) {
    halves    <- suppressWarnings( apply(posterior, 2, split,
                                         f=c(rep(1, floor(nrow(posterior)/2)),
                                         rep(2, nrow(posterior)-floor(nrow(posterior)/2))),
                                         drop=T) )
    GD        <- prsf(lapply(1:2, function(R) {
                    temp <- sapply(1:ncol(posterior), function(g) halves[[g]][[R]])
                    temp <- mcmc(as.matrix(temp[1:min(lengths(halves[[1]])),]))
                    colnames(temp) <- colnames(posterior)
                    return(temp)
                 }))
    
  } else {
    GD <- prsf(lapply(post_list, as.mcmc))
  }
  ## Information regarding the run
  mcmc.info <- list(algorithm=method, n.iter=nrow(post_list[[1]]),
                    n.burnin=burnin, n.thin=thinning, n.adapt=adapt,
                    n.chains=nchains, elapsed.mins=elapsedTime/60)
  ## Final list of results
  Result    <- list(posterior=posterior,
                    yhat=do.call("rbind",ppred),
                    LP=unlist(logPost),
                    Monitor=if(ncol(Monitor[[1]]) == 1) {
                      unlist(Monitor)
                    } else { do.call("rbind",Monitor) },
                    DIC=DIC,
                    acc=accRate,
                    ESS=ESS,
                    PSRF=GD$psrf[,1],
                    MPSRF=GD$mpsrf,
                    mcmc.info=mcmc.info)
  class(Result) <- "YABS"
  return(Result)
}

'%!in%' <- function(x,y){ !('%in%'(x,y)) }
