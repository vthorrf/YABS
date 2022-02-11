MCMC <- function(Model, Data, Initial.Values=NULL, iterations=NULL,
                 burnin=NULL, status=NULL, thinning=NULL,
                 algo=c("harmwg","harm","sharm")) {
  ### Initial settings
  if(length(algo) > 1)        algo           <- "harmwg"
  if(is.null(Initial.Values)) Initial.Values <- Data$PGF(Data)
  if(is.null(iterations))     iterations     <- 100
  if(iterations <= 0)         iterations     <- 100
  if(is.null(burnin))         burnin         <- 0
  if(burnin < 0)              burnin         <- 0
  if(is.null(status))         status         <- round(iterations/10,0)
  if(status <= 0)             status         <- round(iterations/10,0)
  if(status == 0)             status         <- 1
  status <- round(status,0)
  h <- 1e-6; ITER   <- iterations + burnin
  if(is.null(thinning))       thinning       <- 1
  if(thinning <= 0)           thinning       <- 1
  liv        <- length(Initial.Values)
  acceptance <- 0#rep(0, liv)
  MO0        <- Model(Initial.Values, Data)
  thinned    <- matrix(Initial.Values, floor(ITER/thinning)+1,
                       length(Initial.Values), byrow=TRUE)
  DEV  <- matrix(MO0[["Dev"]], floor(ITER/thinning)+1, 1)
  MON  <- matrix(MO0[["Monitor"]], floor(ITER/thinning)+1,
                 length(MO0[["Monitor"]]), byrow=TRUE)
  BURN <- floor(burnin/thinning)+1
  
  ### Fit model
  if(algo == "harmwg") {
    ##############=============== Hit-and-Run Metropolis-within-Gibbs
    method = "HARM-WG"
    cat("Algorithm: Hit-and-Run Metropolis-within-Gibbs\n\n")
    startTime = proc.time()
    fit <- harmwg(Model, Data, ITER, status, thinning, acceptance,
                  DEV, liv, MON, MO0, thinned)
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else if(algo == "harm") {
    ##############=============== Hit-and-Run Metropolis
    method = "HARM"
    cat("Algorithm: Hit-and-Run Metropolis\n\n")
    startTime = proc.time()
    fit <- harm(Model, Data, ITER, status, thinning, acceptance,
                DEV, liv, MON, MO0, thinned)
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else if(algo == "sharm") {
    ##############=============== Steepest Hit-and-Run Metropolis
    method = "SHARM"
    cat("Algorithm: Steepest Hit-and-Run Metropolis\n\n")
    startTime = proc.time()
    fit <- sharm(Model, Data, ITER, status, thinning, acceptance,
                 DEV, h, liv, MON, MO0, thinned)
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("\n")
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish.\n", sep="")
  } else stop("Unkown MCMC algorithm. Please, check documentation.")
  
  ### Results
  post      <- as.matrix(fit$thinned[{1:nrow(fit$thinned)} %!in% seq_len(BURN),])
  colnames(post) <- Data$parm.names
  logPost   <- fit$Mon[{1:length(fit$Mon)} %!in% seq_len(BURN)]
  deviance  <- fit$Dev[{1:length(fit$Dev)} %!in% seq_len(BURN)]
  aic       <- {{2 * liv} + deviance}
  posterior <- data.frame(post, "deviance"=deviance, "aic"=aic)
  Dbar      <- mean(deviance)
  pD        <- var(deviance)/2
  DIC       <- list(DIC=Dbar + pD, Dbar=Dbar, pD=pD)
  accRate   <- fit$Acceptance / ITER
  ESS       <- apply(posterior,2,effectiveSize)
  halves    <- suppressWarnings( apply(posterior, 2, split, c(1,2), drop=T) )
  GD        <- prsf(lapply(1:2, function(R) {
                  temp <- sapply(1:ncol(posterior), function(g) halves[[g]][[R]])
                  temp <- mcmc(as.matrix(temp[1:min(lengths(halves[[1]])),]))
                  colnames(temp) <- c(Data$parm.names, "deviance", "aic")
                  return(temp)
               }))
  mcmc.info <- list(algorithm=method, n.iter=ITER, n.burnin=burnin,
                    n.thin=thinning, elapsed.mins=elapsedTime/60)
  Result    <- list(posterior=posterior, LP=logPost, DIC=DIC,
                    acc=accRate, ESS=ESS, PSRF=GD$psrf[,1],
                    MPSRF=GD$mpsrf, mcmc.info=mcmc.info)
  class(Result) <- "YABS"
  return(Result)
}

'%!in%' <- function(x,y){ !('%in%'(x,y)) }

#### Summary Method
summary.YABS <- function(oop) {
  ### Header
  cat("YABS output generated with ", oop$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("Estimates based on 1 chain of ", oop$mcmc.info$n.iter," iterations,\n",sep="")
  cat("burn-in = ", oop$mcmc.info$n.burnin, " iterations and thin rate = ",
      oop$mcmc.info$n.thin,",\n",sep="")
  cat("yielding ",nrow(oop$posterior)," total samples from the joint posterior.\n",sep="")
  cat("MCMC ran for ",sprintf("%.3f",oop$mcmc.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(oop$posterior))
  DP     <- sprintf("%.3f",apply(oop$posterior,2,sd))
  HDI    <- t(apply(oop$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(oop$posterior),2,sign(apply(oop$posterior,2,median)),"=="))
  ESS    <- oop$ESS
  PSRF   <- oop$PSRF
  params <- colnames(oop$posterior)
  stats <- rbind( c("", c("EAP","sd","2.5%","97.5%","overlap0","f","ESS","PSRF")),
                  cbind(params, EAP, DP, t(matrix(sprintf("%.3f",t(HDI)),nrow=2)),
                        rowSums(sign(HDI)) == 0, sprintf("%.3f",f),
                        sprintf("%.2f",ESS), sprintf("%.3f",PSRF)) )
  align <- apply(nchar(stats), 2, max) + 2
  total <- 101
  for(i in 1:min(total,nrow(stats))) {
    cat(sprintf(paste("%-",align[1],"s",sep=""), stats[i,1]),
        sprintf(paste("%",align[2],"s",sep=""), stats[i,2]),
        sprintf(paste("%",align[3],"s",sep=""), stats[i,3]),
        sprintf(paste("%",align[4],"s",sep=""), stats[i,4]),
        sprintf(paste("%",align[5],"s",sep=""), stats[i,5]),
        sprintf(paste("%",align[6],"s",sep=""), stats[i,6]),
        sprintf(paste("%",align[7],"s",sep=""), stats[i,7]),
        sprintf(paste("%",align[8],"s",sep=""), stats[i,8]),
        sprintf(paste("%",align[9],"s",sep=""), stats[i,9]),"\n",sep="")
  }
  if(nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',nrow(stats) - total - 1,' rows ]\n')
  }
  cat("\n")
  if(sum(PSRF >= 1.1) == 0) {
    cat("Successful convergence based on PSRF (or Rhat) values (all < 1.1).\n")
  } else { cat("**WARNING** PSRF (or Rhat) values indicate convergence failure.\n") }
  cat("PSRF is the potential scale reduction factor (at convergence, PSRF=1).\n")
  cat("ESS is the sample size of each posterior adjusted for autocorrelation.\n")
  cat("\n")
  cat("overlap0 checks if 0 falls in the parameter's 95% credible interval.\n")
  cat("f is the proportion of the posterior with the same sign as the mean;\n")
  cat("i.e., our confidence that the parameter is positive or negative.\n")
  cat("\n")
  cat("DIC info: (pD = var(deviance)/2).\n")
  cat("pD = ",sprintf("%.2f",oop$DIC$pD)," and DIC = ",sprintf("%.2f",oop$DIC$DIC),"\n",sep="")
  cat("DIC is an estimate of expected predictive error (lower is better).\n")
  if({length(oop$MPSRF) != 0} | !is.na(oop$MPSRF)) {
    cat("MPSRF = ",sprintf("%.3f",oop$MPSRF),"\n",sep="")
    cat("MPSRF is the multivariate potential scale reduction factor (at convergence, MPSRF=1).\n")
  }
}

#### Print Method
print.YABS <- function(oop) {
  ### Header
  cat("YABS output generated with ", oop$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("MCMC ran for ",sprintf("%.3f",oop$mcmc.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(oop$posterior))
  DP     <- sprintf("%.3f",apply(oop$posterior,2,sd))
  HDI    <- t(apply(oop$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(oop$posterior),2,sign(apply(oop$posterior,2,median)),"=="))
  ESS    <- oop$ESS
  PSRF   <- oop$PSRF
  params <- colnames(oop$posterior)
  stats <- rbind( c("", c("EAP","sd","2.5%","97.5%","overlap0","f","ESS","PSRF")),
                  cbind(params, EAP, DP, t(matrix(sprintf("%.3f",t(HDI)),nrow=2)),
                        rowSums(sign(HDI)) == 0, sprintf("%.3f",f),
                        sprintf("%.2f",ESS), sprintf("%.3f",PSRF)) )
  align <- apply(nchar(stats), 2, max) + 2
  total <- 101
  for(i in 1:min(total,nrow(stats))) {
    cat(sprintf(paste("%-",align[1],"s",sep=""), stats[i,1]),
        sprintf(paste("%",align[2],"s",sep=""), stats[i,2]),
        sprintf(paste("%",align[3],"s",sep=""), stats[i,3]),
        sprintf(paste("%",align[4],"s",sep=""), stats[i,4]),
        sprintf(paste("%",align[5],"s",sep=""), stats[i,5]),
        sprintf(paste("%",align[6],"s",sep=""), stats[i,6]),
        sprintf(paste("%",align[7],"s",sep=""), stats[i,7]),
        sprintf(paste("%",align[8],"s",sep=""), stats[i,8]),
        sprintf(paste("%",align[9],"s",sep=""), stats[i,9]),"\n",sep="")
  }
  if(nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',nrow(stats) - total - 1,' rows ]\n')
  }
  cat("\n")
  if(sum(PSRF >= 1.1) == 0) {
    cat("Successful convergence based on PSRF (or Rhat) values (all < 1.1).\n")
  } else { cat("**WARNING** PSRF (or Rhat) values indicate convergence failure.\n") }
  cat("\n")
}
