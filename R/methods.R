#### Summary Method
summary.YABS <- function(object, ...) {
  ### Header
  cat("YABS output generated with ", object$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("Estimates based on 1 chain of ", object$mcmc.info$n.iter," iterations,\n",sep="")
  cat("burn-in = ", object$mcmc.info$n.burnin, " iterations and thin rate = ",
      object$mcmc.info$n.thin,",\n",sep="")
  cat("yielding ",nrow(object$posterior)," total samples from the joint posterior.\n",sep="")
  cat("MCMC ran for ",sprintf("%.3f",object$mcmc.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(object$posterior))
  DP     <- sprintf("%.3f",apply(object$posterior,2,sd))
  HDI    <- t(apply(object$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(object$posterior),2,sign(apply(object$posterior,2,median)),"=="))
  ESS    <- object$ESS
  PSRF   <- object$PSRF
  params <- colnames(object$posterior)
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
  cat("pD = ",sprintf("%.2f",object$DIC$pD)," and DIC = ",sprintf("%.2f",object$DIC$DIC),"\n",sep="")
  cat("DIC is an estimate of expected predictive error (lower is better).\n")
  if({length(object$MPSRF) != 0} | !is.na(object$MPSRF)) {
    cat("MPSRF = ",sprintf("%.3f",object$MPSRF),"\n",sep="")
    cat("MPSRF is the multivariate potential scale reduction factor (at convergence, MPSRF=1).\n")
  }
  invisible(stats)
}

#### Print Method
print.YABS <- function(object, ...) {
  ### Header
  cat("YABS output generated with ", object$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("MCMC ran for ",sprintf("%.3f",object$mcmc.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(object$posterior))
  DP     <- sprintf("%.3f",apply(object$posterior,2,sd))
  HDI    <- t(apply(object$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(object$posterior),2,sign(apply(object$posterior,2,median)),"=="))
  ESS    <- object$ESS
  PSRF   <- object$PSRF
  params <- colnames(object$posterior)
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
  invisible(object)
}