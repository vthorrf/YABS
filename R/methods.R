#### Summary Methods
### MCMC
summary.YABS_MCMC <- function(object, ...) {
  ### Header
  cat("YABS output generated with ", object$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("Estimates based on ", object$mcmc.info$n.chains,
      if(object$mcmc.info$n.chains == 1) " chain of " else " chains of ",
      object$mcmc.info$n.iter," iterations,\n",sep="")
  cat("burn-in = ", object$mcmc.info$n.burnin, " iterations, adaptation = ",
      object$mcmc.info$n.adapt," iterations, and thin rate = ",
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
### LA
summary.YABS_LA <- function(object, ...) {
  ### Header
  cat("YABS output generated with ", object$la.info$algorithm," algorithm.\n",sep="")
  cat("Estimates based on ", object$la.info$n.iter," iterations from the joint posterior.\n",sep="")
  cat("LA ran for ",sprintf("%.3f",object$la.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(object$posterior))
  DP     <- sprintf("%.3f",apply(object$posterior,2,sd))
  HDI    <- t(apply(object$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(object$posterior),2,sign(apply(object$posterior,2,median)),"=="))
  ESS    <- object$ESS
  params <- colnames(object$posterior)
  stats <- rbind( c("", c("EAP","sd","2.5%","97.5%","overlap0","f","ESS")),
                  cbind(params, EAP, DP, t(matrix(sprintf("%.3f",t(HDI)),nrow=2)),
                        rowSums(sign(HDI)) == 0, sprintf("%.3f",f),
                        sprintf("%.2f",ESS)) )
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
        sprintf(paste("%",align[8],"s",sep=""), stats[i,8]),"\n",sep="")
  }
  if(nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',nrow(stats) - total - 1,' rows ]\n')
  }
  cat("\n")
  cat("ESS is the sample size of each posterior adjusted for autocorrelation.\n")
  cat("\n")
  cat("overlap0 checks if 0 falls in the parameter's 95% credible interval.\n")
  cat("f is the proportion of the posterior with the same sign as the mean;\n")
  cat("i.e., our confidence that the parameter is positive or negative.\n")
  cat("\n")
  cat("DIC info: (pD = var(deviance)/2).\n")
  cat("pD = ",sprintf("%.2f",object$DIC$pD)," and DIC = ",sprintf("%.2f",object$DIC$DIC),"\n",sep="")
  cat("DIC is an estimate of expected predictive error (lower is better).\n")
  invisible(stats)
}
### VB
summary.YABS_VB <- function(object, ...) {
  ### Algorithm label
  method <- object$vb.info$method
  algorithm <- switch(
    method,
    sagva   = "Stochastic Approximation for Gaussian Variational Approximation (SAGVA)",
    qnsagva = "Quasi-Newton Stochastic Approximation for Gaussian Variational Approximation (QN-SAGVA)",
    mcvi    = "Markov Chain Variational Inference (MCVI)",
    svgd    = "Stein Variational Gradient Descent (SVGD)",
    as.character(method)
  )

  ### Header
  cat("YABS output generated with ", algorithm, " algorithm.\n", sep = "")
  cat("Estimates based on ", object$vb.info$optimization.iterations,
      " optimization iterations,\n", sep = "")
  cat("yielding ", nrow(object$posterior),
      " draws from the variational approximation", sep = "")
  if (identical(method, "svgd") && is.finite(object$vb.info$n.particles)) {
    cat(" based on ", object$vb.info$n.particles, " particles", sep = "")
  }
  cat(".\n")
  cat("Minimum iterations = ", object$vb.info$min.iter,
      ", patience = ", object$vb.info$patience,
      ", and thin rate = ", object$vb.info$thinning, ".\n", sep = "")
  cat("VB ran for ", sprintf("%.3f", object$vb.info$elapsed.mins[1]),
      " minutes.\n\n", sep = "")

  ### Summary table
  posterior <- as.matrix(object$posterior)
  EAP    <- sprintf("%.3f", colMeans(posterior, na.rm = TRUE))
  DP     <- sprintf("%.3f", apply(posterior, 2, sd, na.rm = TRUE))
  HDI    <- t(apply(posterior, 2, quantile, probs = c(.025, .975),
                    na.rm = TRUE, names = FALSE))
  medsign <- sign(apply(posterior, 2, median, na.rm = TRUE))
  f      <- colMeans(sweep(sign(posterior), 2, medsign, "=="), na.rm = TRUE)
  ESS    <- object$ESS
  if (is.null(ESS) || length(ESS) != ncol(posterior)) {
    ESS <- rep(NA_real_, ncol(posterior))
  }
  params <- colnames(posterior)
  if (is.null(params)) params <- paste0("parm", seq_len(ncol(posterior)))

  stats <- rbind(
    c("", c("EAP", "sd", "2.5%", "97.5%", "overlap0", "f", "ESS")),
    cbind(params, EAP, DP,
          t(matrix(sprintf("%.3f", t(HDI)), nrow = 2)),
          rowSums(sign(HDI)) == 0,
          sprintf("%.3f", f),
          ifelse(is.finite(ESS), sprintf("%.2f", ESS), "NA"))
  )
  align <- apply(nchar(stats), 2, max) + 2
  total <- 101
  for (i in 1:min(total, nrow(stats))) {
    cat(sprintf(paste("%-", align[1], "s", sep = ""), stats[i, 1]),
        sprintf(paste("%", align[2], "s", sep = ""), stats[i, 2]),
        sprintf(paste("%", align[3], "s", sep = ""), stats[i, 3]),
        sprintf(paste("%", align[4], "s", sep = ""), stats[i, 4]),
        sprintf(paste("%", align[5], "s", sep = ""), stats[i, 5]),
        sprintf(paste("%", align[6], "s", sep = ""), stats[i, 6]),
        sprintf(paste("%", align[7], "s", sep = ""), stats[i, 7]),
        sprintf(paste("%", align[8], "s", sep = ""), stats[i, 8]), "\n", sep = "")
  }
  if (nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',
        nrow(stats) - total - 1, ' rows ]\n')
  }

  ### Convergence information
  cat("\n")
  converged <- !is.null(object$diagnostics$converged) &&
    isTRUE(object$diagnostics$converged)
  if (converged) {
    cat("The automatic stabilization threshold was reached.\n")
  } else {
    cat("The maximum iteration budget was reached before the automatic ",
        "stabilization threshold was satisfied.\n", sep = "")
    cat("This does not by itself indicate that the variational approximation ",
        "is invalid; inspect the objective and parameter stability.\n", sep = "")
  }
  if (!is.null(object$diagnostics$criterion) &&
      is.finite(object$diagnostics$criterion)) {
    cat("Final convergence criterion = ",
        format(object$diagnostics$criterion, digits = 5), ".\n", sep = "")
  }
  if (!is.null(object$diagnostics$invalid.updates) &&
      is.finite(object$diagnostics$invalid.updates) &&
      object$diagnostics$invalid.updates > 0) {
    cat("Invalid numerical updates or trajectories = ",
        object$diagnostics$invalid.updates, ".\n", sep = "")
  }
  cat("Convergence criteria are method-specific and do not constitute an exact test of posterior approximation quality.\n")

  cat("\n")
  if (identical(method, "svgd")) {
    cat("ESS is not defined for the deterministic final SVGD particle population.\n")
  } else {
    cat("ESS is the nominal number of generated draws from the fitted variational approximation.\n")
  }

  cat("\n")
  cat("overlap0 checks if 0 falls in the parameter's 95% credible interval.\n")
  cat("f is the proportion of the approximation with the same sign as the median;\n")
  cat("i.e., the estimated probability that the parameter has its predominant sign.\n")

  ### Objective information
  if (!is.null(object$vb.info$ELBO) && is.finite(object$vb.info$ELBO)) {
    cat("\n")
    cat("Estimated evidence lower bound (ELBO) = ",
        sprintf("%.3f", object$vb.info$ELBO), ".\n", sep = "")
    cat("Larger ELBO values indicate a better variational objective for fits of the same model and parameterization.\n")
  } else if (identical(method, "svgd")) {
    cat("\n")
    cat("An ELBO is not reported for SVGD because it represents the posterior with an empirical particle distribution.\n")
  }

  ### DIC information
  if (!is.null(object$DIC)) {
    cat("\n")
    cat("DIC info: (pD = var(deviance)/2).\n")
    cat("pD = ", sprintf("%.2f", object$DIC$pD),
        " and DIC = ", sprintf("%.2f", object$DIC$DIC), "\n", sep = "")
    cat("DIC is an estimate of expected predictive error (lower is better).\n")
  }

  ### Method-specific information
  if (identical(method, "qnsagva") &&
      !is.null(object$accepted.curvature.pairs)) {
    cat("Accepted quasi-Newton curvature pairs = ",
        object$accepted.curvature.pairs, ".\n", sep = "")
  }
  if (identical(method, "mcvi") && !is.null(object$hmc.step.size) &&
      is.finite(object$hmc.step.size)) {
    cat("Final HMC step size = ",
        format(object$hmc.step.size, digits = 5), ".\n", sep = "")
  }
  if (identical(method, "svgd") &&
      !is.null(object$kernel.bandwidth) &&
      is.finite(object$kernel.bandwidth)) {
    cat("Final SVGD kernel bandwidth = ",
        format(object$kernel.bandwidth, digits = 5), ".\n", sep = "")
  }

  invisible(stats)
}

#### Print Methods
### MCMC
print.YABS_MCMC <- function(x, ...) {
  ### Header
  cat("YABS output generated with ", x$mcmc.info$algorithm," algorithm.\n",sep="")
  cat("MCMC ran for ",sprintf("%.3f",x$mcmc.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(x$posterior))
  DP     <- sprintf("%.3f",apply(x$posterior,2,sd))
  HDI    <- t(apply(x$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(x$posterior),2,sign(apply(x$posterior,2,median)),"=="))
  ESS    <- x$ESS
  PSRF   <- x$PSRF
  params <- colnames(x$posterior)
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
  invisible(x)
}
### LA
print.YABS_LA <- function(x, ...) {
  ### Header
  cat("YABS output generated with ", x$la.info$algorithm," algorithm.\n",sep="")
  cat("LA ran for ",sprintf("%.3f",x$la.info$elapsed.mins[1])," minutes.\n\n",sep="")
  
  ### Summary table
  EAP    <- sprintf("%.3f",colMeans(x$posterior))
  DP     <- sprintf("%.3f",apply(x$posterior,2,sd))
  HDI    <- t(apply(x$posterior,2,quantile,c(.025,.975)))
  f      <- colMeans(sweep(sign(x$posterior),2,sign(apply(x$posterior,2,median)),"=="))
  ESS    <- x$ESS
  params <- colnames(x$posterior)
  stats <- rbind( c("", c("EAP","sd","2.5%","97.5%","overlap0","f","ESS")),
                  cbind(params, EAP, DP, t(matrix(sprintf("%.3f",t(HDI)),nrow=2)),
                        rowSums(sign(HDI)) == 0, sprintf("%.3f",f),
                        sprintf("%.2f",ESS)) )
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
        sprintf(paste("%",align[8],"s",sep=""), stats[i,8]),"\n",sep="")
  }
  if(nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',nrow(stats) - total - 1,' rows ]\n')
  }
  cat("\n")
  invisible(x)
}
### VB
print.YABS_VB <- function(x, ...) {
  ### Algorithm label
  method <- x$vb.info$method
  algorithm <- switch(
    method,
    sagva   = "Stochastic Approximation for Gaussian Variational Approximation (SAGVA)",
    qnsagva = "Quasi-Newton Stochastic Approximation for Gaussian Variational Approximation (QN-SAGVA)",
    mcvi    = "Markov Chain Variational Inference (MCVI)",
    svgd    = "Stein Variational Gradient Descent (SVGD)",
    as.character(method)
  )

  ### Header
  cat("YABS output generated with ", algorithm, " algorithm.\n", sep = "")
  cat("VB ran for ", sprintf("%.3f", x$vb.info$elapsed.mins[1]),
      " minutes.\n\n", sep = "")

  ### Summary table
  posterior <- as.matrix(x$posterior)
  EAP    <- sprintf("%.3f", colMeans(posterior, na.rm = TRUE))
  DP     <- sprintf("%.3f", apply(posterior, 2, sd, na.rm = TRUE))
  HDI    <- t(apply(posterior, 2, quantile, probs = c(.025, .975),
                    na.rm = TRUE, names = FALSE))
  medsign <- sign(apply(posterior, 2, median, na.rm = TRUE))
  f      <- colMeans(sweep(sign(posterior), 2, medsign, "=="), na.rm = TRUE)
  ESS    <- x$ESS
  if (is.null(ESS) || length(ESS) != ncol(posterior)) {
    ESS <- rep(NA_real_, ncol(posterior))
  }
  params <- colnames(posterior)
  if (is.null(params)) params <- paste0("parm", seq_len(ncol(posterior)))

  stats <- rbind(
    c("", c("EAP", "sd", "2.5%", "97.5%", "overlap0", "f", "ESS")),
    cbind(params, EAP, DP,
          t(matrix(sprintf("%.3f", t(HDI)), nrow = 2)),
          rowSums(sign(HDI)) == 0,
          sprintf("%.3f", f),
          ifelse(is.finite(ESS), sprintf("%.2f", ESS), "NA"))
  )
  align <- apply(nchar(stats), 2, max) + 2
  total <- 101
  for (i in 1:min(total, nrow(stats))) {
    cat(sprintf(paste("%-", align[1], "s", sep = ""), stats[i, 1]),
        sprintf(paste("%", align[2], "s", sep = ""), stats[i, 2]),
        sprintf(paste("%", align[3], "s", sep = ""), stats[i, 3]),
        sprintf(paste("%", align[4], "s", sep = ""), stats[i, 4]),
        sprintf(paste("%", align[5], "s", sep = ""), stats[i, 5]),
        sprintf(paste("%", align[6], "s", sep = ""), stats[i, 6]),
        sprintf(paste("%", align[7], "s", sep = ""), stats[i, 7]),
        sprintf(paste("%", align[8], "s", sep = ""), stats[i, 8]), "\n", sep = "")
  }
  if (nrow(stats) > total) {
    cat(' [ reached getOption("max.print") -- omitted ',
        nrow(stats) - total - 1, ' rows ]\n')
  }

  ### Convergence information
  cat("\n")
  converged <- !is.null(x$diagnostics$converged) &&
    isTRUE(x$diagnostics$converged)
  if (converged) {
    cat("The automatic stabilization threshold was reached.\n")
  } else {
    cat("The maximum iteration budget was reached before the automatic ",
        "stabilization threshold was satisfied.\n", sep = "")
    cat("This does not by itself indicate that the variational approximation ",
        "is invalid; inspect the objective and parameter stability.\n", sep = "")
  }
  if (!is.null(x$diagnostics$criterion) && is.finite(x$diagnostics$criterion)) {
    cat("Final convergence criterion = ",
        format(x$diagnostics$criterion, digits = 5), ".\n", sep = "")
  }
  if (!is.null(x$diagnostics$invalid.updates) &&
      is.finite(x$diagnostics$invalid.updates) &&
      x$diagnostics$invalid.updates > 0) {
    cat("Invalid numerical updates or trajectories = ",
        x$diagnostics$invalid.updates, ".\n", sep = "")
  }
  cat("\n")

  invisible(x)
}
