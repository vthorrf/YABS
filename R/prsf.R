prsf <- function (x, confidence = 0.95) {
  x <- as.mcmc.list(x)
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  
  ## Estimate mean within-chain variance (W) and between-chain variance
  ## (B/Niter), and calculate sampling variances and covariance of the
  ## estimates (varW, varB, covWB)
  ##
  ## Multivariate (upper case)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, var, simplify=TRUE), dim=c(Nvar,Nvar,Nchain))
  W <- apply(S2, c(1,2), mean)
  xbar <- matrix(sapply(x, apply, 2, mean, simplify=TRUE), nrow=Nvar,
                 ncol=Nchain)
  B <- Niter * var(t(xbar))

  if(Nvar > 1) {
      if (is.R()) {
          CW <- chol(as.matrix(nearPD(W)$mat))
          emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose=TRUE)),
                                  transpose=TRUE),
                        symmetric=TRUE, only.values=TRUE)$values[1]
      }
      else {
          emax <- eigen(qr.solve(W,B), symmetric=FALSE, only.values=TRUE)$values
      }
      mpsrf <- sqrt( (1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter )
  }
  else
    mpsrf <- NULL
  ## Univariate (lower case)
  w <- diag(W)
  b <- diag(B)

  s2 <- matrix(apply(S2, 3, diag), nrow=Nvar, ncol=Nchain)
  muhat <- apply(xbar,1,mean)
  var.w <- apply(s2, 1, var)/Nchain              
  var.b <- (2 * b^2)/(Nchain - 1)      
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) -
                              2 * muhat * var(t(s2), t(xbar)))
  
  ## Posterior interval combines all uncertainties in a t interval with
  ## center muhat, scale sqrt(V), and df.V degrees of freedom.
  V <- (Niter - 1) * w / Niter  + (1 + 1/Nchain) * b/ Niter
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * 
            var.b + 2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
  df.V <- (2 * V^2)/var.V
  
  ## Potential scale reduction factor (that would be achieved by
  ## continuing simulations forever) is estimated by 
  ##   R = sqrt(V/W) * df.adj
  ## where df.adj is a degrees of freedom adjustment for the width
  ## of the t-interval.
  ##
  ## To calculate upper confidence interval we divide R2 = R^2 into two
  ## parts, fixed and random.  The upper limit of the random part is
  ## calculated assuming that B/W has an F distribution.
  ##
  df.adj <- (df.V + 3)/(df.V + 1)
  B.df <- Nchain - 1
  W.df <- (2 * w^2)/var.w
  R2.fixed <- (Niter - 1)/Niter
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
  R2.estimate <- R2.fixed + R2.random
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * R2.random
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  
  out <- list(psrf = psrf, mpsrf=mpsrf)
  class(out) <- "gelman.diag"
  out
}
