\name{infoCrit}
\alias{infoCrit}
\title{Information Criteria for Model Evaluation}
\description{
This function computes a set of classical and complexity-based information criteria for model evaluation using deviance and posterior samples.
}
\usage{
infoCrit(posterior, deviance, sampleSize, infoC = TRUE)
}
\arguments{
  \item{posterior}{This required argument accepts a \eqn{S \times D} matrix of posterior estimates for \eqn{D} parameters and \eqn{S} samples.}
  \item{deviance}{This required argument accepts a vector of length \eqn{S} of the estimates of deviance.}
  \item{sampleSize}{The sample size \eqn{n} used to fit the model.}
  \item{infoC}{Logical. If \code{TRUE} (default), additional complexity-based criteria are computed using the posterior covariance matrix. If \code{FALSE}, only classical information criteria are returned.}
}
\value{
A data frame containing the computed information criteria:
  \describe{
    \item{\code{aic}}{Akaike Information Criterion.}
    \item{\code{aicc}}{Finite-sample corrected AIC.}
    \item{\code{bc}}{Bridge Criterion.}
    \item{\code{bic}}{Bayesian Information Criterion (Schwarz).}
    \item{\code{hqc}}{Hannan–Quinn Criterion.}
    \item{\code{fcomp}}{Frobenius-norm complexity (if \code{infoC = TRUE}).}
    \item{\code{icomp}}{Information complexity (if \code{infoC = TRUE}).}
    \item{\code{ifim}}{Cramér–Rao complexity (if \code{infoC = TRUE}).}
    \item{\code{peu}}{Posterior expected utility (if \code{infoC = TRUE}).}
    \item{\code{peuln}}{Log-scaled sample size PEU (if \code{infoC = TRUE}).}
    \item{\code{peums}}{PEU under model misspecification (if \code{infoC = TRUE}).}
    \item{\code{mdl}}{Minimum description length (if \code{infoC = TRUE}).}
    \item{\code{sc2}}{Stochastic complexity including covariance determinant (if \code{infoC = TRUE}).}
    \item{\code{sc3}}{Simplified stochastic complexity (if \code{infoC = TRUE}).}
  }
}
\references{
Akaike, H. (1974). A new look at the statistical model identification. IEEE Transactions on Automatic Control, 19(6), 716–723.

Bozdogan, H. (1990). On the information-based measure of covariance complexity and its application to the evaluation of multivariate linear models. Communications in Statistics - Theory and Methods, 19(1), 221–278.

Hannan, E. J., & Quinn, B. G. (1979). The determination of the order of an autoregression. Journal of the Royal Statistical Society: Series B, 41(2), 190–195.

Rissanen, J. (1986). Stochastic complexity and modeling. The Annals of Statistics, 14(3), 1080-1100.

Schwarz, G. (1978). Estimating the dimension of a model. The Annals of Statistics, 6(2), 461–464.
}
\examples{##### Multiple Linear Regression with 5 Predictors
### Packages and functions====
require(YABS)
require(compiler)

### Random data====
seed <- 1234; N <- 200; V <- 5
set.seed(seed)
X <- sapply(1:V, function(g) rnorm(N))
betas <- runif(V+1, .3, .7)
y <- rnorm(N, cbind(1,X) %*% betas, 1)

### DATA LIST====
parm.names <- c(paste("beta",0:V,sep=""),"sigma")
pos.beta   <- grep("beta", parm.names)
pos.sigma  <- grep("sigma", parm.names)
PGF <- function(Data) {
  beta  <- rnorm(ncol(Data$X)+1)
  sigma <- rnorm(1)
  return(c(beta, sigma))
}
Data <- list( X=X, y=y, n=N, parm.names=parm.names, pos.beta=pos.beta,
              pos.sigma=pos.sigma, PGF=PGF, mon.names="LP" )
Initial.Values <- PGF(Data)

### MODEL====
Model <- function(parm, Data){
  ### Parameters
  beta  <- parm[Data$pos.beta]
  sigma <- exp(parm[Data$pos.sigma])
  
  ### Log-Prior
  beta.prior  <- sum(dnorm(beta, 0, 1, log=TRUE))
  sigma.prior <- sum(dgamma(sigma, 1e-2, 1e-2, log=TRUE))
  Lp <- beta.prior + sigma.prior
  
  ### Log-Likelihood
  mu <- cbind(1,Data$X) %*% beta
  pwLL <- dnorm(Data$y, mu, sigma, log=TRUE)
  LL <- sum(pwLL)
  
  ### Log-Posterior
  LP <- LL + Lp
  
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=pwLL, yhat=rnorm(Data$n, mu, sigma), parm=c(beta, sigma))
  return(Modelout)
}
Model <- compiler::cmpfun(Model)

### FIT====
fit <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
            burnin=100, status=110, thinning=1, algo="mwg")

### Compute information criteria====
infoCrit(fit$posterior, fit$Dev, sampleSize = N)
}