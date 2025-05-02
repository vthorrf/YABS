\name{WBIC}
\alias{WBIC}
\title{Widely Applicable Bayesian Information Criterion}
\description{
This function calculates an approximation of the Widely Applicable Bayesian Information Criterion (WBIC; Watanabe, 2013).
}
\usage{
WBIC(x)
}
\arguments{
   \item{x}{This required argument accepts a 𝑆×𝑁 matrix of pointwise log-likelihood (LL) for N records and S samples.}
}
\value{
The Widely Applicable Bayesian Information Criterion.
}
\references{
Watanabe, S. (2013). A widely applicable Bayesian information criterion. The Journal of Machine Learning Research, 14(1), 867-897. https://dl.acm.org/doi/abs/10.5555/2567709.2502609
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
y <- rnorm(N, cbind(1,X) \%*\% betas, 1)

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
  mu <- cbind(1,Data$X) \%*\% beta
  pwLL <- dnorm(Data$y, mu, sigma, log=TRUE)
  LL <- sum(pwLL)
  
  ### Log-Posterior
  LP <- LL + Lp
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=pwLL, yhat=rnorm(Data$n, mu, sigma), parm=parm)
  return(Modelout)
}
Model <- compiler::cmpfun(Model)

### FIT====
fit <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
            burnin=100, status=110, thinning=1, algo="mwg")
WBIC(fit$Monitor)
}
