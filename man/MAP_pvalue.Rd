\name{MAP_pvalue}
\alias{MAP_pvalue}
\title{MAP-based p-value}
\description{
The MAP-based p-value is defined as the ratio of the posterior density of a reference value to the Maximum A Posteriori (MAP) estimate of the posterior density. Therefore, it indicates how much more likely is the probability of the reference value under the posterior in comparison to the MAP estimate under the same posterior.
}
\usage{
MAP_pvalue <- function(x, reference.value=0, MAP.ID=NULL)
}
\arguments{
   \item{x}{A vector of the posterior samples of a specific parameter.}
   \item{reference.value}{The reference value to be tested. Defaults to 0.}
   \item{MAP.ID}{The position of the MAP estimate in x. If not provided, the EAP is used instead.}
}
\value{
The ratio of the posterior density of a reference value to the posterior density of the MAP estimate.
}
\references{
Makowski, D., Ben-Shachar, M. S., Chen, S. A., & Lüdecke, D. (2019). Indices of effect existence and significance in the Bayesian framework. Frontiers in Psychology, 10, 2767. https://doi.org/10.3389/fpsyg.2019.02767
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
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  
  ### Log-Posterior
  LP <- LL + Lp
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(Data$n, mu, sigma), parm=parm)
  return(Modelout)
}
Model <- compiler::cmpfun(Model)

### FIT====
fit <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
            burnin=100, status=110, thinning=1, algo="mwg")
MAP_pvalue(fit$posterior[,1], reference.value=0)
}
