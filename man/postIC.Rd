\name{postIC}
\alias{postIC}
\title{Fully Bayesian Information Criterion}
\description{
This function computes multiple Bayesian model evaluation criteria.
}
\usage{
postIC(x, posterior = NULL)
}
\arguments{
   \item{x}{A required \eqn{S \times N} matrix of pointwise log-likelihood values, where \eqn{S} is the number of posterior samples and \eqn{N} is the number of observations.}
   \item{posterior}{An optional \eqn{S \times d} matrix of posterior samples for the model parameters. This argument is required for the computation of the Bayesian MDL.}
}
\value{
A list containing the following elements:
\itemize{
  \item{\code{WBIC}}{The Widely Applicable Bayesian Information Criterion.}
  \item{\code{WAIC}}{The Widely Applicable Information Criterion.}
  \item{\code{p_WAIC}}{The effective number of parameters estimated via WAIC.}
  \item{\code{lppd}}{The log pointwise predictive density.}
  \item{\code{MDL}}{The fully Bayesian Minimum Description Length (returned on the deviance scale). Returns \code{NA} if \code{posterior} is not provided.}
}
}
\details{
WBIC is defined as a weighted average of the log-likelihood evaluated at an inverse temperature of \eqn{\beta = 1 / \log(n)}, where \eqn{n} is the sample size.

WAIC is computed as:
\deqn{WAIC = -2 (lppd - p_{WAIC})}
where \eqn{lppd = \sum_i \log \left( \frac{1}{S} \sum_s \exp(\ell_{si}) \right)} and \eqn{p_{WAIC}} is the sum of the posterior variances of the pointwise log-likelihood.

The Bayesian MDL is computed using a Laplace-type approximation:
\deqn{
MDL \approx -2 \left[
\mathbb{E}_{post}[\log f(y|\theta)]
- \frac{d}{2} \log n
+ \frac{1}{2} \log |\Sigma|
+ \frac{d}{2} \log(2\pi e)
\right]
}
where \eqn{d} is the number of parameters and \eqn{\Sigma} is the posterior covariance matrix.

All criteria are returned on the deviance scale, so smaller values indicate better model fit.
}
\references{
Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571–3594.

Watanabe, S. (2013). A widely applicable Bayesian information criterion. Journal of Machine Learning Research, 14(1), 867–897.

Rissanen, J. (1989). Stochastic complexity in statistical inquiry. World Scientific.

Grünwald, P. D. (2007). The Minimum Description Length Principle. MIT Press.
}
\examples{
##### Multiple Linear Regression with 5 Predictors
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

Data <- list(X=X, y=y, n=N, parm.names=parm.names,
             pos.beta=pos.beta, pos.sigma=pos.sigma,
             PGF=PGF, mon.names="LP")

Initial.Values <- PGF(Data)

### MODEL====
Model <- function(parm, Data){
  beta  <- parm[Data$pos.beta]
  sigma <- exp(parm[Data$pos.sigma])
  
  beta.prior  <- sum(dnorm(beta, 0, 1, log=TRUE))
  sigma.prior <- sum(dgamma(sigma, 1e-2, 1e-2, log=TRUE))
  Lp <- beta.prior + sigma.prior
  
  mu <- cbind(1,Data$X) \%*\% beta
  pwLL <- dnorm(Data$y, mu, sigma, log=TRUE)
  LL <- sum(pwLL)
  
  LP <- LL + Lp
  
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=pwLL, yhat=rnorm(Data$n, mu, sigma), parm=c(beta, sigma))
  return(Modelout)
}

Model <- compiler::cmpfun(Model)

### FIT====
fit <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
            burnin=100, status=110, thinning=1, algo="mwg")

### Criteria====
postIC(fit$Monitor, posterior = fit$posterior)
}