\name{prsf}
\alias{prsf}
\title{Gelman and Rubin's convergence diagnostic}
\description{
Internal function. A copy of the `gelman.diag` function from the package coda. The main differences are: (i) it will compute the nearest positive definite matrix of the Cholesky decomposition of the posterior, whenever it fails to be PD; (ii) there is no transformation of the parameters; (iii) there is no autoburnin; and (iv) the multivariate PRSF will be calculate whenever there is more than one parameter. See the `gelman.diag` function from the package coda for more details.
}
\usage{
prsf(x, confidence = 0.95)
}
\arguments{
   \item{x}{An mcmc.list object with more than one chain, and with starting values that are overdispersed with respect to the posterior distribution.}
   \item{confidence}{The coverage probability of the confidence interval for the potential scale reduction factor.}
}
\examples{
# This function should be used internally only. See the `MCMC` function for more details.
}
