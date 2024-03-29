\name{LA}
\alias{LA}
\title{Laplace Approximation}
\description{
The `LA` function is a basic user interface for calling Laplace Approximation method. The user provides a Model, a Data list, and Initial values (optional if PGF is provided). The function compiles the information and sends it to the `optim` function of the base package `stats`. The results are consolidated and summarized to be output as an object of class YABS_LA.
}
\usage{
LA(Model, Data, Initial.Values=NULL, lower = -Inf, upper = Inf,
   control = list(), hessian = FALSE, par.cov=NULL, SIR=TRUE, 
   method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
   iterations=NULL, nearPD=TRUE, check.convergence=FALSE)
}
\arguments{
   \item{Model}{The user provided model. It should be a function that always output a list containing the following named objects: (i) LP, the logposterior; (ii) Dev, the deviance; (iii) Monitor, a single value or a vector of values to be monitored; (iv) yhat, values simulated from the posterior; and (v) parm, the value or vector with the parameters.}
   \item{Data}{A named list of the data objects required by the model. Should always include: (i) parm.names, the name of the parameters; (ii) pos.parameter, the positions of each group of parameters (see the Example); (iii) PGF, the Parameter-Generating Function (PGF) if initial values for the parameters are not provided; and (iv) mon.names, the names of values monitored (see the Example).}
   \item{Initial.Values}{A vector of initial values equal in length to the number of parameters. If it is not provided, it will use the PGF provided in Data to generated valid initial values.}
   \item{lower}{Lower bound on the parameters for the "L-BFGS-B" method, or bounds in which to search for method "Brent".}
  \item{upper}{Upper bound on the parameters for the "L-BFGS-B" method, or bounds in which to search for method "Brent".}
  \item{control}{a list of control parameters. See ‘Details’ of the `optim` function.}
  \item{hessian}{Logical. Should a numerically differentiated Hessian matrix be use to estimate the variance of the parameters?}
  \item{par.cov}{An optional covariance matrix of the parameters (i.e., the diagonal represents the variance of the parameters to be estimated).}
  \item{SIR}{This logical argument indicates whether or not Sampling-Importance Resampling (SIR) is conducted to draw independent posterior samples. This argument defaults to TRUE.}
  \item{method}{The method to be used. See ‘Details’ of the `optim` function. Can be abbreviated.}
  \item{iterations}{The number of samples to be drawn if SIR is TRUE.}
  \item{nearPD}{This logical argument indicates whether or not the covariance matrix of the parameters should be transformed for the nearest positive-definite matrix. If this argument is set to FALSE and the estimated Hessian is not positive-definite, the function will throw an error.}
  \item{check.convergence}{This logical argument indicates whether the function should check if convergence was achieved before doing SIR. This argument defaults to FALSE. If it is set to TRUE and the optimization routine does not achieve convergence, the user will have to indicate if they desire to continue with SIR anyway, if that argument is set to TRUE.}
}
\value{
A list of class YABS_LA containing:
1. MAP, the maximum a posteriori estimates of the parameters;
2. yhat, the yhat as calculated by the Model;
3. LP, the logposterior of the MAP estimates;
4. Monitor, the monitor values as calculated by the Model;
5. Dev, the deviance of the MAP estimates;
6. AIC, the Akaike Information criterion for the MAP estimates;
7. convergence, a logical indicator if the optimization routine converged.

If SIR is set to TRUE, it also contains:
1. posterior, the posterior distribution of the parameters;
2. DIC, the calculated DIC and its components (Dbar and pD);
3. ESS, Effective Sample Size estimates for each parameter.
And yhat, LP, Monitor, Dev, and AIC are matrices/vectors of the posterior samples.
}
\references{
Azevedo-Filho, A., & Shachter, R. D. (1994, January). Laplace's method approximations for probabilistic inference in belief networks with continuous variables. In Uncertainty in Artificial Intelligence (pp. 28-36). Morgan Kaufmann.
Laplace, P. S. (1986). Mémoires de Mathématique et de Physique, Tome Sixieme [Memoir on the probability of causes of events.]. Stat. Sci., 366-367.
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
fit <- LA(Model, Data, Initial.Values=NULL, iterations=1000,
          hessian=TRUE, method="BFGS")
plot.ts(fit$posterior)
fit
}
