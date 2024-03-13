\name{MCMC}
\alias{MCMC}
\title{Call the desired sampler}
\description{
The `MCMC` function is a basic user interface for calling the samplers provided with the YABS package. The user provides a Model, a Data list, and Initial values (optional if PGF is provided). The function compiles the information and sends it to internal samplers, written in C++. The results are consolidated and summarized to be output as an object of class YABS.
}
\usage{
MCMC(Model, Data, Initial.Values=NULL, iterations=NULL, burnin=NULL, status=NULL,
     thinning=NULL, adapt=NULL, nchains=1, parallel=FALSE, cores=NULL, update.progress=NULL,
     pckgs=NULL, opt.init=TRUE, par.cov=NULL, hessian=FALSE, control=list(fnscale=-1),
     algo=c("mwg","rwm","barker","ohss"))
}
\arguments{
   \item{Model}{The user provided model. It should be a function that always output a list containing the following named objects: (i) LP, the logposterior; (ii) Dev, the deviance; (iii) Monitor, a single value or a vector of values to be monitored; (iv) yhat, values simulated from the posterior; and (v) parm, the value or vector with the parameters.}
   \item{Data}{A named list of the data objects required by the model. Should always include: (i) parm.names, the name of the parameters; (ii) pos.parameter, the positions of each group of parameters (see the Example); (iii) PGF, the Parameter-Generating Function (PGF) if initial values for the parameters are not provided; and (iv) mon.names, the names of values monitored (see the Example).}
   \item{Initial.Values}{A vector of initial values equal in length to the number of parameters. If it is not provided, it will use the PGF provided in Data to generated valid initial values.}
   \item{iterations}{The number of iterations to be saved (excluding burning and thinning).}
   \item{burnin}{The number of initial iterations to be discarded (excluding thinning).}
   \item{status}{Indicates how often, in iterations, the user would like the status printed to the screen. Should be an integer between 1 and the number of iterations (including burning and thinning).}
   \item{thinning}{Indicates that every nth iteration will be retained, while the other iterations are discarded. Must be a positive integer.}
   \item{adapt}{The number of iterations to run in the adaptive phase.}
   \item{nchains}{The number of Markov chains to run.}
   \item{parallel}{If TRUE, run MCMC chains in parallel on multiple CPU cores.}
   \item{cores}{If parallel=TRUE, specify the number of CPU cores used.}
   \item{update.progress}{How many times should the progress bar be updated when running in parallel? Note that setting this to a large value should result in the worse performance, due to additional overhead communicating among the parallel processes.}
   \item{pckgs}{If parallel=TRUE and if `Model` uses non-base-R functions, this is a character vector with the name of the packages from where the functions come from.}
   \item{opt.init}{If TRUE, the initial values will be improved (i.e., set closer to the maximum a posteriori) and the initial step sizes will be estimated with the negative inverse Hessian of the MAP estimates.}
   \item{par.cov}{An optional covariance matrix of the parameters (i.e., the diagonal represents the step sizes for the parameters to be estimated using the MCMC algorithms).}
   \item{hessian}{If opt.init = TRUE, this argument indicates if a numerically differentiated Hessian matrix be used for the initial step sizes. Otherwise, only the diagonal matrix of the second-order partial derivatives will be used for the initial step sizes.}
   \item{control}{If opt.init = TRUE, a list of control parameters to be passed to the `optim` function from the stats package. Should always include fnscale = -1.}
   \item{algo}{The MCMC algorithm to be used. Currently, the following are implemented: (i) "mwg", Metropolis-Within-Gibbs; (ii) "rwm", Random-walk Metropolis; (iii) "barker", Barker Proposal Metropolis; and (iv) "ohss", Oblique Hyperrectangle Slice Sampler. If no (or more than one) algorithm is selected, defaults to "rwm".}
}
\value{
A list of class YABS containing:
1. posterior, the posterior distribution of the parameters, as well as the deviance and AIC;
2. LP, the posterior distribution of the monitored values;
3. DIC, the calculated DIC and its components (Dbar and pD);
4. acc, the achieved acceptance rate;
5. ESS, Effective Sample Size estimates for each parameter;
6. PSRF, the 'potential scale reduction factor' for each parameter;
7. MPSRF, the 'multivariate potential scale reduction factor', if the model contains more than one parameter; and
8. mcmc.info, a list containing: (i) algorithm, the name of the MCMC algorithm used; (ii) n.iter, the total number of iterations after including burning and excluding thinning; (iii) n.burnin, the number of burning samples; (iv) n.thin, the magnitude of thinning; and (v) elapsed.mins, the time the model took to run in minutes.
}
\references{
Livingstone, S., & Zanella, G. (2019). The Barker proposal: combining robustness and efficiency in gradient-based MCMC. JRSS: Series B (Statistical Methodology), 84(2), 496-523.
Thompson, M. D. (2011). Slice Sampling with Multivariate Steps. http://hdl.handle.net/1807/31955
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
fit1 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="mwg")
fit2 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="rwm")
fit3 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="barker")
fit4 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="ohss")
plot.ts(fit1$posterior)
plot.ts(fit2$posterior)
plot.ts(fit3$posterior)
plot.ts(fit4$posterior)
fit1
fit2
fit3
fit4
}
