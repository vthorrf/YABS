\name{VB}
\alias{VB}
\title{Call the desired variational inference algorithm}
\description{
The `VB` function is a basic user interface for calling the variational inference algorithms provided with the YABS package. The user provides a Model, a Data list, and Initial values (optional if PGF is provided). The function optionally improves the initial values through multi-start maximum a posteriori estimation, estimates or regularizes an initial covariance matrix, and sends the information to internal algorithms written in C++. The results are consolidated and summarized to be output as an object of class YABS_VB.
}
\usage{
VB(Model, Data, Initial.Values = NULL, iterations = NULL,
   status = NULL, thinning = NULL, adapt = NULL,
   method = c("sagva", "qnsagva", "mcvi", "svgd"),
   n.particles = 50, par.cov = NULL,
   w.min = 1e-8, w.max = 1.0,
   armijo.c1 = 1e-4, wolfe.c2 = 0.9,
   ls.max.iter = 20,
   h = 1e-4, stop.tolerance = 1e-4,
   learning.rate = NULL, min.iter = NULL, patience = 50,
   memory = 10, damping = 0.2,
   T.mcmc = 3, L.leapfrog = 3, k.rb = 5,
   mcvi.spsa.scale = 1e-3, mcvi.epsilon = NULL,
   step.size = 0.01, h.scale = 1.0, use.adam = TRUE,
   opt.init = TRUE, n.starts = 3, init.jitter = 0.5,
   particle.scale = 1.5, retry.on.failure = TRUE,
   control = list(), seed = NULL)
}
\arguments{
   \item{Model}{The user-provided model. It should be a function that always outputs a list containing the following named objects: (i) LP, the log posterior; (ii) Dev, the deviance; (iii) Monitor, a single value or a vector of values to be monitored; (iv) yhat, values simulated from the posterior predictive distribution; and (v) parm, the value or vector with the parameters on the scale to be returned to the user.}
   \item{Data}{A named list of the data objects required by the model. It should include: (i) parm.names, the names of the parameters returned in Model$parm; (ii) PGF, the Parameter-Generating Function if initial values are not provided; and, optionally, (iii) latent.parm.names, the names of the parameters on the unconstrained scale supplied to Model. Other model-specific data objects can also be included.}
   \item{Initial.Values}{A finite numeric vector of initial values equal in length to the number of latent parameters. If it is not provided, the function uses the PGF supplied in Data to generate valid initial values.}
   \item{iterations}{The maximum number of optimization iterations. For SAGVA, QN-SAGVA, and MCVI, the number of returned posterior draws is approximately iterations divided by thinning. For SVGD, the returned posterior approximation consists of the final particles.}
   \item{status}{Indicates how often, in optimization iterations, the user would like the status printed to the screen. It must be a positive integer.}
   \item{thinning}{For SAGVA, QN-SAGVA, and MCVI, determines the number of posterior draws generated from the final approximation, with approximately one draw for every nth optimization iteration. It also controls internal storage for SVGD. It must be a positive integer.}
   \item{adapt}{A backward-compatible control used to define the default minimum number of optimization iterations. If min.iter is not provided, min.iter defaults to the larger of 50 and adapt.}
   \item{method}{The variational inference algorithm to be used. Currently, the following are implemented: (i) "sagva", Stochastic Approximation for Gaussian Variational Approximation; (ii) "qnsagva", an experimental limited-memory quasi-Newton version of SAGVA; (iii) "mcvi", Hamiltonian Markov Chain Variational Inference; and (iv) "svgd", Stein Variational Gradient Descent. If more than one method is supplied, the first valid option selected by match.arg is used; the default is "sagva".}
   \item{n.particles}{The number of particles used by SVGD. The final particles constitute the returned empirical posterior approximation. The minimum is two.}
   \item{par.cov}{An optional covariance matrix matching the length of Initial.Values. It is regularized to be positive definite and is used to initialize the Gaussian approximation, the MCVI initial distribution, or the SVGD particle population. If omitted, a Laplace covariance is estimated from the Hessian at the best available MAP solution.}
   \item{w.min}{Retained for backward compatibility with the preliminary SAGVA interface. The corrected SAGVA algorithm does not use a line search, and this argument is currently ignored.}
   \item{w.max}{Retained for backward compatibility with the preliminary SAGVA interface. The corrected SAGVA algorithm does not use a line search, and this argument is currently ignored.}
   \item{armijo.c1}{Retained for backward compatibility with the preliminary SAGVA interface. The corrected SAGVA algorithm does not use an Armijo line search, and this argument is currently ignored.}
   \item{wolfe.c2}{Retained for backward compatibility with the preliminary SAGVA interface. The corrected SAGVA algorithm does not use Wolfe conditions, and this argument is currently ignored.}
   \item{ls.max.iter}{Retained for backward compatibility with the preliminary SAGVA interface. The corrected SAGVA algorithm does not use a line search, and this argument is currently ignored.}
   \item{h}{The base finite-difference step used to calculate numerical gradients and, when required, Hessians. The internal C++ routines scale the step according to the magnitude of each parameter.}
   \item{stop.tolerance}{The convergence tolerance. Early stopping is considered only after min.iter and requires the convergence criterion to remain below this value for patience consecutive checks. Set stop.tolerance = 0 to disable tolerance-based early stopping and use the complete optimization budget.}
   \item{learning.rate}{An optional positive learning rate. If NULL, a method-specific default or decreasing schedule is used. For MCVI, the wrapper uses a default value of 0.005 when no value is supplied.}
   \item{min.iter}{The minimum number of optimization iterations before early stopping can occur. If NULL, it defaults to the larger of 50 and adapt, bounded by iterations.}
   \item{patience}{The number of consecutive convergence checks that must satisfy stop.tolerance before estimation is stopped. This prevents termination after a single unusually small update.}
   \item{memory}{For QN-SAGVA, the maximum number of stochastic secant pairs retained by the limited-memory quasi-Newton approximation.}
   \item{damping}{For QN-SAGVA, the Powell damping constant used to stabilize curvature updates and preserve a positive-definite approximation.}
   \item{T.mcmc}{For MCVI, the number of sequential Hamiltonian transformations applied to each draw from the initial variational distribution.}
   \item{L.leapfrog}{For MCVI, the number of leapfrog steps used within each Hamiltonian transformation.}
   \item{k.rb}{For MCVI, the number of Monte Carlo trajectories used at each optimization iteration to estimate the auxiliary variational lower bound and its updates.}
   \item{mcvi.spsa.scale}{For MCVI, the perturbation magnitude used by the common-random-number Simultaneous Perturbation Stochastic Approximation update. It must be positive.}
   \item{mcvi.epsilon}{An optional positive initial leapfrog step size for MCVI. If NULL, a reasonable initial value is estimated automatically. When numerical recovery is required, the step size is reduced.}
   \item{step.size}{The SVGD update step size. When use.adam = TRUE, this scales the Adam update.}
   \item{h.scale}{A positive multiplier applied to the median-heuristic RBF kernel bandwidth used by SVGD.}
   \item{use.adam}{If TRUE, SVGD uses Adam to scale the particle updates. If FALSE, a direct fixed-step update is used.}
   \item{opt.init}{If TRUE, the initial values are improved using multi-start MAP estimation. BFGS is attempted first, with a Nelder-Mead fallback when necessary. The best finite solution is used to initialize variational inference.}
   \item{n.starts}{The number of starting values used during multi-start MAP estimation. The first is Initial.Values and the remaining starts are randomly perturbed versions of it.}
   \item{init.jitter}{The scale of the random perturbations used to generate additional MAP starting values.}
   \item{particle.scale}{For SVGD, the scale applied to the initial covariance when drawing particles around the MAP estimate. Most particles are drawn from this scale and a smaller group is initialized more broadly to reduce sensitivity to local modes.}
   \item{retry.on.failure}{If TRUE, a numerically unstable first variational run is retried using a contracted initial covariance. For MCVI and SVGD, the leapfrog or particle step size is also reduced when appropriate.}
   \item{control}{If opt.init = TRUE, a list of control parameters passed to the optim function from the stats package. The wrapper minimizes the negative log posterior internally, so fnscale should not be changed to -1.}
   \item{seed}{An optional integer seed passed to set.seed before initialization and estimation.}
}
\value{
A list of class YABS_VB containing:
1. posterior, posterior draws of the parameters returned by Model$parm; for SVGD, these are the final particles;
2. yhat, posterior predictive values associated with the returned approximation;
3. LP, the log posterior evaluated for each returned draw or particle;
4. Monitor, the monitored model quantities associated with the returned approximation;
5. DIC, the variance-based DIC and its components Dbar and pD;
6. ESS, the number of generated variational draws for Gaussian and MCVI approximations, and NA for SVGD because the final particles are not a Markov chain;
7. PSRF and MPSRF, returned as NA because MCMC potential scale reduction factors are not defined for a single variational approximation;
8. vb.info, a list containing: (i) method, the variational algorithm used; (ii) optimization.iterations, the number of optimization iterations completed; (iii) posterior.draws, the number of returned draws or particles; (iv) thinning; (v) min.iter; (vi) patience; (vii) n.particles for SVGD; (viii) ELBO, when available; and (ix) elapsed.mins;
9. latent.posterior, the posterior approximation on the latent or unconstrained parameter scale supplied to Model;
10. variational.mean and variational.VarCov, the estimated mean and covariance when available;
11. MAP, information about the selected MAP initialization when opt.init = TRUE;
12. diagnostics, a list containing the convergence flag, final convergence criterion, completed optimization iterations, and number of invalid numerical updates or trajectories;
13. accepted.curvature.pairs, for QN-SAGVA, the number of accepted quasi-Newton curvature pairs;
14. q0.mean, q0.VarCov, inverse.model, hmc.step.size, and ELBO.history for MCVI; and
15. final.particles, latent.particles, particles.lp, kernel.bandwidth, and Stein.update.history for SVGD.
}
\details{
SAGVA estimates a full Gaussian variational approximation using stochastic approximation based on posterior gradients and negative Hessian information. Because a full numerical Hessian is calculated, its computational cost increases rapidly with the number of parameters.

QN-SAGVA is an experimental extension implemented in YABS. It replaces the repeated full Hessian calculation with a Powell-damped limited-memory quasi-Newton approximation. It can be more efficient in larger models, but should be validated carefully because it is not one of the three literature-based algorithms.

MCVI uses a diagonal-Gaussian initial approximation followed by a sequence of volume-preserving Hamiltonian transformations. A learned reverse model is included in the auxiliary variational lower bound. Because Model is supplied as a black-box R function rather than an automatic-differentiation graph, some MCVI parameters are optimized using common-random-number SPSA.

SVGD represents the posterior with an empirical population of interacting particles. It therefore does not provide a normalized variational density or a directly comparable ELBO. The returned posterior for this method is the final particle population, rather than particles collected at different optimization iterations.
}
\references{
Salimans, T., & Knowles, D. A. (2013). Fixed-form variational posterior approximation through stochastic linear regression. Bayesian Analysis, 8(4), 837-882.

Salimans, T., Kingma, D. P., & Welling, M. (2015). Markov Chain Monte Carlo and Variational Inference: Bridging the Gap. Proceedings of the 32nd International Conference on Machine Learning, 37, 1218-1226.

Liu, Q., & Wang, D. (2016). Stein Variational Gradient Descent: A General Purpose Bayesian Inference Algorithm. Advances in Neural Information Processing Systems, 29, 2378-2386.

Nocedal, J. (1980). Updating quasi-Newton matrices with limited storage. Mathematics of Computation, 35(151), 773-782.
}
\examples{##### Multiple Linear Regression with 5 Predictors
### Packages and functions====
require(YABS)
require(compiler)

### Random data====
seed <- 1234; N <- 200; V <- 5
set.seed(seed)
X <- sapply(1:V, function(g) rnorm(N))
betas <- runif(V + 1, .3, .7)
y <- rnorm(N, cbind(1, X) \%*\% betas, 1)

### DATA LIST====
parm.names <- c(paste("beta", 0:V, sep = ""), "sigma")
latent.parm.names <- c(paste("beta", 0:V, sep = ""), "log_sigma")
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
  beta <- rnorm(ncol(Data$X) + 1)
  log_sigma <- rnorm(1)
  return(c(beta, log_sigma))
}
Data <- list(X = X, y = y, n = N, parm.names = parm.names,
             latent.parm.names = latent.parm.names,
             pos.beta = pos.beta, pos.sigma = pos.sigma,
             PGF = PGF, mon.names = "LP")
Initial.Values <- PGF(Data)

### MODEL====
Model <- function(parm, Data) {
  ### Parameters
  beta <- parm[Data$pos.beta]
  sigma <- exp(parm[Data$pos.sigma])

  ### Log-Prior
  beta.prior <- sum(dnorm(beta, 0, 1, log = TRUE))
  sigma.prior <- sum(dgamma(sigma, 1e-2, 1e-2, log = TRUE)) +
    parm[Data$pos.sigma]
  Lp <- beta.prior + sigma.prior

  ### Log-Likelihood
  mu <- cbind(1, Data$X) \%*\% beta
  pw_ll <- dnorm(Data$y, mu, sigma, log = TRUE)
  LL <- sum(pw_ll)

  ### Log-Posterior
  LP <- LL + Lp
  Modelout <- list(LP = LP, Dev = -2 * LL, Monitor = pw_ll,
                   yhat = rnorm(Data$n, mu, sigma),
                   parm = c(beta, sigma))
  return(Modelout)
}
Model <- compiler::cmpfun(Model)

### FIT====
\donttest{
# Stochastic Approximation for Gaussian Variational Approximation
fit1 <- VB(Model, Data, Initial.Values = NULL, iterations = 500,
           status = 100, thinning = 2, method = "sagva", seed = 1234)

# Experimental limited-memory quasi-Newton SAGVA
fit2 <- VB(Model, Data, Initial.Values = NULL, iterations = 500,
           status = 100, thinning = 2, method = "qnsagva",
           memory = 10, seed = 1234)

# Markov Chain Variational Inference
fit3 <- VB(Model, Data, Initial.Values = NULL, iterations = 300,
           status = 60, thinning = 2, method = "mcvi",
           T.mcmc = 2, L.leapfrog = 3, k.rb = 3, seed = 1234)

# Stein Variational Gradient Descent
fit4 <- VB(Model, Data, Initial.Values = NULL, iterations = 500,
           status = 100, method = "svgd", n.particles = 50,
           step.size = 0.01, seed = 1234)

summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
}
}
