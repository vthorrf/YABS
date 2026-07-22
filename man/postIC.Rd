\name{postIC}
\alias{postIC}
\title{Fully Bayesian Information Criterion}
\description{
This function computes multiple Bayesian model information criteria.
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
  \item{\code{WBIC}}{Widely Applicable Bayesian Information Criterion (Watanabe, 2013).}
  \item{\code{WAIC}}{Widely Applicable Information Criterion (Watanabe, 2010).}
  \item{\code{p_WAIC}}{Effective number of parameters from WAIC.}
  \item{\code{lppd}}{Log-pointwise predictive density.}
  \item{\code{TIC}}{Takeuchi Information Criterion (Takeuchi, 1976).}
  \item{\code{EIC}}{Extended Information Criterion (Ishiguro et al., 1997).}
  \item{\code{MDL}}{Fully Bayesian Minimum Description Length (Grünwald, 2007; Rissanen, 1989) on the deviance scale. Returns \code{NA} if \code{posterior} is not supplied.}
  \item{\code{DIC}}{The classic Deviance Information Criterion (Spiegelhalter et al., 2002).}
  \item{\code{p_DIC}}{Classic effective number of parameters for DIC.}
  \item{\code{DIC_p}}{Alternative DIC that uses the variance‑based penalty (Gelman et al., 2014).}
  \item{\code{p_V}}{Variance‑based effective number of parameters (Gelman et al., 2014).}
  \item{\code{DIC_i}}{Parameterization‑invariant DIC proposed by Xiao and Rabe-Hesketh (2026).}
}
}
\details{
  All criteria are returned on the \emph{deviance scale}, so smaller values
  indicate better model fit.

  \subsection{WBIC}{
    The Widely Applicable Bayesian Information Criterion (Watanabe, 2013) is
    defined as a weighted average of the log-likelihood at an inverse
    temperature \eqn{\beta = 1/\log(N)}:
    \deqn{\mathrm{WBIC} = -2 \sum_{s=1}^{S} w_s \log f(\mathbf{y}\mid\theta^{(s)}),}
    where the weights are \eqn{w_s \propto \exp\big((\beta-1)\log
    f(\mathbf{y}\mid\theta^{(s)})\big)}. WBIC is consistent for model
    selection even in singular learning problems.
  }

  \subsection{WAIC}{
    The Widely Applicable Information Criterion (Watanabe, 2010) is:
    \deqn{\mathrm{WAIC} = -2\big( \mathrm{lppd} - p_{\mathrm{WAIC}} \big)}
    with
    \eqn{\mathrm{lppd} = \sum_{i=1}^{N} \log \left( \frac{1}{S}
    \sum_{s=1}^{S} \exp(\ell_{si}) \right)} (the log pointwise predictive
    density) and
    \eqn{p_{\mathrm{WAIC}} = \sum_{i=1}^{N} \mathrm{Var}_s[\ell_{si}]}
    (the effective number of parameters). WAIC is asymptotically equivalent
    to leave-one-out cross-validation.
  }

  \subsection{TIC}{
    The Takeuchi Information Criterion (Takeuchi, 1976) is an extension of
    AIC that does not assume the true model is among the candidates. It
    uses the empirical Fisher information and the covariance of the score
    function. In the Bayesian context, we use the posterior means of the
    parameters \eqn{\hat{\theta} = \mathbb{E}_{\theta|\mathbf{y}}[\theta]}.
    The deviance is \eqn{D(\hat{\theta}) = -2\sum_{i=1}^{N}\ell_i(\hat{\theta})},
    where \ell_i(\hat{\theta}) = \log f(y_i \mid \hat{\theta}). The penalty
    is twice the sum over observations of the posterior variance of the
    pointwise log-likelihood:
    \deqn{p_{\mathrm{TIC}} = \sum_{i=1}^{N} \mathrm{Var}_{\theta|\mathbf{y}}
    [\ell_i(\theta)].}
    Thus,
    \deqn{\mathrm{TIC} = D(\hat{\theta}) + 2 p_{\mathrm{TIC}}.}
    In our implementation, we use the posterior means of the pointwise
    log-likelihoods rather than plugging \eqn{\hat{\theta}} into the model,
    which is computationally convenient and asymptotically equivalent.
  }

  \subsection{EIC}{
    The Extended Information Criterion (Ishiguro et al., 1997)
    is a bootstrap-based variant of TIC. It estimates the optimism by
    resampling from the posterior predictive distribution. Our implementation
    uses a simple parametric bootstrap where for each MCMC draw \eqn{\theta^{(s)}},
    we compute the average log-likelihood of the observed data
    \eqn{\bar{\ell}_s = \frac{1}{N}\sum_{i=1}^{N}\ell_i(\theta^{(s)})} and
    the overall mean \eqn{\bar{\ell} = \frac{1}{S}\sum_{s=1}^{S}\bar{\ell}_s}.
    The bias estimate is the average deviation:
    \deqn{\widehat{\mathrm{bias}} = \frac{1}{S}\sum_{s=1}^{S}
    \big( \bar{\ell}_s - \bar{\ell} \big).}
    The deviance at the posterior mean pointwise log-likelihoods is
    \eqn{D(\hat{\ell}) = -2\sum_{i=1}^{N}\bar{\ell}_i}, where
    \eqn{\bar{\ell}_i = \frac{1}{S}\sum_{s=1}^{S}\ell_i(\theta^{(s)})}.
    Then,
    \deqn{\mathrm{EIC} = D(\hat{\ell}) - 2N \cdot \widehat{\mathrm{bias}}.}
    EIC tends to be more robust than TIC when the model is misspecified.
  }

  \subsection{Bayesian MDL}{
    The Bayesian Minimum Description Length (Rissanen, 1989; Grünwald, 2007)
    is approximated using a Laplace expansion:
    \deqn{\mathrm{MDL} \approx -2\Bigg[ \mathbb{E}_{\theta|\mathbf{y}}[\log
    f(\mathbf{y}\mid\theta)] - \frac{d}{2}\log N + \frac{1}{2}\log|\Sigma|
    + \frac{d}{2}\log(2\pi e) \Bigg],}
    where \eqn{d} is the number of parameters and \eqn{\Sigma} is the posterior
    covariance matrix. The term \eqn{\mathbb{E}_{\theta|\mathbf{y}}[\log
    f(\mathbf{y}\mid\theta)]} is approximated by the average of the total
    log-likelihood across MCMC draws. The criterion is returned on the
    deviance scale (multiplied by \eqn{-2}) to be comparable with other
    criteria.
  }

  \subsection{DIC family}{
    The classic DIC (Spiegelhalter et al., 2002) is:
    \deqn{\mathrm{DIC} = D(\bar{\theta}) + 2p_{\mathrm{DIC}},}
    where \eqn{p_{\mathrm{DIC}} = \mathbb{E}[D(\theta)] - D(\bar{\theta})}.
    It relies on the plug-in deviance \eqn{D(\bar{\theta})} and can be
    unstable or negative when the posterior is multimodal (e.g., label
    switching, sign switching, degenerate nonidentifiability).

    Gelman, Hwang, & Vehtari (2014) proposed a variance-based penalty:
    \deqn{p_V = \frac{1}{2}\mathrm{Var}_{\theta|\mathbf{y}}[D(\theta)],}
    which is always non‑negative and invariant to reparameterization.
    They defined \eqn{\mathrm{DIC}_p = D(\bar{\theta}) + 2p_V}.

    Xiao & Rabe-Hesketh (2026) introduced \eqn{\mathrm{DIC}_i}, which
    removes the plug-in deviance entirely:
    \deqn{\mathrm{DIC}_i = \mathbb{E}_{\theta|\mathbf{y}}[D(\theta)] + p_V.}
    This criterion is parameterization‑invariant, never produces a negative
    penalty, and is asymptotically equivalent to WAIC. It is recommended
    when multimodality is present or when the likelihood is not easily
    factorizable.
  }
}
\references{
Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive information criteria for Bayesian models. \emph{Statistics and Computing, 24}, 997–1016.

Grünwald, P. D. (2007). \emph{The Minimum Description Length Principle}. MIT Press.

Ishiguro, M., Sakamoto, Y., & Kitagawa, G. (1997). Bootstrapping loglikelihood and EIC, an extension of AIC. \emph{Annals of the Institute of Statistical Mathematics, 49}, 411–434.

Rissanen, J. (1989). \emph{Stochastic Complexity in Statistical Inquiry}. World Scientific.

Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & van der Linde, A. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society, Series B, 64}, 583–639.

Takeuchi, K. (1976). Distribution of information statistics and a criterion for model selection. \emph{Surikaisekikenkyusho Kokyuroku (in Japanese), 243}, 1–18.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. \emph{Journal of Machine Learning Research, 11}, 3571–3594.

Watanabe, S. (2013). A widely applicable Bayesian information criterion. \emph{Journal of Machine Learning Research, 14}(1), 867–897.

Xiao, X., & Rabe-Hesketh, S. (2026). A Parameterization-Invariant DIC. \emph{arXiv preprint}, \url{https://arxiv.org/html/2605.27844v1}.
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