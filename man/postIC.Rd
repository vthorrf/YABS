\name{postIC}
\alias{postIC}
\title{Posterior Information Criteria, ICOMP Analogues, and Predictive Measures}
\description{
Computes Bayesian predictive criteria, posterior deviance summaries,
importance-sampling diagnostics, and posterior correlation-based analogues of
ICOMP criteria from pointwise posterior log-likelihood values and, optionally,
posterior parameter draws.
}
\usage{
postIC(x, posterior = NULL)
}
\arguments{
  \item{x}{A required numeric \eqn{S \times N} matrix of pointwise
  log-likelihood values. Rows correspond to posterior draws and columns
  correspond to observations or independent predictive units.}

  \item{posterior}{An optional numeric \eqn{S \times d} matrix of posterior
  parameter draws. It must have the same number of rows as \code{x}. When
  supplied, it is used to calculate the posterior parameter correlation
  complexity and the posterior ICOMP analogues. Fixed or deterministic
  parameter columns should not be included.}
}
\value{
A list containing:
\itemize{
  \item{\code{Dbar}}{Posterior mean deviance.}

  \item{\code{WBIC}}{WBIC estimated by self-normalized importance sampling
  from the ordinary posterior to the tempered posterior.}
  \item{\code{WBIC_beta}}{Inverse temperature \eqn{1/\log(N)}.}
  \item{\code{WBIC_ESS}}{Effective sample size of the normalized WBIC
  importance weights.}
  \item{\code{WBIC_relative_ESS}}{WBIC effective sample size divided by the
  number of posterior draws.}
  \item{\code{WBIC_max_weight}}{Largest normalized WBIC importance weight.}
  \item{\code{WBIC_perplexity}}{Perplexity of the normalized WBIC weights.}
  \item{\code{WBIC_relative_perplexity}}{WBIC perplexity divided by the
  number of posterior draws.}

  \item{\code{lppd}}{Log pointwise predictive density.}
  \item{\code{WAIC}}{Backward-compatible alias for \code{WAIC2}.}
  \item{\code{WAIC1}}{WAIC using the first common effective-parameter
  formulation.}
  \item{\code{p_WAIC1}}{Effective number of parameters used by
  \code{WAIC1}.}
  \item{\code{WAIC2}}{WAIC using the posterior variance of pointwise
  log-likelihoods.}
  \item{\code{p_WAIC2}}{Effective number of parameters used by
  \code{WAIC2}.}

  \item{\code{CPO}}{Raw importance-sampling conditional predictive
  ordinates.}
  \item{\code{log_CPO}}{Log conditional predictive ordinates. This is
  numerically preferable to \code{CPO}.}
  \item{\code{LPML}}{Log pseudo-marginal likelihood, equal to the sum of
  \code{log_CPO}. Larger values indicate better predictive performance.}
  \item{\code{LOOIC_IS}}{Raw importance-sampling leave-one-out information
  criterion, equal to \eqn{-2\,\mathrm{LPML}}.}
  \item{\code{p_LOO_IS}}{Raw importance-sampling estimate of the effective
  number of parameters, \eqn{\mathrm{lppd}-\mathrm{LPML}}.}

  \item{\code{LOOIC}}{PSIS-LOO information criterion when the optional
  \pkg{loo} package is installed; otherwise \code{NA}.}
  \item{\code{SE_LOOIC}}{Estimated standard error of \code{LOOIC}.}
  \item{\code{elpd_LOO}}{Expected log predictive density from PSIS-LOO.}
  \item{\code{SE_elpd_LOO}}{Estimated standard error of
  \code{elpd_LOO}.}
  \item{\code{p_LOO}}{Effective number of parameters from PSIS-LOO.}
  \item{\code{SE_p_LOO}}{Estimated standard error of \code{p_LOO}.}
  \item{\code{Pareto_k}}{Observation-level Pareto shape diagnostics when
  \pkg{loo} is installed; otherwise \code{NULL}.}
  \item{\code{PSIS_LOO}}{The complete object returned by \code{loo::loo()},
  or \code{NULL} when \pkg{loo} is unavailable.}

  \item{\code{DIC_i}}{Parameterization-invariant, plug-in-free DIC.}
  \item{\code{p_V}}{Variance-based effective number of parameters,
  \eqn{\tfrac{1}{2}\mathrm{Var}\{D(\theta)\}}.}

  \item{\code{ICOMP}}{Posterior correlational analogue of
  ICOMP based on inverse Fisher information.}
  \item{\code{ICOMP_PEU}}{Posterior correlational analogue of the
  posterior-expected-utility form of ICOMP. This is the principal ICOMP
  result returned by the function.}
  \item{\code{ICOMP_PEU_LN}}{Posterior correlational analogue of the
  log-sample-size-weighted ICOMP-PEU criterion.}
  \item{\code{CICOMP}}{Posterior correlational analogue of consistent
  ICOMP.}
  \item{\code{C_R}}{Posterior correlational complexity,
  \eqn{-\tfrac{1}{2}\log|R_{\mathrm{post}}|}, where
  \eqn{R_{\mathrm{post}}} is the posterior parameter correlation matrix.}
}
}
\details{
The input must contain the same pointwise likelihood factorization for every
model being compared. For hierarchical or latent-variable models, a column of
\code{x} should represent the predictive unit that would be left out during
cross-validation, such as an entire participant or cluster when appropriate.

Criteria expressed on the deviance scale, including \code{WBIC}, \code{WAIC},
\code{LOOIC_IS}, \code{LOOIC}, \code{DIC_i}, \code{ICOMP},
\code{ICOMP_PEU}, \code{ICOMP_PEU_LN}, and \code{CICOMP}, are minimized.
Log predictive quantities such as \code{lppd}, \code{LPML}, and
\code{elpd_LOO} are maximized.

\subsection{WBIC}{
The widely applicable Bayesian information criterion is defined as the
posterior expectation of the negative twice log-likelihood at inverse
temperature
\deqn{\beta=1/\log(N).}

This function receives draws from the ordinary posterior, not draws generated
directly from the tempered posterior. It therefore estimates WBIC by
self-normalized importance sampling with weights
\deqn{w_s \propto \exp\{(\beta-1)\ell(\theta^{(s)})\}.}

The approximation may be unreliable when the tempered posterior is poorly
represented by ordinary posterior draws. Users should inspect
\code{WBIC_relative_ESS}, \code{WBIC_max_weight}, and
\code{WBIC_relative_perplexity}. Concentrated weights indicate that the WBIC
estimate is based effectively on only a small fraction of the posterior draws.
For \eqn{N=1}, WBIC and its diagnostics are returned as \code{NA}.
}

\subsection{WAIC}{
Let \eqn{\ell_{si}} be the log-likelihood contribution of predictive unit
\eqn{i} at posterior draw \eqn{s}. The log pointwise predictive density is
\deqn{\mathrm{lppd}=\sum_{i=1}^{N}\log\left\{\frac{1}{S}
\sum_{s=1}^{S}\exp(\ell_{si})\right\}.}

The function returns both common WAIC penalty estimates:
\deqn{p_{\mathrm{WAIC},1}=2\sum_{i=1}^{N}\left[
\log E_{\theta\mid y}\{p(y_i\mid\theta)\}-
E_{\theta\mid y}\{\log p(y_i\mid\theta)\}\right],}
\deqn{p_{\mathrm{WAIC},2}=\sum_{i=1}^{N}
\mathrm{Var}_{\theta\mid y}\{\ell_i(\theta)\}.}

The corresponding criteria are
\deqn{\mathrm{WAIC}_j=-2(\mathrm{lppd}-p_{\mathrm{WAIC},j}).}
The element named \code{WAIC} remains an alias for the variance-based
\code{WAIC2} to preserve compatibility with earlier versions of YABS.
}

\subsection{CPO, LPML, and raw importance-sampling LOO}{
The conditional predictive ordinate is the leave-one-out predictive density
\deqn{\mathrm{CPO}_i=p(y_i\mid y_{-i}).}
Using draws from the full posterior, it can be estimated through the
harmonic-mean identity
\deqn{\widehat{\mathrm{CPO}}_i=
\left\{\frac{1}{S}\sum_{s=1}^{S}
\frac{1}{p(y_i\mid\theta^{(s)})}\right\}^{-1}.}

The log pseudo-marginal likelihood is
\deqn{\mathrm{LPML}=\sum_{i=1}^{N}\log(\mathrm{CPO}_i),}
and \code{LOOIC_IS} is \eqn{-2\,\mathrm{LPML}}. These quantities target
leave-one-out predictive performance, but the raw harmonic-mean importance
ratios can have high or infinite variance. They should therefore be interpreted
cautiously, especially when they differ substantially from PSIS-LOO or WAIC.
}

\subsection{PSIS-LOO}{
When the optional \pkg{loo} package is installed, \code{postIC()} also calls
\code{loo::loo()} on the pointwise log-likelihood matrix. Pareto-smoothed
importance sampling stabilizes the largest importance ratios and supplies
observation-level Pareto \eqn{\hat{k}} diagnostics.

Because the current \code{postIC()} interface does not receive chain
identifiers, the call uses \code{r_eff = 1}. This assumption mainly affects
Monte Carlo effective-sample-size and standard-error calculations when MCMC
draws are autocorrelated. Users requiring chain-aware diagnostics should call
\code{loo::loo()} directly with an appropriate relative-efficiency estimate.
The complete \pkg{loo} object is retained in \code{PSIS_LOO}.
}

\subsection{Parameterization-invariant DIC}{
For posterior deviance draws \eqn{D(\theta^{(s)})}, define
\deqn{\bar D=E_{\theta\mid y}\{D(\theta)\},\qquad
p_V=\frac{1}{2}\mathrm{Var}_{\theta\mid y}\{D(\theta)\}.}
The plug-in-free criterion is
\deqn{\mathrm{DIC}_i=\bar D+p_V.}
It is invariant to one-to-one reparameterizations and avoids evaluating the
deviance at the posterior mean parameters.
}

\subsection{Posterior correlational ICOMP analogues}{
Let
\deqn{R_{\mathrm{post}}=\mathrm{Cor}_{\theta\mid y}(\theta)}
be the correlation matrix estimated from the posterior parameter draws. The
correlational complexity used by this function is
\deqn{C_R(R_{\mathrm{post}})
=-\frac{1}{2}\log|R_{\mathrm{post}}|.}

This complexity is zero when the posterior parameters are uncorrelated and
increases as their joint posterior approaches linear dependence. Because it
uses a correlation rather than a covariance matrix, it is invariant to
separate shifts, changes of units, rescalings, and sign reversals of the
parameters. It is not invariant to arbitrary nonlinear transformations or to
transformations that mix parameters. It also discards differences in the
marginal posterior variances and measures only posterior interdependence.

The function returns the following posterior analogues:
\deqn{\mathrm{ICOMP}=\bar D+2C_R,}
\deqn{\mathrm{ICOMP}_{\mathrm{PEU}}=\bar D+d+2C_R,}
\deqn{\mathrm{ICOMP}_{\mathrm{PEU\_LN}}
=\bar D+d+2\log(N)C_R,}
\deqn{\mathrm{CICOMP}
=\bar D+d\{\log(N)+1\}+2C_R.}

The element \code{ICOMP_PEU} is recommended as the main posterior ICOMP
summary because it combines posterior mean fit, an explicit dimensional
penalty, and a penalty for posterior parameter dependence. The basic
\code{ICOMP} element has no direct penalty for adding posteriorly independent
parameters. \code{ICOMP_PEU_LN} increases the dependence penalty with sample
size, while \code{CICOMP} combines a CAIC-like dimensional penalty with
correlational complexity.

These quantities are not the exact published ICOMP criteria. The original
forms use minus twice the maximized log-likelihood and a complexity measure
of the estimated inverse Fisher information matrix. In \code{postIC()},
\eqn{\bar D} replaces the maximized plug-in deviance and
\eqn{R_{\mathrm{post}}} replaces the correlation matrix derived from inverse
Fisher information. The posterior distribution and its prior therefore affect
the complexity term. The original asymptotic or consistency properties should
not be assumed to transfer automatically to these posterior analogues.

When a posterior column has zero variance, the ICOMP quantities are returned
as \code{NA}, because a correlation matrix cannot be defined. When the
posterior correlation matrix is singular, \eqn{|R_{\mathrm{post}}|=0} and the
correlational complexity diverges; \code{C_R} and all four ICOMP quantities
are therefore returned as \code{Inf}.
}
}
\references{
Bozdogan, H. (1990). On the information-based measure of covariance
complexity and its application to the evaluation of multivariate linear
models. \emph{Communications in Statistics---Theory and Methods, 19},
221--278.

Bozdogan, H. (2000). Akaike's information criterion and recent developments
in information complexity. \emph{Journal of Mathematical Psychology, 44},
62--91.

Geisser, S., & Eddy, W. F. (1979). A predictive approach to model selection.
\emph{Journal of the American Statistical Association, 74}, 153--160.

Gelfand, A. E., Dey, D. K., & Chang, H. (1992). Model determination using
predictive distributions with implementation via sampling-based methods. In
J. M. Bernardo et al. (Eds.), \emph{Bayesian Statistics 4} (pp. 147--167).
Oxford University Press.

Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
information criteria for Bayesian models. \emph{Statistics and Computing, 24},
997--1016.

Kartal-Koc, E., & Bozdogan, H. (2015). Model selection in multivariate
adaptive regression splines using information complexity as the fitness
function. \emph{Machine Learning, 101}, 35--58.

Pamukcu, E., Bozdogan, H., & Calik, S. (2015). A novel hybrid dimension
reduction technique for undersized high dimensional gene expression data sets
using information complexity criterion for cancer classification.
\emph{Computational and Mathematical Methods in Medicine, 2015}, 370640.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC.
\emph{Statistics and Computing, 27}, 1413--1432.

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., & Gabry, J. (2024). Pareto
smoothed importance sampling. \emph{Journal of Machine Learning Research,
25}(72), 1--58.

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
widely applicable information criterion in singular learning theory.
\emph{Journal of Machine Learning Research, 11}, 3571--3594.

Watanabe, S. (2013). A widely applicable Bayesian information criterion.
\emph{Journal of Machine Learning Research, 14}, 867--897.

Xiao, X., & Rabe-Hesketh, S. (2026). A parameterization-invariant DIC.
\emph{arXiv preprint arXiv:2605.27844}.
}
\examples{
## Normal model with an analytically sampled posterior for the mean
set.seed(123)
y <- stats::rnorm(30, mean = 0.4, sd = 1)

## Prior: mu ~ N(0, 10^2); sampling SD is known and equal to 1
post_var <- 1 / (length(y) + 1 / 100)
post_mean <- post_var * sum(y)
mu_draws <- stats::rnorm(1000, post_mean, sqrt(post_var))

log_lik <- t(vapply(
  mu_draws,
  function(mu) stats::dnorm(y, mean = mu, sd = 1, log = TRUE),
  numeric(length(y))
))

result <- postIC(log_lik, posterior = matrix(mu_draws, ncol = 1))
result$WAIC
result$LPML
result$DIC_i
result$ICOMP_PEU

## The dependence penalty is nonzero with correlated posterior parameters.
posterior_two <- cbind(
  theta1 = mu_draws,
  theta2 = 0.8 * mu_draws + stats::rnorm(length(mu_draws), sd = 0.1)
)
result_two <- postIC(log_lik, posterior = posterior_two)
result_two$C_R
result_two$ICOMP
result_two$ICOMP_PEU
result_two$ICOMP_PEU_LN
result_two$CICOMP

## PSIS-LOO fields are populated when the optional loo package is installed.
result$LOOIC
}
