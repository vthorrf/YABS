\name{nuts}
\alias{nuts}
\title{No-U-Turn Sampler}
\description{
Internal function. Samples new values from the log-posterior of the user provided model using the No-U-Turn Sampler (NUTS), an adaptive Hamiltonian Monte Carlo (HMC) algorithm that automatically tunes trajectory length via recursive tree building and step-size adaptation.
}
\usage{
nuts(Model, Data, Iterations, Status, InitialValues, Thinning,
     thinned, postpred, Dev, Mon, Mo0, step_size = 0.01,
     max_depth = 10, Adapt = 1000, target = 0.65)
}
\arguments{
   \item{Model}{The user provided model.}
   \item{Data}{A named list of the data objects required by the model. Should also include the Parameter-Generating Function (PGF) if initial values for the parameters are not provided.}
   \item{Iterations}{Total number of iterations (including burn-in).}
   \item{Status}{The number of iterations to run before updating the status of the sampler.}
   \item{InitialValues}{A numeric vector of initial parameter values.}
   \item{Thinning}{Thinning rate. Must be a positive integer.}
   \item{thinned}{The (empty) matrix of thinned parameter samples.}
   \item{postpred}{The (empty) matrix of posterior predictive samples.}
   \item{Dev}{The (empty) matrix of deviances.}
   \item{Mon}{The (empty) matrix of monitored variables/parameters (other than the model's parameters).}
   \item{Mo0}{A list with the output of \code{Model} evaluated at the initial values.}
   \item{step_size}{Initial step size for the leapfrog integrator.}
   \item{max_depth}{Maximum tree depth for trajectory expansion.}
   \item{Adapt}{Number of iterations used for step-size adaptation via dual averaging.}
   \item{target}{Target acceptance probability used in dual averaging adaptation.}
}
\value{
A list containing:
\itemize{
  \item{\code{Acceptance}}{The average acceptance probability across iterations.}
  \item{\code{Dev}}{The posterior distribution of the deviance.}
  \item{\code{Mon}}{The posterior distribution of the monitored variables/parameters.}
  \item{\code{thinned}}{The posterior distribution of all model parameters.}
  \item{\code{postpred}}{Posterior predictive samples.}
  \item{\code{step_size}}{Final adapted step size used by the sampler.}
}
}
\examples{
# This function should be used internally only. See the `MCMC` function for more details.
}