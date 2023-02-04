\name{harm}
\alias{harm}
\title{Random-walk Metropolis Sampler}
\description{
Internal function. Samples new values from the logposterior of the user provided model using the Random-walk Metropolis MCMC algorithm.
}
\usage{
harm(Model, Data, Iterations, Status, Thinning, ACC, DevianceMat,
     LIV, Monitor, Mo0, samples, PPD, Adapt, Sigma)
}
\arguments{
   \item{Model}{The user provided model.}
   \item{Data}{A named list of the data objects required by the model. Should also include the Parameter-Generating Function (PGF) if initial values for the parameters are not provided.}
   \item{Iterations}{Total number of iterations (including burn-in).}
   \item{Status}{The number of iterations to run before updating the status of the sampler.}
   \item{Thinning}{Thinning rate. Must be a positive integer.}
   \item{ACC}{The initial acceptance rate.}
   \item{DevianceMat}{The (empty) matrix of deviances.}
   \item{LIV}{The length of initial values (i.e., the number of parameters).}
   \item{Monitor}{The (empty) matrix of variables/parameters to be be monitored (other than the model's parameters).}
   \item{Mo0}{A list with the output of Model for the Initial Values.}
   \item{samples}{The (empty) matrix of thinned values.}
   \item{PPD}{The (empty) matrix of posterior criterion varible samples.}
   \item{Adapt}{Total number of adaptive iterations.}
   \item{Sigma}{The initial covariance matrix between parameters.}
}
\value{
A list containing: (i) "Acceptance", the achieved acceptance rate; (ii) "Dev", the posterior distribution of the deviance; (iii) "Mon", the posterior distribution of the monitored variables/parameters; and (iv) "thinned", the posterior distribution of all parameters of the model.
}
\examples{
# This function should be used internally only. See the `MCMC` function for more details.
}
