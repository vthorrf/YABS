\name{harmwg}
\alias{harmwg}
\title{Metropolis-Within Gibbs Sampler}
\description{
Internal function. Samples new values from the logposterior of the user provided model using the Metropolis-Within Gibbs MCMC algorithm.
}
\usage{
harmwg(Model, Data, Iterations, Status, Thinning, Acceptance, Dev, LIV,
       Mon, Mo0, thinned, adapt, epsilon)
}
\arguments{
   \item{Model}{The user provided model.}
   \item{Data}{A named list of the data objects required by the model. Should also include the Parameter-Generating Function (PGF) if initial values for the parameters are not provided.}
   \item{Iterations}{Total number of iterations (including burn-in).}
   \item{Status}{The number of iterations to run before updating the status of the sampler.}
   \item{Thinning}{Thinning rate. Must be a positive integer.}
   \item{Acceptance}{The initial acceptance rate.}
   \item{Dev}{The (empty) matrix of deviances.}
   \item{LIV}{The length of initial values (i.e., the number of parameters).}
   \item{Mon}{The (empty) matrix of variables/parameters to be be monitored (other than the model's parameters).}
   \item{Mo0}{A list with the output of Model for the Initial Values.}
   \item{thinned}{The (empty) matrix of thinned values.}
   \item{adapt}{Total number of adaptive iterations.}
   \item{epsilon}{The initial covariance matrix between parameters.}
}
\value{
A list containing: (i) "Acceptance", the achieved acceptance rate; (ii) "Dev", the posterior distribution of the deviance; (iii) "Mon", the posterior distribution of the monitored variables/parameters; and (iv) "thinned", the posterior distribution of all parameters of the model.
}
\examples{
# This function should be used internally only. See the `MCMC` function for more details.
}
