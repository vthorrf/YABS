\name{gcharm}
\alias{gcharm}
\title{Gradient-Corrected Hit-and-Run Metropolis Sampler}
\description{
Internal function. Samples new values from the logposterior of the user provided model using the Gradient-Corrected Hit-and-Run Metropolis MCMC algorithm.
}
\usage{
gcharm(Model, Data, Iterations, Status, Thinning, Acceptance, Dev, h,
       LIV, Mon, Mo0, thinned)
}
\arguments{
   \item{Model}{The user provided model.}
   \item{Data}{A named list of the data objects required by the model. Should also include the Parameter-Generating Function (PGF) if initial values for the parameters are not provided.}
   \item{Iterations}{Total number of iterations (including burn-in).}
   \item{Status}{The number of iterations to run before updating the status of the sampler.}
   \item{Thinning}{Thinning rate. Must be a positive integer.}
   \item{Acceptance}{The initial acceptance rate.}
   \item{Dev}{The (empty) matrix of deviances.}
   \item{h}{The interval of numeric differencing.}
   \item{LIV}{The length of initial values (i.e., the number of parameters).}
   \item{Mon}{The (empty) matrix of variables/parameters to be be monitored (other than the model's parameters).}
   \item{Mo0}{A list with the output of Model for the Initial Values.}
   \item{thinned}{The (empty) matrix of thinned values.}
}
\value{
A list containing: (i) "Acceptance", the achieved acceptance rate; (ii) "Dev", the posterior distribution of the deviance; (iii) "Mon", the posterior distribution of the monitored variables/parameters; and (iv) "thinned", the posterior distribution of all parameters of the model.
}
\examples{
# This function should be used internally only. See the `MCMC` function for more details.
}