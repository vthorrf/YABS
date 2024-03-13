dmnorm    <- function(x, mu=NULL, sigma, log=FALSE){
  if(is.null(mu)) mu <- colMeans(x)
  d <- ncol(x)
  Mu <- matrix(mu, byrow=T, nrow=nrow(x), ncol=ncol(x))
  z <- rowSums({{x - Mu} %*% solve(sigma)} * {x - Mu})
  logDensity <- {{-d/2}*log(2*pi)} +
    {-.5 * log(det(sigma))} +
    {-.5 * z}
  if(log) {
    return(logDensity)
  } else {
    return(exp(logDensity))
  }
}
