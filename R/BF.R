BF <- function(x, reference.value=0, prior.density, ...) {
  density(x, from=reference.value)$y[1] / prior.density(reference.value, ...)
}