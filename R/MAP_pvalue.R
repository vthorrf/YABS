MAP_pvalue <- function(x, reference.value=0, MAP.ID=NULL) {
  if(!is.null(MAP.ID)) {
    MAP <- x[MAP.ID]
    from_to <- sort(c(reference.value, MAP))
    order_test <- order(c(reference.value, MAP))
    dens <- density(x, from=from_to[1], to=from_to[2], n=2)$y
    return(dens[order_test[1]]/dens[order_test[2]])
  } else {
    MAP <- mean(x)
    from_to <- sort(c(reference.value, MAP))
    order_test <- order(c(reference.value, MAP))
    dens <- density(x, from=from_to[1], to=from_to[2], n=2)$y
    return(dens[order_test[1]]/dens[order_test[2]])
  }
}
