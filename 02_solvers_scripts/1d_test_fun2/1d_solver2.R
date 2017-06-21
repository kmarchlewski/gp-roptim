
ObjectFunction <- function (x, vectorized = FALSE) {
  if (vectorized == TRUE) {
    a1 <- 3.0 * (x - 0.75) ^ 2
    a2 <- -3.5 * exp(-(25 * (x - 0.3)) ^ 2)
    a3 <- -2.0 * exp(-(4 * (x - 0.9)) ^ 2)
    a4 <- -0.1 / (x - 1.025)
    
    val <- a1 + a2 + a3 + a4
    state <- matrix("converged", nrow(x), 1)
  } else {
    a1 <- 3.0 * (x - 0.75) ^ 2
    a2 <- -3.5 * exp(-(25 * (x - 0.3)) ^ 2)
    a3 <- -2.0 * exp(-(4 * (x - 0.9)) ^ 2)
    a4 <- -0.1 / (x - 1.015)
    
    val <- a1 + a2 + a3 + a4
    state <- "converged"
  }

  return (list(val = val, state = state))
}
