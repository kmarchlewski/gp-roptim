
ObjectFunction = function (x, vectorized = FALSE) {
  if (vectorized == TRUE) {
    a1 <- 3 * (x - 0.75) ^ 2
    a2 <- -3 * exp(-(12 * (x - 0.3)) ^ 2)
    a3 <- -1.6 * exp(-(4 * (x - 0.9)) ^ 2)
    a4 <- -2 * exp(-(7 * (x - 1.6)) ^ 2)
    val <- a1 + a2 + a3 + a4
    state <- matrix("converged", nrow(x), 1)
  } else {
    a1 <- 3 * (x - 0.75) ^ 2
    a2 <- -3 * exp(-(12 * (x - 0.3)) ^ 2)
    a3 <- -1.6 * exp(-(4 * (x - 0.9)) ^ 2)
    a4 <- -2 * exp(-(7 * (x - 1.6)) ^ 2)
    val <- a1 + a2 + a3 + a4
    state <- "converged"
  }

  return (list(val = val, state = state))
}
