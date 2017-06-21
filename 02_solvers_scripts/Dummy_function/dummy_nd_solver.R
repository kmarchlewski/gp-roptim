
ObjectFunction <- function (x, vectorized = FALSE, fun = dummy) {
  if (vectorized == TRUE) {
    val <- matrix(apply(x, 1, fun), nrow(x), 1)
    state <- matrix("converged", nrow(x), 1)
  } else {
    val <- fun(x)
    state <- "converged"
  }

  return (list(val = val, state = state))
}

dummy <- function (x) {
  return (450)
}

