
ObjectFunction <- function (x, vectorized = FALSE, fun = RM2d2) {
  if (vectorized == TRUE) {
    val <- matrix(apply(x, 1, fun), nrow(x), 1)
    state <- matrix("converged", nrow(x), 1)
  } else {
    val <- fun(x)
    state <- "converged"
  }
  
  return (list(val = val, state = state))
}

RM2d <- function (x, a = 3, b = 0.75, c = 1) {
  # alpha <- -0.5 * x[2] + 1
  # beta <- 0.5 * x[2]
  alpha <- exp(-x[2])
  beta <- 1 - exp(-x[2])
  
  val <- alpha * (b * (x[1] - 1) ^ 2 + c) +
         beta * (a * sqrt(abs(x[1] - 1)))
  
  return (val)
}

RM2d2 <- function (x, a = 5, b = 2, c = 1) {
  alpha <- 0.25 * (x[2] - 2) ^ 2
  beta <- 1 - alpha

  val <- alpha * (b * (x[1] - 1) ^ 4 + c) +
         beta * (a * (x[1] - 1) ^ 2)
  
  return (val)
}
