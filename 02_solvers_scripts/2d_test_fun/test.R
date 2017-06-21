
#
# File is used to test new 2d functions.
#

library(rgl)
library(fields)

L <- 2
n <- 50

x.min <- c(0, 0)
x.max <- c(2, 2)

x.1 <- seq(x.min[1], x.max[1], len = n)
x.2 <- seq(x.min[2], x.max[2], len = n)
X <- as.matrix(expand.grid(x.1, x.2))

ObjectFunction <- function (x, a = 10, b = 3, c = 2) {
  alpha <- ((a - b) / 2) * x[2] + b
  beta <- -(c / 2) * x[2] + c
  val <- alpha * (x[1] - 1) ^ 2 + beta
  
  return (list(val = val))
}

ObjectFunction2 <- function (x, a = 3, b = 2, c = 1) {
  # alpha <- -0.5 * x[2] + 1
  # beta <- 0.5 * x[2]
  # alpha <- exp(-x[2])
  # beta <- 1 - exp(-x[2])
  alpha <- 0.25 * (x[2] - 2) ^ 2
  beta <- 1 - alpha
  
#   val <- alpha * (b * (x[1] - 1) ^ 2 + c) +
#          beta * (a * sqrt(abs(x[1] - 1)))
  val <- alpha * (b * (x[1] - 1) ^ 4 + c) +
    beta * (5 * (x[1] - 1) ^ 2)
  
  return (list(val = val))
}

X.1 <- matrix(x.1, n, n)
X.2 <- t(matrix(x.2, n, n))
Y <- matrix(apply(X, 1, function (x) ObjectFunction(x)$val), n, n)
Y.2 <- matrix(apply(X, 1, function (x) ObjectFunction2(x)$val), n, n)

persp3d(X.1, X.2, Y.2,
        col = "blue", alpha = 0.6, expand = 0.75,
        xlab = "", ylab = "", zlab = "")

image.plot(Y.2)
contour(Y.2, add = TRUE)
