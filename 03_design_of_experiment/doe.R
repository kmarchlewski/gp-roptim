
library(lhs)

DoE <- function (N, L, x.min, x.max, add.points = NULL,
                 Function = ObjectFunction,
                 method = "geneticLHS", export = TRUE) {
  if (method == "geneticLHS") {
    
    if (is.null(add.points)) {
      X.add <- NULL
      N.add <- 0
    } else if (is.matrix(add.points) & ncol(add.points) == L) {
      X.add <- add.points
      N.add <- nrow(add.points)
    } else {
      stop("Wrong \"add.points\" value.")
    }
    
    X <- t(x.min + t(geneticLHS(N, L)) * (x.max - x.min))
    X <- rbind(X, X.add)
    f.aug <- Function(X, vectorized = TRUE)
    f <- matrix(f.aug$val, N + N.add, 1)
    state <- f.aug$state
    
  } else if (method == "maximinLHS") {
    
    if (is.null(add.points)) {
      X.add <- NULL
      N.add <- 0
    } else if (is.matrix(add.points) & ncol(add.points) == L) {
      X.add <- add.points
      N.add <- nrow(add.points)
    } else {
      stop("Wrong \"add.points\" value.")
    }
    
    X <- t(x.min + t(maximinLHS(N, L)) * (x.max - x.min))
    X <- rbind(X, X.add)
    f.aug <- Function(X, vectorized = TRUE)
    f <- matrix(f.aug$val, N + N.add, 1)
    state <- f.aug$state
    
  } else if (method == "FullFactorial") {
    
    n <- N ^ (1 / L)
    
    if ((n %% 1) == 0) {
      n <- n
    } else if ((n %% 1) != 0) {
      n <- ceiling(n)
      N <- n ^ L
      warning(paste0("Demanded number of points is not adequate. ",
              N, " points used."))
    }

    if (is.null(add.points)) {
      X.add <- NULL
      N.add <- 0
    } else if (is.matrix(add.points) & ncol(add.points) == L) {
      X.add <- add.points
      N.add <- nrow(add.points)
    } else {
      stop("Wrong \"add.points\" value.")
    }
    
    X <- expand.grid(lapply(1:L, function (x) seq(x.min[x], x.max[x], len = n)))
    X <- as.matrix(X)
    X <- rbind(X, X.add)
    f.aug <- Function(X, vectorized = TRUE)
    f <- matrix(f.aug$val, N + N.add, 1)
    state <- f.aug$state
  } else {
    stop("Wrong \"method\" value.")
  }
  
  if (export == TRUE) {
    write.table(cbind(X, f, state), file = "doe_data.csv",
                col.names = c(sapply(1:L, function (i) paste0("x_", i)),
                              "f", "state"),
                row.names = FALSE, sep = ";", dec = ".", quote = FALSE)
    return (list(coord = X, val = f, state = state))
  } else if (export == FALSE) {
    return (list(coord = X, val = f, state = state))
  } else {
    return (list(coord = X, val = f, state = state))
    stop("Wrong \"export\" value.")
  }
}
