
ObjectFunction = function (x, vectorized = FALSE) {
  l.lim.uc <- 0.0
  u.lim.uc <- 0.4
  
  if (vectorized == TRUE) {
    val <- matrix(NA, nrow(x), 1)
    state <- matrix(NA, nrow(x), 1)
    
    for (i in 1:nrow(x)) {
     x.row <- x[i, , drop = FALSE]
     
      if (all(x.row > as.matrix(l.lim.uc) &
          x.row < as.matrix(u.lim.uc))) {
        val[i,] <- as.matrix(-10)
        state[i,] <- as.matrix("unconverged")
      } else {
        a1 <- 3 * (x.row - 0.75) ^ 2
        a2 <- -3 * exp(-(12 * (x.row - 0.3)) ^ 2)
        a3 <- -1.6 * exp(-(4 * (x.row - 0.9)) ^ 2)
        a4 <- -2 * exp(-(7 * (x.row - 1.6)) ^ 2)
        val[i,] <- a1 + a2 + a3 + a4
        state[i,] <- as.matrix("converged")
      }
    }
    
  } else {
    
    if (all(x > l.lim.uc & x < u.lim.uc)) {
      val <- -10
      state <- "unconverged"
    } else {
      a1 <- 3 * (x - 0.75) ^ 2
      a2 <- -3 * exp(-(12 * (x - 0.3)) ^ 2)
      a3 <- -1.6 * exp(-(4 * (x - 0.9)) ^ 2)
      a4 <- -2 * exp(-(7 * (x - 1.6)) ^ 2)
      val <- a1 + a2 + a3 + a4
      state <- "converged" 
    }
    
  }

  return (list(val = val, state = state))
}
