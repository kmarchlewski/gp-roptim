
# Sprawdzanie minimalnej liczby pkt. niezbednych do poprawnego dopasowania
# modelu powinno byc w wielu fun.
# Punktow powinno byc co najmniej (bez uwzgl. wyrazow mieszanych):
# poly.deg * L + 1 + L + 2
# Trend wielomianowy powinien uwzgledniac wyrazy mieszane (np. x_1 * x_2)
# poprawic ladowanie gaussAlgebra tak zeby komunikaty byly wyswietlane na czerwono
# zmienic listy na data.frame tam gdzie mamy same liczby
# FunctionOptim sprawdza czy w wyniku optymalizacji nie wychodzimy poza obszar
# poszukiwan optimum

TruncGvec <- function (g.vec, lim = .Machine$double.eps) {
  library(abind)
  
  tab <- g.vec@vec[[1]]@tab
  d <- dim(tab)
  sel.poly <- apply(matrix(is.infinite(tab[, 1, ]), d[1], d[3]), 2, any)
  tab.g <- tab[, , !sel.poly, drop = FALSE]
  tab.p <- tab[, , sel.poly, drop = FALSE]
  
  d <- dim(tab.g)
  l <- 1:(d[2] - 2) - 1
  # zmienione na coeff bylo a.max
  coeff <- (l / exp(1)) ^ (l / 2)
  sigmas <- tab.g[, rep(1, d[2] - 2), , drop = FALSE]
  
  # a.max <- sapply(1:d[3], function (i) {
  #                           t(t(matrix(sigmas[, , i], d[1], d[2] - 2)) ^ l)
  #                         }, simplify = "array")
  
  a.max <- array(sapply(1:d[3], function (i) {
                                  t(coeff * t(matrix(sigmas[, , i], d[1], d[2] - 2)) ^ l)
                                }), dim = c(d[1], d[2] - 2, d[3]))
  
  gAlg.max <- sapply(1:d[3], function (i) {
                               prod(rowSums(matrix(abs(tab.g[, 3:d[2], i]) * a.max[, , i],
                                                   d[1], d[2] - 2)))
                             }, simplify = "array")
  
  sel.gAlg <- (gAlg.max < lim)
  tab.g.trunc <- tab.g[, , !sel.gAlg, drop = FALSE]
  tab.trunc <- abind(tab.g.trunc, tab.p, along = 3)
  
  g.vec.trunc <- as.gvector(new("gAlg", tab = tab.trunc))
  return (list(gVec = g.vec.trunc))
}

NumZeroSubstitute <- function (num, lim = .Machine$double.eps) {
  N <- nrow(num)
  
  if (any(num <= -lim)) {
    n <- which(num <= -lim)
    
    PrintNumber <- function (i) {
      warning (paste0("Number ", sprintf("%.10f", num[i, ]),
                      " is smaller then limit value ", lim))
    }
    
    invisible(sapply(n, PrintNumber))
  }
  
  return (ifelse(num < lim, rep(0, N), num))
}

ErrorGauss <- function (st.dev) {
  L <- length(st.dev)
  ErrAlg <- V(Gauss(st.dev))
  
  ErrNum <- function (delta.x) {
    N <- nrow(delta.x)
    coeff <- 1 / (sqrt(((2 * pi) ^ L)) * prod(st.dev))
    val <- coeff * exp(- 0.5 * rowSums(t(t(delta.x) / st.dev) ^ 2))
    return (matrix(val, N, 1))
  }
  
  return (list(ErrAlg = ErrAlg, ErrNum = ErrNum, ErrName = "Gauss"))
}

ErrorBeta22 <- function (delta) {
  library(abind)
  
  L <- length(delta)
  M <- 5
  
  sigma <- c(0.0607890643500889, 0.100500633192404, 0.14988703175021,
             0.0607890643500889, 0.100500633192404)
  sigma <- as.matrix(expand.grid(lapply(1:L, function (i) sigma)))
  sigma <- t(2 * delta * t(sigma)) ^ 2
   
  mu <- c(0.3783410796803280, 0.2391008580028830, 0,
          -0.378341079680328, -0.239100858002883)
  mu <- as.matrix(expand.grid(lapply(1:L, function (i) mu)))
  mu <- t(2 * delta * t(mu))
  
  a <- c(0.316140901699594, 0.738565413648224, 1.42555841080486,
         0.316140901699593, 0.738565413648225)
  a <- as.matrix(expand.grid(lapply(1:L, function (i) a)))
  a <- t(t(a) / (2 * delta))
  
  tab <- abind(aperm(array(sigma, c(M ^ L, L, 1)), c(2, 3, 1)),
               aperm(array(mu, c(M ^ L, L, 1)), c(2, 3, 1)),
               aperm(array(a, c(M ^ L, L, 1)), c(2, 3, 1)),
               along = 2)
  
  ErrAlg <-  as.gvector(new("gAlg", tab = tab))
  ErrAlg.t <- TruncGvec(ErrAlg)$gVec
  #
  if (any(dim(ErrAlg.t@vec[[1]]@tab) != dim(ErrAlg@vec[[1]]@tab))) {
    print(dim(ErrAlg.t@vec[[1]]@tab))
  }
  #
  
  ErrNum <- function (x) {
    N <- nrow(x)
    L <- ncol(x)
    y <- 0.5 * t(t(x) / delta + 1)

    val <- matrix(sapply(1:L, function (i) {
                                dbeta(y[, i, drop = FALSE], 2, 2) / (2 * delta[i])
                              }), N, L)
    val <- matrix(apply(val, 1, prod), N, 1)
    
    return (val)
  }
  
  return (list(ErrAlg = ErrAlg.t, ErrNum = ErrNum, ErrName = "Beta-22"))
}

### temp beta dist. check
#library(rgl)
#delta <- c(0.5, 1.5)
#x <- matrix(NA, 50, 2)
#x[, 1] <- seq(-1, 1, len = 50)
#x[, 2] <- seq(-2, 2, len = 50)
#
#val.plot <- ErrorBeta22(delta)$ErrNum(expand.grid(x[, 1], x[, 2]))
#val.plot.2 <- calc(ErrorBeta22(delta)$ErrAlg, as.matrix(expand.grid(x[, 1], x[, 2])))
#range(val.plot - val.plot.2)
#
#persp3d(matrix(x[, 1], nrow(x), nrow(x)), t(matrix(x[, 2], nrow(x), nrow(x))),
#        matrix(val.plot, nrow(x), nrow(x)), col = "blue")
#persp3d(matrix(x[, 1], nrow(x), nrow(x)), t(matrix(x[, 2], nrow(x), nrow(x))),
#        matrix(val.plot.2, nrow(x), nrow(x)), col = "red", add = TRUE)
#
#delta = 1
#x = matrix(seq(-3, 3, len = 100), 100, 1)
#val.plot <- ErrorBeta22(delta)$ErrNum(x)
#val.plot.2 <- calc(ErrorBeta22(delta)$ErrAlg, x)
#plot(x, val.plot, type = "l")
#lines(x, val.plot.2, col ="red")
#lines(x, ErrorGauss(delta)$ErrNum(x), col="blue")
###

CovarianceMatern <- function (theta) {
  L <- length(theta)
  M <- 7

  sigma <- c(0.471944930499580, 0.286464200868313, 0.145940602273492,
             1.812927886021040, 0.706857808096680, 1.350991423953320,
             0.995758253911648)
  sigma <- matrix(sigma, M, 1)
  
  a <- c(0.0886027859347606, 0.0132081206640146, 0.000774086748169311,
         0.0238626993273815, 0.2776228537222050, 0.205267100887735000,
         0.3906602990061930)
  a <- matrix(a, M, 1)
  
  CovAlg <- V(sapply(1:M, function (i) {
                            prod(theta) * Gauss(rep(sigma[i, 1], L) * theta)
                          }
             )) %*% (a * (sqrt(2 * pi) * sigma) ^ L)
  
  CovNum <- function (delta.x) {
    N <- nrow(delta.x)
    R <- matrix(sqrt(rowSums(t(t(delta.x) / theta) ^ 2 )), N, 1)
    c.1 <- sqrt(5) * abs(R)
    c.2 <- 5 / 3 * R ^ 2
    
    coeff <- 1 + c.1 + c.2
    val <- coeff * exp(-sqrt(5) * abs(R))
    
    return (val)
  }
  
  return (list(CovAlg = CovAlg, CovNum = CovNum, CovName = "Matern 5/2"))
}

CovarianceGauss <- function (theta) {
  L <- length(theta)
  coeff <- sqrt(((2 * pi) ^ L) * prod(theta ^ 2))
  CovAlg <- V(coeff * Gauss(theta))
  
  CovNum <- function (delta.x) {
    N <- nrow(delta.x)
    val <- exp(- 0.5 * rowSums(t(t(delta.x) / theta) ^ 2))
    return (matrix(val, N, 1))
  }
 
  return (list(CovAlg = CovAlg, CovNum = CovNum, CovName = "Gauss"))
}

CovarianceMatrix <- function (X.1, X.2, theta,
                              CovarianceFunction = krig.settings$CovarianceType) {
  L <- length(theta)
  n.1 <- nrow(X.1)
  n.2 <- nrow(X.2)
  V <- expand.grid(i = 1:n.1, j = 1:n.2)
  DIFFS <- X.1[V$i, , drop = FALSE] - X.2[V$j, , drop = FALSE]
  C <- CovarianceFunction(theta)$CovAlg
  
  return (matrix(calc(C, DIFFS), n.1, n.2))
}

CovarianceMatrixNum <- function (X.1, X.2, theta,
                                 CovarianceFunction = krig.settings$CovarianceType) {
  L <- length(theta)
  n.1 <- nrow(X.1)
  n.2 <- nrow(X.2)
  V <- expand.grid(i = 1:n.1, j = 1:n.2)
  DIFFS <- X.1[V$i, , drop = FALSE] - X.2[V$j, , drop = FALSE]
  C <- CovarianceFunction(theta)$CovNum
  
  return (matrix(C(DIFFS), n.1, n.2))
}

Poly <- function (L, poly.deg = krig.settings$poly.deg) {
  
  if (poly.deg == 0) {
    H.0 <- const.gAlg(1, L)
    
    return (V(H.0))
  } else {
    Term <- function (x, L) {
      H.1 <- 1
      n <- ceiling(x / L)
      x <- x - L * (n - 1)
      for (i in 1:n) {
        H.1 <- H.1 * linear.gAlg(x, L)
      }
      
      return (H.1)
    }
    
    H.0 <- const.gAlg(1, L)
    H.deg <- lapply(1:(L * poly.deg), function (x) Term(x, L))
    
    return (V(H.0, H.deg))
  }
}

PolyNum <- function (X, poly.deg = krig.settings$poly.deg) {
  
  if (poly.deg == 0) {
    H <- matrix(1, nrow(X), 1)
  } else {
    X <- lapply(1:poly.deg, function (i) X ^ i)
    X <- do.call(cbind, X)
    H <- cbind(1, X) 
  }
  
  return (H)
}

SafeInvert <- function (M, eigen.rat = krig.settings$eigen.rat) {
  M.inv <- try(solve(M), silent = TRUE)
  N.mod <- dim(M)[1]
  val <- NULL
  
  if (class(M.inv) == "try-error") {
    warning("\"solve()\" failed to invert the matrix. \"eigen()\" used.",
            noBreaks. = TRUE)
    M.mod <- eigen(M, symmetric = TRUE)
    sel <- M.mod$val > (max(M.mod$val) / eigen.rat)
    N.mod <- sum(sel)
    val <- M.mod$val[sel]
    M.inv <- M.mod$vec[, sel, drop = FALSE] %*%
             ((1 / M.mod$val[sel]) * diag(nrow = N.mod)) %*%
             t(M.mod$vec[, sel, drop = FALSE])
  }
  
  return (list(matrix = M.inv, dim = N.mod, val = val))
}

KrigingParOpt <- function (X, f, settings = krig.settings) {
  N <- nrow(X)
  L <- ncol(X)
  par.fun.type <- settings$par.fun.type
  theta.min <- settings$theta.min
  theta.max <- settings$theta.max
  alpha.min <- settings$alpha.min
  alpha.max <- settings$alpha.max
  opt.search.type <- settings$opt.search.type
  par.min <- c(rep(-1, L), 0)
  par.max <- c(rep(4, L), 4)
  
  if (par.fun.type == "LogLikelihood") {
    ParFun = LogLikelihood
    ParFunOpt = LogLikelihoodOpt
  } else if (par.fun.type == "LOOCV") {
    ParFun = LOOCV
    ParFunOpt = LOOCVOpt
  } else {
    stop("Wrong \"par.fun.type\" value.")
  }
  
  if (opt.search.type == "genetic") {
    gen.enhance <- settings$gen.enhance
    gen.size <- settings$gen.size
    pop.size <- settings$pop.size
    
    result.gen <- nsga2(function (x) {
                          ParFunOpt(X, f, x, par.min, par.max, settings)$val
                        },
                        L + 1, 1, generations = gen.size, popsize = pop.size,    
                        lower.bounds = par.min,
                        upper.bounds = par.max)
    sel.gen <- which.min(result.gen$val)
    
    if (gen.enhance == TRUE) {
      result.grad <- optim(result.gen$par[sel.gen, ],
                           function (x) {
                             ParFunOpt(X, f, x, par.min, par.max, settings)$val
                           },
                           method = "L-BFGS-B",
                           lower = par.min,
                           upper = par.max)
      
      if (result.grad$convergence == 0) {
        theta <- ScaleTheta(result.grad$par[1:L], theta.min, theta.max,
                            par.min[1:L], par.max[1:L])
        alpha <- ScaleAlpha(result.grad$par[L + 1], alpha.min, alpha.max,
                            par.min[L + 1], par.max[L + 1])
      } else {
        warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                noBreaks. = TRUE)
        theta <- ScaleTheta(result.gen$par[sel.gen, 1:L], theta.min, theta.max,
                            par.min[1:L], par.max[1:L])
        alpha <- ScaleAlpha(result.gen$par[sel.gen, L + 1], alpha.min, alpha.max,
                            par.min[L + 1], par.max[L + 1])
      }
    } else if (gen.enhance == FALSE) {
      theta <- ScaleTheta(result.gen$par[sel.gen, 1:L], theta.min, theta.max,
                          par.min[1:L], par.max[1:L])
      alpha <- ScaleAlpha(result.gen$par[sel.gen, L + 1], alpha.min, alpha.max,
                          par.min[L + 1], par.max[L + 1])
    } else {
      stop("Wrong \"gen.enhance\" value.")
    }
    par.fun <- ParFun(X, f, theta, alpha, settings)
    
  } else if (opt.search.type == "gradient") {
    lhs.add.points <- settings$lhs.add.points
    lhs.point.number <- settings$lhs.point.number
    lhs.enhance <- settings$lhs.enhance
    
    if (is.null(lhs.add.points)) {
      par.vec <- maximinLHS(lhs.point.number, L + 1)
    } else if (is.matrix(lhs.add.points)) {
      
      if (all(signif(apply(lhs.add.points[, 1:L, drop = FALSE], 2, min),
                     SignifPlaces(theta.min) + 5) < theta.min) ||
          all(signif(apply(lhs.add.points[, 1:L, drop = FALSE], 2, max),
                     SignifPlaces(theta.max) + 5) > theta.max) ||
          signif(min(lhs.add.points[, L + 1]),
                 SignifPlaces(alpha.min) + 5) < alpha.min ||
          signif(max(lhs.add.points[, L + 1]),
                 SignifPlaces(alpha.max) + 5) > alpha.max) {

        warning ("Additional points weren't used. Points out of range.")
                       
        lhs.add.points <- matrix(c(theta.min + (theta.max - theta.min) / 2,
                                   alpha.min + (alpha.max - alpha.min) / 2),
                                 1, L + 1)
      }
      
      lhs.add.points <- cbind(t((t(lhs.add.points[, 1:L, drop = FALSE]) - theta.min) *
                               (1 / (theta.max - theta.min))),
                              (lhs.add.points[, L + 1] - alpha.min) *
                              (1 / (alpha.max - alpha.min)))

      par.vec <- rbind(maximinLHS(lhs.point.number, L + 1), lhs.add.points)
    } else {
      stop ("Wrong \"lhs.add.points\" value.")
    }
    
    par.vec <- t(par.min + t(par.vec) * (par.max - par.min))
    par.fun <- apply(par.vec, 1,
                         function (x) {
                           ParFunOpt(X, f, x, par.min, par.max, settings)$val
                         })
    
    if (lhs.enhance == TRUE) {
      lhs.enhance.number <- settings$lhs.enhance.number
      lhs.enhance.points <- ceiling(lhs.enhance.number * lhs.point.number)
      sel <- order(par.fun)[1:lhs.enhance.points]
      theta.enhance <- matrix(0, lhs.enhance.points, L)
      alpha.enhance <- rep(0, lhs.enhance.points)
      par.fun.val.enhance <- rep(0, lhs.enhance.points)
      
      for (i in 1:lhs.enhance.points) {
        result.grad <- optim(par.vec[sel[i], ],
                             function (x) {
                               ParFunOpt(X, f, x, par.min, par.max, settings)$val
                             },
                             method = "L-BFGS-B",
                             lower = par.min,
                             upper = par.max)
        
        if (result.grad$convergence == 0) {
          theta.enhance[i,] <- ScaleTheta(result.grad$par[1:L], theta.min, theta.max,
                                          par.min[1:L], par.max[1:L])
          alpha.enhance[i] <- ScaleAlpha(result.grad$par[L + 1], alpha.min, alpha.max,
                                         par.min[L + 1], par.max[L + 1])
          
          par.fun.val.enhance[i] <- ParFun(X, f,
                                           theta.enhance[i,], alpha.enhance[i],
                                           settings)$val
        } else {
          warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                  noBreaks. = TRUE)
          theta.enhance[i,] <- rep(NA, L)
          alpha.enhance[i] <- NA
          par.fun.val.enhance[i] <- NA
        }
      }
      if (sum(par.fun.val.enhance, na.rm = TRUE) == 0) {
        sel <- which.min(par.fun)
        theta <- ScaleTheta(par.vec[sel, 1:L], theta.min, theta.max,
                            par.min[1:L], par.max[1:L])
        alpha <- ScaleAlpha(par.vec[sel, L + 1], alpha.min, alpha.max,
                            par.min[L + 1], par.max[L + 1])
      } else {
        sel <- which.min(par.fun.val.enhance)
        theta <- theta.enhance[sel,]
        alpha <- alpha.enhance[sel]
      }
    } else if (lhs.enhance == FALSE) {
      sel <- which.min(par.fun)
      theta <- ScaleTheta(par.vec[sel, 1:L], theta.min, theta.max,
                          par.min[1:L], par.max[1:L])
      alpha <- ScaleAlpha(par.vec[sel, L + 1], alpha.min, alpha.max,
                          par.min[L + 1], par.max[L + 1])
    } else {
      stop("Wrong \"lhs.enhance\" value.")
    }
    par.fun <- ParFun(X, f, theta, alpha, settings)
  } else if (opt.search.type == "disabled") {

    theta <- theta.min + (theta.max - theta.min) / 2
    alpha <- alpha.max
    par.fun <- ParFun(X, f, theta, alpha, settings)
    
  } else {
    stop ("Wrong \"opt.search.type\" value.")
  }
  
  return (list(theta.opt = theta, alpha.opt = alpha, v.opt = par.fun$v,
               par.fun.opt = par.fun$val, type = par.fun.type))
}

LogLikelihood <- function (X, f, theta, alpha = 0.9999, settings = krig.settings) {
  N <- nrow(X)
  L <- ncol(X)
  poly.deg <- settings$poly.deg
  CovarianceFunction <- settings$CovarianceType
  eigen.rat <- settings$eigen.rat
  
  R <- CovarianceMatrixNum(X, X, theta, CovarianceFunction)
  R.alpha <- alpha * R + (1 - alpha) * diag(nrow = N)
  H <- PolyNum(X, poly.deg)
  
  R.alpha.inv <- SafeInvert(R.alpha, eigen.rat)
  N.mod <- R.alpha.inv$dim
  val <- R.alpha.inv$val
  R.alpha.inv <- R.alpha.inv$matrix
  
  H.R.alpha.H <- SafeInvert(t(H) %*% R.alpha.inv %*% H, eigen.rat)
  H.R.alpha.H <- H.R.alpha.H$matrix
  
  b <- H.R.alpha.H %*% t(H) %*% R.alpha.inv %*% f
  f.H.b <- f - H %*% b
  
  if (N.mod != N) {
      v <- ((t(f.H.b) %*% R.alpha.inv %*% f.H.b) / N.mod)[1]
      log.lik.val <- (N.mod * log(2 * pi) + N.mod * log(v) + sum(val) + N.mod) / 2
  } else { 
    v <- ((t(f.H.b) %*% R.alpha.inv %*% f.H.b) / N)[1]
    log.det <- determinant(R.alpha, logarithm = TRUE)
    log.lik.val <- (N * log(2 * pi) + N * log(v) +
                    (log.det$modulus[1] * log.det$sign) + N) / 2
  }
    
  return (list(v = v, val = log.lik.val))
}

LOOCV <- function (X, f, theta, alpha = 1, settings = krig.settings) {
  N <- nrow(X)
  L <- ncol(X)
  poly.deg <- settings$poly.deg
  CovarianceFunction <- settings$CovarianceType
  eigen.rat <- settings$eigen.rat
  
  R <- CovarianceMatrixNum(X, X, theta, CovarianceFunction)
  R.alpha <- alpha * R + (1 - alpha) * diag(nrow = N)
  H <- PolyNum(X, poly.deg)

  R.alpha.inv <- SafeInvert(R.alpha, eigen.rat)
  N.mod <- R.alpha.inv$dim
  val <- R.alpha.inv$val
  R.alpha.inv <- R.alpha.inv$matrix
  
  H.R.alpha.H <- SafeInvert(t(H) %*% R.alpha.inv %*% H, eigen.rat)
  H.R.alpha.H <- H.R.alpha.H$matrix

  b <- H.R.alpha.H %*% t(H) %*% R.alpha.inv %*% f
  f.H.b <- f - H %*% b

  if (N.mod != N) {
    v <- ((t(f.H.b) %*% R.alpha.inv %*% f.H.b) / N.mod)[1]
  } else { 
    v <- ((t(f.H.b) %*% R.alpha.inv %*% f.H.b) / N)[1]   
  }

  H.R.alpha <- t(H) %*% R.alpha.inv
  A <- R.alpha.inv %*% (diag(nrow = N) - H %*% H.R.alpha.H %*% H.R.alpha)
  diffs <- (t(f) %*% A) / diag(A)
  loocv.val <- sum(diffs ^ 2)
  
  return (list(v = v, val = loocv.val))
}

LogLikelihoodOpt <- function (X, f, par.vec, x.min, x.max,
                              settings = krig.settings) {
  theta.min <- settings$theta.min
  theta.max <- settings$theta.max
  alpha.min <- settings$alpha.min
  alpha.max <- settings$alpha.max
  
  theta <- ScaleTheta(par.vec[1:L], theta.min, theta.max,
                      x.min[1:L], x.max[1:L])
  alpha <- ScaleAlpha(par.vec[L + 1], alpha.min, alpha.max,
                      x.min[L + 1], x.max[L + 1])
  
  log.lik <- LogLikelihood(X, f, theta, alpha, settings)
  
  return (log.lik)
}

LOOCVOpt <- function (X, f, par.vec, x.min, x.max, settings = krig.settings) {
  theta.min <- settings$theta.min
  theta.max <- settings$theta.max
  alpha.min <- settings$alpha.min
  alpha.max <- settings$alpha.max
  
  theta <- ScaleTheta(par.vec[1:L], theta.min, theta.max,
                      x.min[1:L], x.max[1:L])
  alpha <- ScaleAlpha(par.vec[L + 1], alpha.min, alpha.max,
                      x.min[L + 1], x.max[L + 1])
  
  loocv <- LOOCV(X, f, theta, alpha, settings)
  
  return (loocv)
}

KrigingModel <- function (X, f, theta, v, alpha = 0.9999,
                          settings = krig.settings) {
  poly.deg <- settings$poly.deg
  CovarianceFunction <- settings$CovarianceType
  
  N <- nrow(X)
  L <- ncol(X)
  eigen.rat <- settings$eigen.rat
  sigma.sq <- alpha * v
  
  R <- CovarianceMatrixNum(X, X, theta, CovarianceFunction)
  R.alpha <- alpha * R + (1 - alpha) * diag(nrow = N)
  SIGMA <- v * R.alpha
  H <- PolyNum(X, poly.deg)
  
  SIGMA.inv <- SafeInvert(SIGMA, eigen.rat)
  N.mod <- SIGMA.inv$dim
  SIGMA.inv <- SIGMA.inv$matrix
  
  H.SIGMA.inv.H <- SafeInvert(t(H) %*% SIGMA.inv %*% H, eigen.rat)
  H.SIGMA.inv.H <- H.SIGMA.inv.H$matrix
  H.SIGMA.inv <- t(H) %*% SIGMA.inv
  
  A <- SIGMA.inv %*% (diag(nrow = N.mod) - H %*% H.SIGMA.inv.H %*% H.SIGMA.inv)
  B <- SIGMA.inv %*% H %*% H.SIGMA.inv.H
  # C <- H.SIGMA.inv.H %*% H.SIGMA.inv
  C <- t(B)
  D <- -H.SIGMA.inv.H
  M.inv <- cbind(rbind(A, C), rbind(B, D))
  
  r <- lag(CovarianceFunction(theta)$CovAlg, X)
  h <- Poly(L, poly.deg)
  
  EST <- t(rbind(f, matrix(0, poly.deg * L + 1))) %*% M.inv %*%
         V(sigma.sq * r, h)
  # VAR <- const.gAlg(calc(V(sigma.sq * CovarianceFunction(theta)$CovAlg),
  #                        matrix(0, 1, L)), L) -
  #        V(sigma.sq * r, h) %*% M.inv %*% V(sigma.sq * r, h)
  
  VAR <- const.gAlg(calc(V(sigma.sq * CovarianceFunction(theta)$CovAlg),
                          matrix(0, 1, L)), L) -
         V(sigma.sq * r, h) %*%
         (M.inv * (2 * upper.tri(M.inv) + diag(nrow = N.mod + poly.deg * L + 1))) %*%
         V(sigma.sq * r, h)
  
  VAR.t <- TruncGvec(VAR)$gVec
  EST.t <- TruncGvec(EST)$gVec
  
  ## TEST: gAlg object - trend
  beta <- H.SIGMA.inv.H %*% t(H) %*% SIGMA.inv %*% f
  MU <- h %*% beta 
  ##
  
  return (list(EST = EST.t, VAR = VAR.t, M.inv = M.inv, MU = MU))
}

RobustKrigingModel <- function (X, f, theta, v, alpha = 0.9999,
                                settings = krig.settings,
                                u.settings = uncertainty.settings) {
  
  ErrorType <- u.settings$ErrorType
  u.st.dev <- u.settings$st.dev
  u.delta <- u.settings$delta

  poly.deg <- settings$poly.deg
  CovarianceFunction <- settings$CovarianceType
  
  N <- nrow(X)
  L <- ncol(X)
  eigen.rat <- settings$eigen.rat
  sigma.sq <- alpha * v
  
  R <- CovarianceMatrixNum(X, X, theta, CovarianceFunction)
  R.alpha <- alpha * R + (1 - alpha) * diag(nrow = N)
  SIGMA <- v * R.alpha
  H <- PolyNum(X, poly.deg)
  
  SIGMA.inv <- SafeInvert(SIGMA, eigen.rat)
  N.mod <- SIGMA.inv$dim
  SIGMA.inv <- SIGMA.inv$matrix
  
  H.SIGMA.inv.H <- SafeInvert(t(H) %*% SIGMA.inv %*% H, eigen.rat)
  H.SIGMA.inv.H <- H.SIGMA.inv.H$matrix
  H.SIGMA.inv <- t(H) %*% SIGMA.inv
  
  A <- SIGMA.inv %*% (diag(nrow = N.mod) - H %*% H.SIGMA.inv.H %*% H.SIGMA.inv)
  B <- SIGMA.inv %*% H %*% H.SIGMA.inv.H
  # C <- H.SIGMA.inv.H %*% H.SIGMA.inv
  C <- t(B)
  D <- - H.SIGMA.inv.H
  M.inv <- cbind(rbind(A, C), rbind(B, D))
  
  if (ErrorType(1)$ErrName == "Gauss") {
    X.ERR <- ErrorType(u.st.dev)$ErrAlg
  } else if (ErrorType(1)$ErrName == "Beta-22") {
    X.ERR <- ErrorType(u.delta)$ErrAlg
  } else {
    stop("Wrong \"ErrName\" value.")
  }

  r.err <- lag(CovarianceFunction(theta)$CovAlg %% X.ERR, X)
  h.err <- Poly(L, poly.deg) %% X.ERR
  
  EST <- t(rbind(f, matrix(0, poly.deg * L + 1))) %*% M.inv %*%
         V(sigma.sq * r.err, h.err)
  # VAR <- const.gAlg(calc(V(X.ERR %% (sigma.sq * CovarianceFunction(theta)$CovAlg)
  #                          %% X.ERR), t(rep(0, L))), L) -
  #        V(sigma.sq * r.err, h.err) %*% M.inv %*% V(sigma.sq * r.err, h.err)
  
  VAR <- const.gAlg(calc(V(X.ERR %% (sigma.sq * CovarianceFunction(theta)$CovAlg)
                           %% X.ERR), t(rep(0, L))), L) -
         V(sigma.sq * r.err, h.err) %*%
         (M.inv * (2 * upper.tri(M.inv) + diag(nrow = N.mod + poly.deg * L + 1))) %*%
         V(sigma.sq * r.err, h.err)
  
  VAR.t <- TruncGvec(VAR)$gVec
  EST.t <- TruncGvec(EST)$gVec
  
  ## TEST: gAlg object - trend
  h <- Poly(L, poly.deg)
  beta <- H.SIGMA.inv.H %*% t(H) %*% SIGMA.inv %*% f
  MU <- h %*% beta 
  ##
  
  return (list(EST = EST.t, VAR = VAR.t, M.inv = M.inv, MU = MU))
}

FunctionProjection <- function (ObjectFunction, EST, X, f,
                                x.min, x.max, x.sec, dim.comb, n = 20) {
  x.1 <- seq(x.min[dim.comb[1]], x.max[dim.comb[1]], len = n)
  x.2 <- seq(x.min[dim.comb[2]], x.max[dim.comb[2]], len = n)
  
  M <- expand.grid(X1 = x.1, X2 = x.2)
  X.SEC <- data.frame(matrix(x.sec, n * n, L, byrow = TRUE))
  X.SEC[, dim.comb[1]] <- M$X1
  X.SEC[, dim.comb[2]] <- M$X2
  
  persp3d(matrix(x.1, n, n), t(matrix(x.2, n, n)),
          matrix(apply(X.SEC, 1, function (x) ObjectFunction(x)$val), n, n),
          col = "blue", alpha = 0.6, expand = 0.75,
          xlab = "", ylab = "", zlab = "")
  
  persp3d(matrix(x.1, n, n), t(matrix(x.2, n, n)),
          matrix(calc(EST, as.matrix(X.SEC)), n, n),
          col = "green", alpha = 0.6, expand = 0.75,
          xlab = "", ylab = "", zlab = "", add = TRUE)
  
  sel <- (X[,-dim.comb] == matrix(x.sec[-dim.comb], nrow(X), L - 2, byrow = TRUE))
  sel <- apply(sel, 1, function (x) all(x))
  
  points <- cbind(X[, c(dim.comb[1], dim.comb[2])], f)[sel, , drop = FALSE]
  
  points3d(points, size = 8)
}

ImprOptimum2 <- function (ImprFunction, x.min, x.max,
                          impr.settings, opt.settings) {
  
  opt.minimum <- opt.settings$find.minimum
  impr.minimum <- impr.settings$find.minimum
  impr.type <- impr.settings$imp.type
  gen.size <- impr.settings$gen.size
  pop.size <- impr.settings$pop.size
  gen.enhance <- impr.settings$gen.enhance
  
  L <- length(x.min)
  p.min <- rep(-10, L)
  p.max <- rep(10, L)
  
  if (impr.type == "Expected Improvement") {
    result.gen <- nsga2(function (P) {
      X <- ScaleDomainMat(P, x.min, x.max, p.min, p.max)
      f <- ImprFunction(X)
      f.1 <- (-1) ^ (opt.minimum - 1) * f$est
      f.2 <- -f$var
      f.3 <- (-1) ^ (impr.minimum - 1) * f$ei
      
      return (t(cbind(f.1, f.2, f.3)))
    }, L, 3,
    generations = gen.size, popsize = pop.size,
    lower.bounds = p.min, upper.bounds = p.max,
    vectorized = TRUE)
    
    sel <- which.min(result.gen$value[, 3])
    
    if (gen.enhance == TRUE) {
      result.grad <- optim(result.gen$par[sel, ],
                           function (p) {
                             X <- ScaleDomainMat(matrix(p, 1, L),
                                                 x.min, x.max, p.min, p.max)
                             (-1) ^ (impr.minimum - 1) * ImprFunction(X)$ei
                           },
                           method = "L-BFGS-B",
                           lower = p.min, upper = p.max,
                           control = list(trace = 0))
      
      if (result.grad$convergence == 0) {
        opt.coord <- ScaleDomainVec(result.grad$par, x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (impr.minimum - 1) * result.grad$val
      } else {
        warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                noBreaks. = TRUE)
        opt.coord <- ScaleDomainVec(result.gen$par[sel, ], x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 3]
      }
    } else if (gen.enhance == FALSE) {
      opt.coord <- ScaleDomainVec(result.gen$par[sel, ], x.min, x.max, p.min, p.max)
      opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 3]
    } else {
      stop("Wrong \"gen.enhance\" value.")
    }
    return (list(opt.coord = opt.coord, opt.value = opt.value))
  } else if (impr.type == "Relative Expected Improvement") {
    #
    #
  } else {
    stop("Wrong \"imp.type\" value.")
  }
}

ImprOptimum1 <- function (ImprFunction, x.min, x.max,
                          impr.settings, opt.settings) {
  
  opt.minimum <- opt.settings$find.minimum
  impr.minimum <- impr.settings$find.minimum
  impr.type <- impr.settings$imp.type
  gen.size <- impr.settings$gen.size
  pop.size <- impr.settings$pop.size
  gen.enhance <- impr.settings$gen.enhance

  x0.min <- impr.settings$x0.min
  x0.max <- impr.settings$x0.max
  x1.min <- impr.settings$x1.min
  x1.max <- impr.settings$x1.max
  
  L <- length(x.min)
  p.min <- rep(-10, L)
  p.max <- rep(10, L)
  
  if (impr.type == "Expected Improvement") {
    result.gen <- nsga2(function (P) {
                          X <- ScaleDomainMat(P, x.min, x.max, p.min, p.max)
                          f <- ImprFunction(X)
                          f.1 <- (-1) ^ (opt.minimum - 1) * f$est
                          f.2 <- (-1) ^ (impr.minimum - 1) * f$ei
                          
                          return (t(cbind(f.1, f.2)))
                        }, L, 2,
                        generations = gen.size, popsize = pop.size,
                        lower.bounds = p.min, upper.bounds = p.max,
                        vectorized = TRUE)
    
    sel <- which.min(result.gen$value[, 2])

    if (gen.enhance == TRUE) {
      result.grad <- optim(result.gen$par[sel, ],
                           function (p) {
                             X <- ScaleDomainMat(matrix(p, 1, L),
                                                 x.min, x.max, p.min, p.max)
                             (-1) ^ (impr.minimum - 1) * ImprFunction(X)$ei
                           },
                           method = "L-BFGS-B",
                           lower = p.min, upper = p.max,
                           control = list(trace = 0))
      
      if (result.grad$convergence == 0) {
        opt.coord <- ScaleDomainVec(result.grad$par, x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (impr.minimum - 1) * result.grad$val
      } else {
        warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                noBreaks. = TRUE)
        opt.coord <- ScaleDomainVec(result.gen$par[sel, ], x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 2]
      }
    } else if (gen.enhance == FALSE) {
      opt.coord <- ScaleDomainVec(result.gen$par[sel, ], x.min, x.max, p.min, p.max)
      opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 2]
    } else {
      stop("Wrong \"gen.enhance\" value.")
    }

    return (list(opt.coord = opt.coord, opt.value = opt.value))

  } else if (impr.type == "Relative Expected Improvement") {

    result.gen <- nsga2(function (P) {
                          X <- ScaleDomainMat(P,
                                              c(x1.min, x0.min), c(x1.max, x0.max),
                                              c(p.min, p.min), c(p.max, p.max))
                          f <- ImprFunction(X)
                          f.1 <- (-1) ^ (opt.minimum - 1) *
                                 (f$est - (-1) ^ (1 - opt.minimum) * 2 * sqrt(f$var))
                          f.2 <- (-1) ^ (impr.minimum - 1) * f$rei
                          
                          return (t(cbind(f.1, f.2)))
                        }, 2 * L, 2,
                        generations = gen.size, popsize = pop.size,
                        lower.bounds = c(p.min, p.min),
                        upper.bounds = c(p.max, p.max),
                        vectorized = TRUE)

    sel <- which.min(result.gen$value[, 2])

    if (gen.enhance == TRUE) {

      result.grad <- optim(result.gen$par[sel, ],
                           function (p) {
                             X <- ScaleDomainMat(matrix(p, 1, 2 * L),
                                                 c(x1.min, x0.min), c(x1.max, x0.max),
                                                 c(p.min, p.min), c(p.max, p.max))

                             (-1) ^ (impr.minimum - 1) * ImprFunction(X)$rei
                           },
                           method = "L-BFGS-B",
                           lower = c(p.min, p.min), upper = c(p.max, p.max),
                           control = list(trace = 0))
      
      if (result.grad$convergence == 0) {

        opt.coord <- ScaleDomainVec(result.grad$par,
                                    c(x1.min, x0.min), c(x1.max, x0.max),
                                    c(p.min, p.min), c(p.max, p.max))
        opt.value <- (-1) ^ (impr.minimum - 1) * result.grad$val

      } else {

        warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                noBreaks. = TRUE)
        opt.coord <- ScaleDomainVec(result.gen$par[sel, ],
                                    c(x1.min, x0.min), c(x1.max, x0.max),
                                    c(p.min, p.min), c(p.max, p.max))
        opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 2]

      }
    } else if (gen.enhance == FALSE) {

        opt.coord <- ScaleDomainVec(result.gen$par[sel, ],
                                    c(x1.min, x0.min), c(x1.max, x0.max),
                                    c(p.min, p.min), c(p.max, p.max))
        opt.value <- (-1) ^ (impr.minimum - 1) * result.gen$value[sel, 2]

    } else {

      stop("Wrong \"gen.enhance\" value.")

    }
    
    if ( (opt.value - .Machine$double.eps) < 0 ) {

      sel.2 <- which.min(result.gen$value[, 1])
      opt.coord <- ScaleDomainVec(result.gen$par[sel.2, ],
                                  c(x1.min, x0.min), c(x1.max, x0.max),
                                  c(p.min, p.min), c(p.max, p.max))
      opt.value = NA
    }

    return (list(opt.coord = opt.coord, opt.value = opt.value))

  } else {
    stop("Wrong \"imp.type\" value.")
  }
}

FunctionOptimum <- function (Function, x.min, x.max, settings) {
  search.type <- settings$opt.search.type
  gen.size <- settings$gen.size
  pop.size <- settings$pop.size
  gen.enhance <- settings$gen.enhance
  minimum <- settings$find.minimum
  
  dim <- length(x.min)
  p.min <- rep(-10, dim)
  p.max <- rep(10, dim)
  
  if (search.type == "genetic") {
    result.gen <- nsga2(function (P) {
                          (-1) ^ (minimum - 1) *
                          t(Function(ScaleDomainMat(P, x.min, x.max, p.min, p.max)))
                        }, dim, 1,
                        generations = gen.size, popsize = pop.size,      
                        lower.bounds = p.min, upper.bounds = p.max,
                        vectorized = TRUE)
    sel.gen <- which.min(result.gen$val)
      
    if (gen.enhance == TRUE) {
      result.grad <- optim(result.gen$par[sel.gen, ],
                           function (p) {
                             (-1) ^ (minimum - 1) *
                             Function(ScaleDomainMat(matrix(p, 1, dim),
                                                     x.min, x.max, p.min, p.max))
                           },
                           method = "L-BFGS-B",
                           lower = p.min, upper = p.max,
                           control = list(trace = 0))
      
      if (result.grad$convergence == 0) {
        opt.coord <- ScaleDomainVec(result.grad$par, x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (minimum - 1) * result.grad$val
      } else {
        warning(paste0("\"optim()\" failed to improve the \"nsga2()\" result.\n"),
                noBreaks. = TRUE)
        opt.coord <- ScaleDomainVec(result.gen$par[sel.gen, ], x.min, x.max, p.min, p.max)
        opt.value <- (-1) ^ (minimum - 1) * result.gen$val[sel.gen, ]
      }
    } else if (gen.enhance == FALSE) {
      opt.coord <- ScaleDomainVec(result.gen$par[sel.gen, ], x.min, x.max, p.min, p.max)
      opt.value <- (-1) ^ (minimum - 1) * result.gen$val[sel.gen, ]
    } else {
      stop("Wrong \"gen.enhance\" value.")
    }
  } else if (search.type == "gradient") {
    lhs.add.points <- settings$lhs.add.points
    lhs.point.number <- settings$lhs.point.number
    lhs.enhance <- settings$lhs.enhance
    
    if (is.null(lhs.add.points)) {
      opt.coord <- maximinLHS(lhs.point.number, dim)
      ##
      if (all(apply(opt.coord, 2, min) < rep(0, dim)) |
          all(apply(opt.coord, 2, max) > rep(1, dim))) {
        stop ("maximinLHS: ERROR!")
      }
      ##
    } else if (is.matrix(lhs.add.points)) {
      
      if (all(signif(apply(lhs.add.points, 2, min),
                     SignifPlaces(x.min) + 5) < x.min) ||
          all(signif(apply(lhs.add.points, 2, max),
                     SignifPlaces(x.max) + 5) > x.max)) {
        
        warning ("Additional points weren't used. Points out of range.")
        lhs.add.points <- matrix(x.min + (x.max - x.min) / 2, 1, dim)
      }
      
      lhs.add.points <- (lhs.add.points - x.min) * (1 / (x.max - x.min))
      opt.coord <- rbind(maximinLHS(lhs.point.number, dim), lhs.add.points)
      ##
      if (all(apply(opt.coord, 2, min) < rep(0, dim)) |
          all(apply(opt.coord, 2, max) > rep(1, dim))) {
        stop ("maximinLHS: add.points, ERROR!")
      }
      ##
    } else {
      stop ("Wrong \"lhs.add.points\" value.")
    }
    
    opt.coord <- t(p.min + t(opt.coord) * (p.max - p.min))
    ##
    if (all(apply(opt.coord, 2, min) < p.min) |
        all(apply(opt.coord, 2, max) > p.max)) {
      stop ("maximinLHS: scaled, ERROR!")
    }
    ##
    opt.value <- apply(opt.coord, 1, function (p) {
                                       x <- ScaleDomainMat(matrix(p, 1, dim), x.min, x.max,
                                                           p.min, p.max)
                                       
                                       ##
                                       if (all(apply(x, 2, min) < p.min) |
                                           all(apply(x, 2, max) > p.max)) {
                                         stop ("ScaleDomainMat: ERROR!")
                                       }
                                       ##
                                       
                                       return ((-1) ^ (minimum - 1) * Function(x))
                                     })
    if (lhs.enhance == TRUE) {
      lhs.enhance.number <- settings$lhs.enhance.number
      lhs.enhance.points <- ceiling(lhs.enhance.number * lhs.point.number)
      sel <- order(opt.value)[1:lhs.enhance.points]
      opt.coord.enhance <- matrix(0, lhs.enhance.points, dim)
      opt.value.enhance <- rep(0, lhs.enhance.points)
      
      for (i in 1:lhs.enhance.points) {
        result.grad <- optim(opt.coord[sel[i], ],
                             function (p) {
                               x <- ScaleDomainMat(matrix(p, 1, dim),
                                                   x.min, x.max, p.min, p.max)
                               
                               ##
                               if (all(apply(x, 2, min) < x.min) |
                                   all(apply(x, 2, max) > x.max)) {
                                 stop ("ScaleDomainMat, optim: ERROR!")
                               }
                               ##
                               
                               return ((-1) ^ (minimum - 1) * Function(x))
                             },
                             method = "L-BFGS-B",
                             lower = p.min,
                             upper = p.max)
        
        if (result.grad$convergence == 0) {
          
          ##
          if (all(result.grad$par[1:dim] < p.min) |
              all(result.grad$par[1:dim] > p.max)) {
            stop ("optim: ERROR!")
          }
          ##
          
          opt.coord.enhance[i,] <- ScaleDomainVec(result.grad$par[1:dim],
                                                  x.min, x.max, p.min, p.max)
          opt.value.enhance[i] <- (-1) ^ (minimum - 1) *
                                  Function(t(opt.coord.enhance[i,]))
          
        } else {
          warning(paste0("\"optim()\" failed to improve the \"LHS\" result.\n"),
                  noBreaks. = TRUE)
          opt.coord.enhance[i,] <- rep(NA, dim)
          opt.value.enhance[i] <- NA
        }
      }

      if (sum(opt.value.enhance, na.rm = TRUE) == 0) {
        sel <- which.min(opt.value)
        opt.coord <- ScaleDomainVec(opt.coord[sel,], x.min, x.max, p.min, p.max)
      } else {
        sel <- which.min(opt.value.enhance)
        opt.coord <- opt.coord.enhance[sel,]
      }
    } else if (lhs.enhance == FALSE) {
      sel <- which.min(opt.value)
      opt.coord <- ScaleDomainVec(opt.coord[sel,], x.min, x.max, p.min, p.max)
    } else {
      stop("Wrong \"lhs.enhance\" value.")
    }
    opt.value <- Function(matrix(opt.coord, 1, dim))
  } else {
    stop ("Wrong \"opt.search.type\" value.")
  }
  
  # results generated by FunctionOptimum are sometimes out of bounds
  # (reason unknown)
  if (all(opt.coord < x.min) |
      all(opt.coord > x.max)) {
    stop ("FunctionOptimum: function optimum coordinates out of range!")
  }
  
  return (list(opt.coord = opt.coord, opt.value = opt.value))
}

CrossValidation <- function (X, f, theta, v, alpha = 1,
                             settings = krig.settings) {
  N <- nrow(X)
  km.est = matrix(0, N, 1)
  km.var = matrix(0, N, 1)

  for (i in 1:N) {
    X.adj <- X[-i, , drop = FALSE]
    f.adj <- f[-i, , drop = FALSE]
    kriging.model <- KrigingModel(X.adj, f.adj, theta, v, alpha,
                                  settings = settings)
    km.est[i, 1] <- calc(kriging.model$EST, X[i, , drop = FALSE])
    km.var[i, 1] <- calc(kriging.model$VAR, X[i, , drop = FALSE])
  }
  
  km.var <- NumZeroSubstitute(km.var)
  # km.var <- ifelse(abs(km.var) < .Machine$double.eps, rep(0, N) , km.var)
  
  cv.err <- (f - km.est) / sqrt(km.var)
  
  return (list(cv.est = km.est, cv.var = km.var, cv.err = cv.err))
}

MeanValueAug <- function (R.EST, R.VAR, X, settings) {
  find.minimum <- settings$find.minimum
  N <- nrow(X)
  
  r.est.val <- calc(R.EST, X)
  r.var.val <- calc(R.VAR, X)
  
  r.var.val <- NumZeroSubstitute(r.var.val) 
  # r.var.val <- ifelse(abs(r.var.val) < .Machine$double.eps, rep(0, N), r.var.val)
  
  mean.val <- r.est.val + (-1) ^ (1 - find.minimum) * 2 * sqrt(r.var.val)
  
  return (list(val = mean.val, var = r.var.val))
}

ExpectedImprovementAug <- function (EST, VAR, X, opt.value,
                                    settings = optim.settings) {
  find.minimum <- settings$find.minimum
  N <- nrow(X)
  est.val <- calc(EST, X)
  var.val <- calc(VAR, X)
  var.val <- NumZeroSubstitute(var.val) 
  # var.val <- ifelse(abs(var.val) < .Machine$double.eps, rep(0, N) , var.val)
  
  F <- (-1) ^ (1 - find.minimum) * (matrix(opt.value, N, 1) - est.val)
  K <- sqrt(var.val)
  ei.val <- F * pnorm(F / K) + K * dnorm(F / K)
  
  return (list(ei = ei.val, est = est.val, var = var.val))
}

RelativeExpectedImprovementAug <- function (EST, VarianceFunction,
                                            X, X.1, X.0, M.inv, theta, v, alpha,
                                            opt.value,
                                            o.settings = optim.settings,
                                            k.settings = krig.settings,
                                            u.settings = uncertainty.settings) {
  find.minimum <- o.settings$find.minimum
  N <- nrow(X.1)
  est.val <- calc(EST, X.1)
  var.val <- VarianceFunction(X, X.1, X.0, M.inv, theta, v, alpha,
                              k.settings, u.settings)
  
  var.val <- NumZeroSubstitute(var.val)
  # var.val <- ifelse(abs(var.val) < .Machine$double.eps, rep(0, N), var.val)
  
  F <- (-1) ^ (1 - find.minimum) * (matrix(opt.value, N, 1) - est.val)
  K <- sqrt(var.val)
  rei.val <- F * pnorm(F / K) + K * dnorm(F / K)
  
  return (list(rei = rei.val, est = est.val, var = var.val))
}

InducedVariance <- function (X, X.1, X.0, M.inv, theta, v, alpha,
                             settings = krig.settings,
                             u.settings = uncertainty.settings) {
  
  ErrorType <- u.settings$ErrorType
  u.st.dev <- u.settings$st.dev
  u.delta <- u.settings$delta
  poly.deg <- settings$poly.deg
  CovarianceFunction <- settings$CovarianceType
  
  if (ErrorType(1)$ErrName == "Gauss") {
    X.ERR <- ErrorType(u.st.dev)$ErrAlg
  } else if (ErrorType(1)$ErrName == "Beta-22") {
    X.ERR <- ErrorType(u.delta)$ErrAlg
  } else {
    stop("Wrong \"ErrName\" value.")
  }
  
  sigma.sq <- alpha * v 
  tau.sq <- v * (1 - alpha)

  r <- CovarianceFunction(theta)$CovNum(X.0 - X.0)
  sigma <- sigma.sq * r + tau.sq
  
  k <- sigma.sq * CovarianceMatrixNum(X, X.0, theta, CovarianceFunction)
  l <- t(PolyNum(X.0, poly.deg))
  h <- t(PolyNum(X.1, poly.deg))
  
  r.err.alg <- lag(CovarianceFunction(theta)$CovAlg %% X.ERR, X.1)
  r.err <- sigma.sq * calc(r.err.alg, X)
  r2.err <- as.matrix(sigma.sq * diag(calc(r.err.alg, X.0)))
  
  val.1 <- (r2.err - colSums(rbind(k, l) * M.inv %*% rbind(r.err, h))) ^ 2
  val.2 <- sigma - colSums(rbind(k, l) * M.inv %*% rbind(k, l))
  val <- val.1 / val.2
  
  return (val)
}

SignifPlaces <- function (x) {
  n <- sapply(1:length(x), function (i) {
    match(TRUE, signif(x[i], 1:20) == x[i])
  })
  return (n)
}

ScaleTheta <- function (p, theta.min, theta.max, p.min, p.max) {
  x <- 1 / (1 + exp(-p))
  x.min <- 1 / (1 + exp(-p.min))
  x.max <- 1 / (1 + exp(-p.max))
  
  coeff <- (theta.max - theta.min) / (x.max - x.min)
  const <- (theta.min * x.max - theta.max * x.min) / (x.max - x.min)
  
  theta <- coeff * x + const
  return (theta)
}

ScaleAlpha <- function (p, alpha.min, alpha.max, p.min, p.max) {
  x <- 1 / (1 + exp(-p))
  x.min <- 1 / (1 + exp(-p.min))
  x.max <- 1 / (1 + exp(-p.max))
  
  coeff <- (alpha.max - alpha.min) / (x.max - x.min)
  const <- (alpha.min * x.max - alpha.max * x.min) / (x.max - x.min)
  
  alpha <- coeff * x + const
  return (alpha)
}

ScaleDomainMat <- function(p, x.min, x.max, p.min, p.max) {
  coeff <- (x.max - x.min) / (p.max - p.min)
  const <- (x.min * p.max - x.max * p.min) / (p.max - p.min)
  
  return (t(coeff * t(p) + const))
}

ScaleDomainVec <- function(p, x.min, x.max, p.min, p.max) {
  coeff <- (x.max - x.min) / (p.max - p.min)
  const <- (x.min * p.max - x.max * p.min) / (p.max - p.min)
  
  return (coeff * p + const)
}

REIoptim <- function (X.aug, x.min, x.max, iter.number, ObjectFunction,
                      krig.settings, optim.settings, rei.settings,
                      enable.cv = FALSE,
                      write.to.file = TRUE,
                      continue = FALSE) {

  if (continue == TRUE) {
    X.aug <- read.table(file = "results/km_data.csv", header = TRUE,
                        sep = ";", dec = ".")
    X.aug <- list(coord = as.matrix(X.aug[, 1:L]),
                  val = as.matrix(X.aug[, L + 1]),
                  state = as.matrix(X.aug[, L + 2]))
    
    Par.aug <- read.table(file = "results/km_par_data.csv", header = TRUE,
                          sep = ";", dec = ".")
    Par.aug <- list(Theta = as.matrix(Par.aug[, 1:L]),
                    Alpha = as.matrix(Par.aug[, L + 1]),
                    V = as.matrix(Par.aug[, L + 2]))
    
    Opt.aug <- read.table(file = "results/km_opt_data.csv", header = TRUE,
                          sep = ";", dec = ".")
    Obj.opt.aug <- list(coord = as.matrix(Opt.aug[, 1:L]),
                        value = as.matrix(Opt.aug[, L + 1]),
                        var = as.matrix(Opt.aug[, L + 2]))
    
    Rei.aug <- list(value = as.matrix(Opt.aug[, L + 3]),
                    x1.coord = as.matrix(Opt.aug[, (L + 4):(2 * L + 3)]),
                    x0.coord = as.matrix(Opt.aug[, (2 * L + 4):(3 * L + 3)]))
    
  } else if (continue == FALSE) {
    unlink("results", recursive = TRUE)
    Par.aug <- list(Theta = matrix(0, 0, L),
                    Alpha = matrix(0, 0, 1),
                    V = matrix(0, 0, 1))
    
    Obj.opt.aug <- list(coord = matrix(0, 0, L),
                        value = matrix(0, 0, 1),
                        var = matrix(0, 0, 1))
    
    Rei.aug <- list(value = matrix(0, 0, 1),
                    x1.coord = matrix(0, 0, L),
                    x0.coord = matrix(0, 0, L))
  } else {
    stop("Wrong \"continue\" value.")
  }

  X <- X.aug$coord
  f <- X.aug$val
  state <- X.aug$state
  
  Theta <- Par.aug$Theta
  Alpha <- Par.aug$Alpha
  V <- Par.aug$V
  
  Obj.opt.coord <- Obj.opt.aug$coord  
  Obj.opt.value <- Obj.opt.aug$value
  Obj.opt.var <- Obj.opt.aug$var
  
  Rei.value <- Rei.aug$value
  Rei.x1.coord <- Rei.aug$x1.coord
  Rei.x0.coord <- Rei.aug$x0.coord
  
  PrintSettings(X.aug, x.min, x.max, iter.number,
                krig.settings, optim.settings, rei.settings, uncertainty.settings,
                write.to.file = write.to.file)
  
  for (i in (nrow(X) + 1):(nrow(X) + iter.number)) {
    start.time.km <- proc.time()
    
    if (krig.settings$lhs.prev.points == TRUE) {
      if (nrow(Theta) != 0) {
        krig.settings$lhs.add.points <- tail(cbind(Theta, Alpha),
                                             1, addrownums = FALSE)
      } else {
        krig.settings$lhs.add.points <- NULL
      }
    } else if (krig.settings$lhs.prev.points == FALSE) {
      #
    } else {
      stop("Wrong \"lhs.prev.points\" value.")
    }
   
    X.conv <- X[state == "converged", , drop = FALSE]
    f.conv <- f[state == "converged", , drop = FALSE]
    
    kriging.par <- KrigingParOpt(X.conv, f.conv, krig.settings)
    theta <- kriging.par$theta.opt
    alpha <- kriging.par$alpha.opt
    v <- kriging.par$v.opt
    
    kriging.model.conv <- KrigingModel(X.conv, f.conv, theta, v, alpha)
    km.est.conv <- kriging.model.conv$EST
    
    f.temp <- matrix(NA, nrow(X), 1)

    for (i in 1:nrow(X)) {
      if (state[i, ] == "converged") {
        f.temp[i, ] <- f[i, ]
      } else if (state[i, ] == "unconverged") {
        val.temp <- calc(km.est.conv, X[i, , drop = FALSE])
        
        if (((optim.settings$find.minimum == TRUE) & (val.temp < min(f.conv))) |
            ((optim.settings$find.minimum == FALSE) & (val.temp > max(f.conv)))) {
          val.temp <- mean(f.conv)
          warning("Estimated values for unconverged simulations may cause
                  problem with finding optimum. Mean value used.")
        }
        
        f.temp[i, ] <- val.temp
        
      } else {
        stop ("Wrong \"state\" value.")
      }
    }
    
    Theta <- rbind(Theta, theta[1:L])
    Alpha <- rbind(Alpha, alpha[1])
    V <- rbind(V, v[1])
    
    robust.kriging.model <- RobustKrigingModel(X, f.temp, theta, v, alpha)
    r.km.est <- robust.kriging.model$EST
    r.km.var <- robust.kriging.model$VAR
    km.mu <- robust.kriging.model$MU
    M.inv <- robust.kriging.model$M.inv
    
    end.time.km <- proc.time()
    start.time.cv <- proc.time()
    
    if (enable.cv == TRUE) {
      cross.validation <- CrossValidation(X.conv, f.conv, theta, v, alpha)
      cross.validation$enabled = TRUE
    } else if (enable.cv == FALSE) {
      cross.validation <- list(enabled = FALSE)
    } else {
      stop("Wrong \"enable.cv\" value.")
    }
    
    end.time.cv <- proc.time()
    start.time.fo <- proc.time()
    
    if (optim.settings$find.minimum != TRUE &
        optim.settings$find.minimum != FALSE) {
      stop ("Wrong \"find.minimum\" value.")
    }
    
    mean.opt <- FunctionOptimum(function (x) {
                                  MeanValueAug(r.km.est, r.km.var, x, optim.settings)$val
                                },
                                x.min, x.max, optim.settings)
    obj.opt.value <- mean.opt$opt.value
    obj.opt.coord <- matrix(mean.opt$opt.coord, 1, L)
    
    end.time.fo <- proc.time()
    start.time.rei <- proc.time()
    
    if (rei.settings$opt.search.type == "gen-REI") {

      ImprFunction <- function (x) {
        X.1 <- x[, 1:L, drop = FALSE]
        X.0 <- x[, (L + 1):(2 * L), drop = FALSE]
        
        result <- RelativeExpectedImprovementAug(r.km.est, InducedVariance,
                                                 X, X.1, X.0,
                                                 M.inv, theta, v, alpha,
                                                 obj.opt.value)
        return (result)
      }
      
      rei.max <- ImprOptimum1(ImprFunction,
                              x.min, x.max, rei.settings, optim.settings)
      
      rei.max$opt.sampl <- rei.max$opt.coord[(L + 1):(2 * L)]
      rei.max$opt.coord <- rei.max$opt.coord[1:L]
      
      rei.x0.coord <- rei.max$opt.sampl
      rei.x1.coord <- rei.max$opt.coord
      rei.value <- rei.max$opt.value
      
    } else if (rei.settings$opt.search.type == "gradient" ||
               rei.settings$opt.search.type == "genetic") {
      
      ImprFunction <- function (x) {
        X.1 <- x[, 1:L, drop = FALSE]
        X.0 <- x[, (L + 1):(2 * L), drop = FALSE]
        
        result <- RelativeExpectedImprovementAug(r.km.est, InducedVariance,
                                                 X, X.1, X.0, M.inv, theta, v, alpha,
                                                 obj.opt.value)$rei
        
        return (result)
      }
      
      rei.max <- FunctionOptimum(ImprFunction,
                                 c(rei.settings$x1.min, rei.settings$x0.min),
                                 c(rei.settings$x1.max, rei.settings$x0.max),
                                 settings = rei.settings)

      rei.max$opt.sampl <- rei.max$opt.coord[(L + 1):(2 * L)]
      rei.max$opt.coord <- rei.max$opt.coord[1:L]
      
      rei.x0.coord <- rei.max$opt.sampl
      rei.x1.coord <- rei.max$opt.coord
      rei.value <- rei.max$opt.value
    } else {
      stop ("Wrong \"opt.search.type\" value.")
    }
    
    end.time.rei <- proc.time()
    
    duration <- list(km.duration = end.time.km - start.time.km,
                     cv.duration = end.time.cv - start.time.cv,
                     fo.duration = end.time.fo - start.time.fo,
                     rei.duration = end.time.rei - start.time.rei)
    
    PrintResults(X.aug, kriging.par, cross.validation, robust.kriging.model, rei.max,
                 obj.opt.value, obj.opt.coord, duration,
                 write.to.file = write.to.file)
    
    Obj.opt.coord <- rbind(Obj.opt.coord, obj.opt.coord[1, 1:L])
    Obj.opt.value <- rbind(Obj.opt.value, obj.opt.value[1])
    
    obj.opt.var <- calc(r.km.var, obj.opt.coord)
    
    obj.opt.var <- NumZeroSubstitute(obj.opt.var)
    # obj.opt.var <- ifelse(abs(obj.opt.var) < .Machine$double.eps, 0, obj.opt.var)
    
    Obj.opt.var <- rbind(Obj.opt.var, obj.opt.var)
    
    Rei.x0.coord <- rbind(Rei.x0.coord, rei.x0.coord[1:L])
    Rei.x1.coord <- rbind(Rei.x1.coord, rei.x1.coord[1:L])
    Rei.value <- rbind(Rei.value, rei.value[1])
    
    x.aug <- ObjectFunction(rei.x0.coord)
    X <- rbind(X, rei.x0.coord[1:L])
    f <- rbind(f, x.aug$val)
    state <- rbind(state, x.aug$state)
    X.aug$coord <- X
    X.aug$val <- f
    X.aug$state <- state
    
    if (write.to.file == TRUE) {
      LogIterations(X, f, state, Theta, Alpha, V,
                    Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                    Rei.value, Rei.x1.coord, Rei.x0.coord)
      
      PlotResults(X, f.temp, state, M.inv, x.min, x.max,
                  rei.settings$x0.min, rei.settings$x0.max,
                  rei.settings$x1.min, rei.settings$x1.max,
                  theta, v, alpha,
                  ObjectFunction,
                  r.km.est, r.km.var, km.mu, InducedVariance,
                  Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                  Rei.x1.coord,
                  plot.type, imp.type = "REI")
    }
  }
  
  return(list(X = X, f = f, state = state,
              Theta = Theta, Alpha = Alpha, V = V,
              Rei.value = Rei.value, Rei.x1.coord = Rei.x1.coord,
              Rei.x0.coord = Rei.x0.coord, 
              Obj.opt.coord = Obj.opt.coord, Obj.opt.value = Obj.opt.value,
              Obj.opt.var = Obj.opt.var))
}

EIoptim <- function (X.aug, x.min, x.max, iter.number, ObjectFunction,
                     krig.settings, optim.settings, ei.settings,
                     enable.cv = FALSE,
                     write.to.file = TRUE,
                     continue = FALSE) {
  
  if (continue == TRUE) {
    X.aug <- read.table(file = "results/km_data.csv", header = TRUE,
                        sep = ";", dec = ".")
    X.aug <- list(coord = as.matrix(X.aug[, 1:L]),
                  val = as.matrix(X.aug[, L + 1]),
                  state = as.matrix(X.aug[, L + 2]))
    Par.aug <- read.table(file = "results/km_par_data.csv", header = TRUE,
                          sep = ";", dec = ".")
    Par.aug <- list(Theta = as.matrix(Par.aug[, 1:L]),
                    Alpha = as.matrix(Par.aug[, L + 1]),
                    V = as.matrix(Par.aug[, L + 2]))
    Opt.aug <- read.table(file = "results/km_opt_data.csv", header = TRUE,
                          sep = ";", dec = ".")
    Obj.opt.aug <- list(coord = as.matrix(Opt.aug[, 1:L]),
                        value = as.matrix(Opt.aug[, L + 1]),
                        var = as.matrix(Opt.aug[, L + 2]))
    Ei.aug <- list(value = as.matrix(Opt.aug[, L + 3]),
                   coord = as.matrix(Opt.aug[, (L + 4):(2 * L + 3)]))
  } else if (continue == FALSE) {
    unlink("results", recursive = TRUE)
    Par.aug <- list(Theta = matrix(0, 0, L), Alpha = matrix(0, 0, 1),
                    V = matrix(0, 0, 1))
    
    Obj.opt.aug <- list(value = matrix(0, 0, 1), coord = matrix(0, 0, L),
                        var = matrix(0, 0, 1))
    
    Ei.aug <- list(coord = matrix(0, 0, L), value = matrix(0, 0, 1))
  } else {
    stop("Wrong \"continue\" value.")
  }
  
  X <- X.aug$coord
  f <- X.aug$val
  state <- X.aug$state
  
  Theta <- Par.aug$Theta
  Alpha <- Par.aug$Alpha
  V <- Par.aug$V
  
  Obj.opt.value <- Obj.opt.aug$value
  Obj.opt.coord <- Obj.opt.aug$coord
  Obj.opt.var <- Obj.opt.aug$var
  
  Ei.value <- Ei.aug$value
  Ei.coord <- Ei.aug$coord
  
  PrintSettings(X.aug, x.min, x.max, iter.number,
                krig.settings, optim.settings, ei.settings, uncertainty.settings,
                write.to.file = write.to.file)
  
  for (i in (nrow(X) + 1):(nrow(X) + iter.number)) {
    start.time.km <- proc.time()
    
    if (krig.settings$lhs.prev.points == TRUE) {
        if (nrow(Theta) != 0) {
          krig.settings$lhs.add.points <- tail(cbind(Theta, Alpha),
                                               1, addrownums = FALSE)
        } else {
          krig.settings$lhs.add.points <- NULL
        }
    } else if (krig.settings$lhs.prev.points == FALSE) {
      #
    } else {
      stop("Wrong \"lhs.prev.points\" value.")
    }
    
    X.conv <- X[state == "converged", , drop = FALSE]
    f.conv <- f[state == "converged", , drop = FALSE]
    
    kriging.par <- KrigingParOpt(X.conv, f.conv, krig.settings)
    theta <- kriging.par$theta.opt
    alpha <- kriging.par$alpha.opt
    v <- kriging.par$v.opt
    
    kriging.model.conv <- KrigingModel(X.conv, f.conv, theta, v, alpha)
    km.est.conv <- kriging.model.conv$EST
    f.temp <- matrix(NA, nrow(X), 1)
    
    for (i in 1:nrow(X)) {
      if (state[i, ] == "converged") {
        f.temp[i, ] <- f[i, ]
      } else if (state[i, ] == "unconverged") {
        val.temp <- calc(km.est.conv, X[i, , drop = FALSE])
        
        if (((optim.settings$find.minimum == TRUE) & (val.temp < min(f.conv))) |
            ((optim.settings$find.minimum == FALSE) & (val.temp > max(f.conv)))) {
          val.temp <- mean(f.conv)
          warning("Estimated values for unconverged simulations may cause
                  problem with finding optimum. Mean value used.")
          }

        f.temp[i, ] <- val.temp

      } else {
        stop ("Wrong \"state\" value.")
      }
    }
    
    Theta <- rbind(Theta, theta[1:L])
    Alpha <- rbind(Alpha, alpha[1])
    V <- rbind(V, v[1])
    
    kriging.model <- KrigingModel(X, f.temp, theta, v, alpha)
    km.est <- kriging.model$EST
    km.var <- kriging.model$VAR
    km.mu <- kriging.model$MU
    
    end.time.km <- proc.time()
    start.time.cv <- proc.time()

    if (enable.cv == TRUE) {
      cross.validation <- CrossValidation(X.conv, f.conv, theta, v, alpha)
      cross.validation$enabled = TRUE
    } else if (enable.cv == FALSE) {
      cross.validation <- list(enabled = FALSE)
    } else {
      stop("Wrong \"enable.cv\" value.")
    }
    
    end.time.cv <- proc.time()

    sel.min <- apply(t(t(X.conv) >= x.min), 1, all)
    sel.max <- apply(t(t(X.conv) <= x.max), 1, all)
    if (all(!(sel.min & sel.max))) {
      stop ("Wrong \"x.min\" and \"x.max\" values.")
    }
    X.conv.trunc <- X.conv[sel.min & sel.max, , drop = FALSE]
    f.conv.trunc <- f.conv[sel.min & sel.max, , drop = FALSE]
        
    if (optim.settings$find.minimum == TRUE) {
      obj.opt.value <- min(f.conv.trunc)
      obj.opt.coord <- X.conv.trunc[which.min(f.conv.trunc), , drop = FALSE]
    } else if (optim.settings$find.minimum == FALSE) {
      obj.opt.value <- max(f.conv.trunc)
      obj.opt.coord <- X.conv.trunc[which.max(f.conv.trunc), , drop = FALSE]
    } else {
      stop ("Wrong \"find.minimum\" value.")
    }

    start.time.ei <- proc.time()
    
    EiThreshold <- function (km.est, var.est,
                             ei.max, obj.opt.value, x.min, x.max,
                             i.settings = ei.settings,
                             o.settings = optim.settings) {
      
      ei.value <- ei.max$opt.value
      ei.coord <- ei.max$opt.coord
      
      if (abs(ei.value / obj.opt.value) < i.settings$threshold) {
        
        mv.opt <- FunctionOptimum(function (x) {
          var.est <- calc(km.var, x)
          var.est <- NumZeroSubstitute(var.est)
          
          return (calc(km.est, x) - (-1) ^ (1 - o.settings$find.minimum) *
                  2 * sqrt(var.est))
        },
        x.min, x.max, settings = o.settings)
        
        warning (paste0("\"ei.value / obj.opt.value\" smaller then ",
                        i.settings$threshold, "!\n"))
        return (list(opt.coord = mv.opt$opt.coord, opt.value = NA))
      }
      return (list(opt.coord = ei.coord, opt.value = ei.value))
    }
    
    if (ei.settings$opt.search.type == "gen-EI-1") {
      
      ImprFunction <- function (x) {
        return (ExpectedImprovementAug(km.est, km.var, x, obj.opt.value))
      }
      
      ei.max <- ImprOptimum1(ImprFunction,
                             x.min, x.max, ei.settings, optim.settings)
      
      ei.max <- EiThreshold(km.est, var.est, ei.max, obj.opt.value, x.min, x.max)
      ei.coord <- ei.max$opt.coord
      ei.value <- ei.max$opt.value
      
    } else if (ei.settings$opt.search.type == "gen-EI-2") {
      
      ImprFunction <- function (x) {
        return (ExpectedImprovementAug(km.est, km.var, x, obj.opt.value))
      }
      
      ei.max <- ImprOptimum2(ImprFunction,
                             x.min, x.max, ei.settings, optim.settings)
      
      ei.max <- EiThreshold(km.est, var.est, ei.max, obj.opt.value, x.min, x.max)
      ei.coord <- ei.max$opt.coord
      ei.value <- ei.max$opt.value
      
    } else if (ei.settings$opt.search.type == "gradient" |
               ei.settings$opt.search.type == "genetic") {
      ei.max <- FunctionOptimum(function (x) {
                                  ExpectedImprovementAug(km.est, km.var,
                                  x, obj.opt.value)$ei
                                },
                x.min, x.max, settings = ei.settings)
      
      ei.max <- EiThreshold(km.est, var.est, ei.max, obj.opt.value, x.min, x.max)
      ei.coord <- ei.max$opt.coord
      ei.value <- ei.max$opt.value
      
    } else {
      stop ("Wrong \"opt.search.type\" value.")
    }
    
    end.time.ei <- proc.time()
    
    duration <- list(km.duration = end.time.km - start.time.km,
                     cv.duration = end.time.cv - start.time.cv,
                     ei.duration = end.time.ei - start.time.ei)

    PrintResults(X.aug, kriging.par, cross.validation, kriging.model, ei.max,
                 obj.opt.value, obj.opt.coord, duration,
                 write.to.file = write.to.file)
    
    Obj.opt.coord <- rbind(Obj.opt.coord, obj.opt.coord[1, 1:L])
    Obj.opt.value <- rbind(Obj.opt.value, obj.opt.value[1])
    
    obj.opt.var <- calc(km.var, obj.opt.coord)
    
    obj.opt.var <- NumZeroSubstitute(obj.opt.var)
    #obj.opt.var <- ifelse(abs(obj.opt.var) < .Machine$double.eps, 0, obj.opt.var)
    
    Obj.opt.var <- rbind(Obj.opt.var, obj.opt.var)
    
    Ei.coord <- rbind(Ei.coord, ei.coord[1:L])
    Ei.value <- rbind(Ei.value, ei.value[1])
    
    x.aug <- ObjectFunction(ei.coord)
    X <- rbind(X, ei.coord[1:L])
    f <- rbind(f, x.aug$val)
    state <- rbind(state, x.aug$state)
    X.aug$coord <- X
    X.aug$val <- f
    X.aug$state <- state
    
    if (write.to.file == TRUE) {
      LogIterations(X, f, state, Theta, Alpha, V,
                    Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                    Ei.value, Ei.coord)

      PlotResults(X, f.temp, state, M.inv, x.min, x.max,
                  NULL, NULL, NULL, NULL, theta, v, alpha,
                  ObjectFunction,
                  km.est, km.var, km.mu, InducedVariance,
                  Obj.opt.coord, Obj.opt.value, Obj.opt.var, NULL,
                  plot.type, imp.type = "EI")
    }
  }
  
  return(list(X = X, f = f, state = state,
              Theta = Theta, Alpha = Alpha, V = V,
              Ei.coord = Ei.coord, Ei.value = Ei.value,
              Obj.opt.coord = Obj.opt.coord, Obj.opt.value = Obj.opt.value,
              Obj.opt.var = Obj.opt.var))
}

PrintSettings <- function (X.aug, x.min, x.max, iter.number,
                           krig.settings, optim.settings, i.settings,
                           uncertainty.settings,
                           write.to.file = TRUE) {
  
  messages <- function (X.aug, x.min, x.max, iter.number,
                        krig.settings, optim.settings, i.settings,
                        uncertainty.settings) {
    X <- X.aug$coord
    state <- X.aug$state
    N.conv <- sum(state == "converged")
    N <- nrow(X)
    L <- ncol(X)
    
    k.covariance.name <- krig.settings$CovarianceType(1)$CovName
    k.poly.deg <- krig.settings$poly.deg
    k.theta.min <- krig.settings$theta.min
    k.theta.max <- krig.settings$theta.max
    k.alpha.min <- krig.settings$alpha.min
    k.alpha.max <- krig.settings$alpha.max

    i.imp.type <- i.settings$imp.type
    i.x0.min <- i.settings$x0.min
    i.x0.max <- i.settings$x0.max
    i.x1.min <- i.settings$x1.min
    i.x1.max <- i.settings$x1.max
    
    u.error.name <- uncertainty.settings$ErrorType(1)$ErrName
    u.error.st.dev <- uncertainty.settings$st.dev
    u.error.delta <- uncertainty.settings$delta

    print.vec <- function (vector, dig, col, tab) {
      len <- length(vector)
      for (i in 1:ceiling(len / col)) {
        cat(do.call(paste, as.list(c(rep(" ", tab), sep=""))),
            split(format(vector, digits = dig),
                  ceiling(seq_along(vector) / col))[[i]], "\n")
      }
    }
    
    print.opt.set <- function (opt.settings) {
      find.minimum <- opt.settings$find.minimum
      opt.search.type <- opt.settings$opt.search.type
      gen.enhance <- opt.settings$gen.enhance
      gen.size <- opt.settings$gen.size
      pop.size <- opt.settings$pop.size
      lhs.point.number <- opt.settings$lhs.point.number
      lhs.add.points <- opt.settings$lhs.add.points
      lhs.prev.points <- opt.settings$lhs.prev.points
      lhs.enhance <- opt.settings$lhs.enhance
      lhs.enhance.number <- opt.settings$lhs.enhance.number
    
      if (opt.search.type == "genetic") {
        cat("  Method: genetic algorithm (nsga2)\n")
        if (!is.null(find.minimum)) {
          cat(paste0("  Minimum: ", find.minimum, "\n"))
        }
        cat(paste0("  Generation size: ", gen.size, "\n"))
        cat(paste0("  Population size: ", pop.size, "\n"))
        cat(paste0("  Improving result with gradient algorithm: ",
                   gen.enhance, "\n"))
      } else if (opt.search.type == "gradient") {
        cat("  Method: sampling / gradient algorithm (LHS / L-BFGS-B)\n")
        if (!is.null(find.minimum)) {
          cat(paste0("  Minimum: ", find.minimum, "\n"))
        }
        cat(paste0("  Number of LHS generated points: ", lhs.point.number, "\n"))
        
        if (is.null(lhs.prev.points)) {
          
          if (is.null(lhs.add.points)) {
            cat("  Additional points:  FALSE\n")
          } else if (is.matrix(lhs.add.points)) {
            cat("  Additional points:\n")
            for (j in 1:nrow(lhs.add.points)) {
              print.vec(lhs.add.points[j, ], 4, 8, 2)
            }
          } else {
            stop ("Wrong \"lhs.add.points\" value.")
          }
          
        } else {
         
          if (lhs.prev.points == TRUE) {
            cat("  Additional points: Previous kriging parameters used.\n")
          } else if (lhs.prev.points == FALSE) {
            
            if (is.null(lhs.add.points)) {
              cat("  Additional points:  FALSE\n")
            } else if (is.matrix(lhs.add.points)) {
              cat("  Additional points:\n")
              for (j in 1:nrow(lhs.add.points)) {
                print.vec(lhs.add.points[j, ], 4, 8, 2)
              }
            } else {
              stop ("Wrong \"lhs.add.points\" value.")
            }
          }
           
        }

        cat("  Improving LHS results: ", lhs.enhance,"\n")
        if (lhs.enhance == TRUE) {
          cat("  Percent of best results used: ",
              lhs.enhance.number * 100, "[%]\n")
        }
      } else if (opt.search.type == "gen-EI-1" &&
                 i.imp.type == "Expected Improvement") {
        
        cat("  MOO: Expected Improvement and Objective Function\n")
        cat("  Method: genetic algorithm (nsga2)\n")
        if (is.null(find.minimum)) {
          cat(paste0("  Minimum: ", find.minimum, "\n"))
        }
        cat(paste0("  Generation size: ", gen.size, "\n"))
        cat(paste0("  Population size: ", pop.size, "\n"))
        cat(paste0("  Improving result with gradient algorithm: ",
                   gen.enhance, "\n"))
        
      } else if (opt.search.type == "gen-EI-2" &&
                 i.imp.type == "Expected Improvement") {
        
        cat("  MOO: Expected Improvement, Objective Function and Function Variance\n")
        cat("  Method: genetic algorithm (nsga2)\n")
        if (is.null(find.minimum)) {
          cat(paste0("  Minimum: ", find.minimum, "\n"))
        }
        cat(paste0("  Generation size: ", gen.size, "\n"))
        cat(paste0("  Population size: ", pop.size, "\n"))
        cat(paste0("  Improving result with gradient algorithm: ",
                   gen.enhance, "\n"))
        
      } else if (opt.search.type == "disabled" &&
                 !is.null(opt.settings$CovarianceType)) {
        cat("  Method: search disabled.\n")
        cat("  Arbitrary kriging parameters: theta.mean, alpha.max\n")
      } else if (opt.search.type == "gen-REI" &
                 i.imp.type == "Relative Expected Improvement") {
        
        cat("  MOO: Relative Expected Improvement and Lower/Upper Interval of\n")
        cat("  Mean Response\n")
        cat("  Method: genetic algorithm (nsga2)\n")
        if (is.null(find.minimum)) {
          cat(paste0("  Minimum: ", find.minimum, "\n"))
        }
        cat(paste0("  Generation size: ", gen.size, "\n"))
        cat(paste0("  Population size: ", pop.size, "\n"))
        cat(paste0("  Improving result with gradient algorithm: ",
                   gen.enhance, "\n"))
        
      } else {
        stop ("Wrong \"opt.search.type\" value.")
      }
    }

    cat(paste0("\nRun date:\n", Sys.time(), "\n"))
    
    cat("\n* Kriging model settings:\n")
    cat(paste0("  Covariance function: ", k.covariance.name,"\n"))
    cat(paste0("  Polynomial degree: ", k.poly.deg, "\n"))
    cat("  Kriging parameters limits:\n")
    cat("  theta.min =")
    print.vec(k.theta.min, 4, 8, 0)
    cat("  theta.max =")
    print.vec(k.theta.max, 4, 8, 0)
    cat(paste0("  alpha.min = ", format(k.alpha.min, digits = 6), "\n"))
    cat(paste0("  alpha.max = ", format(k.alpha.max, digits = 6), "\n"))
    cat("\n* Kriging parameters optimization processes:\n")
    print.opt.set(krig.settings)
    
    if (i.imp.type == "Relative Expected Improvement") {
      cat("\n* Error properties:\n")
      cat(paste0("  Error function: ", u.error.name,"\n"))
      
      if (u.error.name == "Gauss") {
        cat("  st.dev =")
        print.vec(u.error.st.dev, 4, 8, 0)
      } else if (u.error.name == "Beta-22") {
        cat("  delta =")
        print.vec(u.error.delta, 4, 8, 0)
      } else {
        stop ("Wrong \"u.error.name\" value.")
      }
    }

    if (i.imp.type == "Expected Improvement") {
      cat("\nFunction optimization processes:\n")
      cat("  Iteration limit: ", iter.number, "\n")
      cat("  Number of initial points: ", N, "\n")
      cat("  Number of converged evaluations: ", N.conv, "\n")
      cat("  Search domain:\n")
      cat("  x.min =")
      print.vec(x.min, 4, 8, 0)
      cat("  x.max =")
      print.vec(x.max, 4, 8, 0)
      cat("\n* Function optimum:\n")
      cat(paste0("  minimum: ", optim.settings$find.minimum, "\n"))
      cat("  method: selection of optimal measured point.\n")
      cat("\n* Expected Improvement maximum:\n")
      print.opt.set(ei.settings)
      
    } else if (i.imp.type == "Relative Expected Improvement") {
      
      cat("\nFunction optimization processes:\n")
      cat("  Iteration limit: ", iter.number, "\n")
      cat("  Number of initial points: ", N, "\n")
      cat("  Number of converged evaluations: ", N.conv, "\n")
      cat("  Search domain:\n")
      cat("  x.min =")
      print.vec(x.min, 4, 8, 0)
      cat("  x.max =")
      print.vec(x.max, 4, 8, 0)
      cat("  x0.min =")
      print.vec(i.x0.min, 4, 8, 0)
      cat("  x0.max =")
      print.vec(i.x0.max, 4, 8, 0)
      cat("  x1.min =")
      print.vec(i.x1.min, 4, 8, 0)
      cat("  x1.max =")
      print.vec(i.x1.max, 4, 8, 0)
      cat("\n* Mean response optimum:\n")
      print.opt.set(optim.settings)
      cat("\n* Relative Expected Improvement maximum:\n")
      print.opt.set(rei.settings)
    } else {
      stop ("Wrong \"imp.type\" value.")
    }
  }
  
  if (write.to.file == TRUE) {
    messages(X.aug, x.min, x.max, iter.number,
             krig.settings, optim.settings, i.settings, uncertainty.settings)
    dir.create("results", showWarnings = FALSE, recursive = FALSE, mode = "755")
    sink(file = "results/settings.log", append = TRUE, type = "output")
    messages(X.aug, x.min, x.max, iter.number,
             krig.settings, optim.settings, i.settings, uncertainty.settings)
    sink()
  } else if (write.to.file == FALSE) {
    messages(X.aug, x.min, x.max, iter.number,
             krig.settings, optim.settings, i.settings, uncertainty.settings)
  } else {
    stop("Wrong \"write.to.file\" value.")
  }
}

PrintResults <- function (X.aug, kriging.par, cross.validation,
                          kriging.model, improvement, f.opt, x.f.opt, duration,
                          write.to.file = TRUE) {
  
  messages <- function (X.aug, kriging.par, cross.validation,
                        kriging.model, improvement, f.opt, x.f.opt, duration) {
    theta <- kriging.par$theta.opt
    alpha <- kriging.par$alpha.opt
    v <- kriging.par$v.opt
    par.fun.type <- kriging.par$type
    par.fun.val <- kriging.par$par.fun.opt
    cv.err <- cross.validation$cv.err
    km.est <- kriging.model$EST
    km.var <- kriging.model$VAR
    i.coord <- improvement$opt.coord
    i.sampl <- improvement$opt.sampl
    i.value <- improvement$opt.value
    
    X <- X.aug$coord
    f <- X.aug$val
    N <- nrow(X)
    state <- X.aug$state
    X.conv <- X[state == "converged", , drop = FALSE]
    f.conv <- f[state == "converged", , drop = FALSE]
    N.conv <- nrow(X.conv)
    ev <- N + 1
    fun.range <- max(f.conv) - min(f.conv)
    fun.range <- as.vector((calc(km.est, X.conv) - f.conv) / fun.range * 100)
    
    if (optim.settings$find.minimum == TRUE) {
      opt.type = "minimum"
    } else if (optim.settings$find.minimum == FALSE) {
      opt.type = "maximum"
    }
    
    print.vec <- function (vector, dig, col, tab) {
      len <- length(vector)
      for (i in 1:ceiling(len / col)) {
        cat(do.call(paste, as.list(c(rep(" ", tab), sep=""))),
            split(format(vector, digits = dig),
                  ceiling(seq_along(vector) / col))[[i]], "\n")
      }
    }
    
    print.time <- function (time, dig) {
      text <- names(time)
      dim <- length(text)
      
      for (i in 1:dim) {
        if (time[i] / 60 >= 1) {
          cat(paste0(" ", text[i], " = "),
              format(time[i] / 60, dig, trim=T), "[min]\n")
        } else {
          cat(paste0(" ", text[i], " = "),
              format(time[i], dig, trim=T), "[s]\n") 
        }
      }
    }
    
    cat(paste0("\nFunction evaluation [", ev, "]:\n\n"))
    
    cat("Objective function convergence summary:\n")
    cat(paste0(" Points number = ", N, "\n"))
    cat(paste0(" Converged function evaluations = ", N.conv, "\n\n"))
    
    cat(paste0("Kriging parameters optimization:\n",
               " Function of parameters type: ",
               par.fun.type, "\n"))
    
    if (par.fun.type == "LogLikelihood") {
      cat(paste0(" Loglikelihood value: ",
                 format(-par.fun.val, digits = 4), "\n"))      
    } else if (par.fun.type == "LOOCV") {
      cat(paste0(" LOOCV value: ",
                 format(par.fun.val, digits = 4), "\n"))
    }

    cat(paste0(" Parameters: \n",
               " theta = "))
    print.vec(theta, 4, 6, 0)
    cat(paste0(" alpha = ", format(alpha, digits = 4), "\n",
               " v = ", format(v, digits = 4), "\n\n"))
    
    cat("Cross - Validation:\n")
    if (cross.validation$enabled == TRUE) {
      cat(" Standardized residuals:\n")
      print.vec(cv.err, 4, 6, 0) 
    } else if (cross.validation$enabled == FALSE) {
      cat(" Check disabled.\n")
    } else {
      stop("Wrong \"enabled\" value.")
    }
    cat("\n")
    
    if (is.null(i.sampl)) {
      cat(paste("Basic kriging model check:\n",
                "Differences between obj. function and its estimates [%]: \n"))
      print.vec(fun.range, 3, 6, 0)
      cat("\n")
      
      cat(paste0("Objective function optimum: \n",
                 " Optimum type: ", opt.type, "\n",
                 " f.opt = ", format(f.opt, digits = 4), "\n",
                 " x.f.opt = "))
      print.vec(x.f.opt, 4, 6, 0)
      cat("\n")
      
      cat(paste0("Expected Improvement maximum: \n",
                 " ei.max = ", format(i.value, digits = 4), "\n",
                 " ei.max/f.opt = ",
                 format(100 * i.value / f.opt, digits = 3), "[%]\n",
                 " x.ei = "))
      print.vec(i.coord, 4, 6, 0)
      
      cat("\nKriging process execution time:\n")
      cat(" creating kriging model:\n")
      print.time(duration$km.duration[3], 3)
      
      if (cross.validation$enabled == TRUE) {
        cat(" cross-validation check:\n")
        print.time(duration$cv.duration[3], 3)
      }
      
      cat(" searching for ei maximum:\n")
      print.time(duration$ei.duration[3], 3)

    } else if (is.vector(i.sampl)) {
      cat(paste0("Mean Response optimum: \n",
                 " Optimum type: ", opt.type, "\n",
                 " mr.opt = ", format(f.opt, digits = 4), "\n",
                 " x.mr.opt = "))
      print.vec(x.f.opt, 4, 6, 0)
      cat("\n")
      
      cat("Relative Expected Improvement maximum: \n")
      if (is.na(i.value)) {
        cat("Positive value or REI not found. Optimum of lower/upper \n")
        cat("confidence interval of Mean Function used. \n")
        cat(paste0(" x1.imf = "))
        print.vec(i.coord, 4, 6, 0)
        cat(paste0(" x0.imf = "))
        print.vec(i.sampl, 4, 6, 0)
      } else {
        cat(paste0(" rei.max = ", format(i.value, digits = 4), "\n",
                   " rei.max/f.opt = ",
                   format(100 * i.value / f.opt, digits = 3), "[%]\n",
                   " x1.rei = "))
        print.vec(i.coord, 4, 6, 0)
        cat(paste0(" x0.rei = "))
        print.vec(i.sampl, 4, 6, 0) 
      }
      
      cat("\nKriging process execution time:\n")
      cat(" creating kriging model:\n")
      print.time(duration$km.duration[3], 3)
      
      if (cross.validation$enabled == TRUE) {
        cat(" cross-validation check:\n")
        print.time(duration$cv.duration[3], 3)
      }
      
      cat(" searching for mean response optimum:\n")
      print.time(duration$fo.duration[3], 3)
      
      cat(" searching for rei maximum:\n")
      print.time(duration$rei.duration[3], 3)
      
    }
  }
  
  if (write.to.file == TRUE) {
    messages(X.aug, kriging.par, cross.validation,
             kriging.model, improvement, f.opt, x.f.opt, duration)
    dir.create("results", showWarnings = FALSE, recursive = FALSE, mode = "755")
    sink(file = "results/kriging.log", append = TRUE, type = "output")
    messages(X.aug, kriging.par, cross.validation,
             kriging.model, improvement, f.opt, x.f.opt, duration)
    sink()
  } else if (write.to.file == FALSE) {
    messages(X.aug, kriging.par, cross.validation,
             kriging.model, improvement, f.opt, x.f.opt, duration)
  } else {
    stop("Wrong \"write.to.file\" value.")
  }
}

LogIterations <- function (X, f, State, Theta, Alpha, V,
                           Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                           I.value, I.coord, I.sampl.coord = NULL) {
  
  write.table(cbind(X, f, State), file = "results/km_data.csv",
              sep = ";", dec = ".",
              col.names = c(sapply(1:L, function (i) paste0("x_", i)),
                            "f", "state"),
              row.names = FALSE, quote = FALSE)
  
  write.table(cbind(Theta, Alpha, V), file = "results/km_par_data.csv",
              sep = ";", dec = ".",
              col.names = c(sapply(1:L, function (i) paste0("theta_", i)),
                            "alpha", "v"),
              row.names = FALSE, quote = FALSE)
  
  if (is.null(I.sampl.coord)) {
    write.table(cbind(Obj.opt.coord, Obj.opt.value, Obj.opt.var, I.value,
                      I.coord),
                file = "results/km_opt_data.csv", sep = ";", dec = ".",
                col.names = c(sapply(1:L, function (i) paste0("x_opt_", i)),
                              "f_opt", "f_opt_var",
                              "ei", sapply(1:L, function (i) paste0("x_ei_", i))),
                row.names = FALSE, quote = FALSE)
  } else if (is.vector(I.sampl.coord) || is.matrix(I.sampl.coord)) {
    write.table(cbind(Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                      I.value, I.coord, I.sampl.coord),
                file = "results/km_opt_data.csv", sep = ";", dec = ".",
                col.names = c(sapply(1:L, function (i) paste0("x_opt_", i)),
                              "mr_opt", "mr_opt_var",
                              "rei", sapply(1:L, function (i) paste0("x_rei_", i)),
                              sapply(1:L, function (i) paste0("x0_rei_", i))),
                row.names = FALSE, quote = FALSE)
  } else {
    stop("Wrong \"I.sampl.coord\" value.")
  }
}

PlotResults <- function (X, f.temp, state, M.inv, x.min, x.max,
                         x0.min, x0.max, x1.min, x1.max, theta, v, alpha,
                         ObjectFunction,
                         km.est, km.var, km.mu, InducedVariance,
                         Obj.opt.coord, Obj.opt.value, Obj.opt.var,
                         Rei.x1.coord,
                         plot.type, imp.type,
                         k.settings = krig.settings,
                         u.settings = uncertainty.settings) {
  L <- ncol(X)
  N <- nrow(X) - 1
  iter.number <- nrow(Obj.opt.value)
  
  custom.grid <- function (lty = 2, ...) {
    double <- function (x) {
      x <- sort(x)
      d <- (x[2] - x[1]) / 2
      seq(x[1] - d, tail(x, 1) + d, d)
    }
    x <- double(axTicks(1));
    abline(v = x, lty = lty, ...)
    x <- double(axTicks(2));
    abline(h = x, lty = lty, ...)
  }
  
  open.dev <- function (plot.type, file.name, N) {
    if (plot.type == "png") {
      png(file = paste0("results/plots/", file.name, "_", N, ".png"),
          width = 640, height = 640)
    } else if (plot.type == "pdf") {
      pdf(file = paste0("results/plots/", file.name, "_", N, ".pdf"),
          onefile = TRUE, useDingbats = FALSE)
    } else {
      stop("Wrong \"plot.type\" value.")
    }
  }
  
  dir.create("results/plots", showWarnings = FALSE,
             recursive = FALSE, mode = "755")
  
  if (imp.type == "EI") {
    
    if (L == 1) {
      N.p <- 200
      X.obj.fun <- matrix(seq(x.min, x.max, len = N.p), N.p, 1)
      f.obj.fun <- ObjectFunction(X.obj.fun)$val
      
      f.est <- calc(km.est, X.obj.fun)
      var.est <- calc(km.var, X.obj.fun)
      
      var.est <- NumZeroSubstitute(var.est)
      # var.est <- ifelse(abs(var.est) < .Machine$double.eps, rep(0, N.p) , var.est)
      
      mu.est <- calc(km.mu, X.obj.fun)
      
      ei <- ExpectedImprovementAug(km.est, km.var, X.obj.fun,
                                tail(Obj.opt.value, 1))$ei

      y.min <- min(rbind(f.obj.fun, f.est, f.est - 2 * sqrt(var.est)))
      y.max <- max(rbind(f.obj.fun, f.est, f.est + 2 * sqrt(var.est)))
      ei.max <- max(ei)
      ei <- ei / ei.max * 0.05 * (y.max - y.min)
      
      colors <- rep("black", N)
      for (i in 1:N) {
        if (state[i,] == "unconverged")
          colors[i] <- "red"
      }
      
      open.dev(plot.type, "kriging", N)
        par(las = 1)
        matplot(X.obj.fun,
                cbind(f.obj.fun, f.est, f.est + 2 * sqrt(var.est), ei, mu.est),
                type = "l", col = 1:5,
                lty = c(1, 1, 2, 1, 1), lwd = c(2, 1.5, 1.5, 1.5, 1.5),
                xlab = "", ylab = "", xaxt = "n", yaxt = "n",
                ylim = c(y.min, y.max))
        matplot(X.obj.fun, f.est - 2 * sqrt(var.est),
                type = "l", col = 3,
                lty = 2, lwd = 1.5,
                add = TRUE)
        polygon(c(X.obj.fun, rev(X.obj.fun)),
                c(f.est + 2 * sqrt(var.est), rev(f.est - 2 * sqrt(var.est))),
                col = rgb(0, 1, 0, 0.1), border = FALSE)
        abline(v = X[N + 1, ], col = 4)
        custom.grid(lwd = 1, col = "gray")
        legend("bottomleft",
               c(expression('F'(x)), expression(widehat('F')(x)),
                 expression(widehat('F')(x)%+-%2*sqrt(widehat(sigma)(x))),
                 expression(EI(x)), expression(mu(x))),
               col = 1:5, lty = c(1, 1, 2, 1, 1),
               lwd = c(2, 1.5, 1.5, 1.5, 1.5), bty = "n")
        points(X[1:N, ], f.temp[1:N, ], pch = 19, col = colors)
        axis(1, at = NULL, tcl = 0.25, lwd = 2, las = 1)
        axis(2, at = NULL, tcl = 0.25, lwd = 2, las = 2)
        par(las = 0)
      dev.off()
      
    } else {
     
      colors <- c(rep("black", N - iter.number + 1),
                  rep("blue", iter.number - 1),
                  "green", "red")
      
      for (i in 1:N) {
        if (state[i,] == "unconverged")
          colors[i] <- "orange"
      }
      
      open.dev(plot.type, "kriging", N)
      par(las = 1)
      pairs(rbind(X, tail(Obj.opt.coord, 1)),
            labels = sapply(1:L, function(i) paste0("x_", i)), 
            main = "Points locations",
            col = colors,
            pch = 19, oma = c(5, 3, 5, 7))
      legend("topright", xpd = TRUE, inset = c(-0.04, 0.05), bty = "n",
             legend = c(expression("f"["LHS"]),
                        expression("f"["EI"]),
                        expression("f"["max(EI)"]),
                        expression("f"["opt"]),
                        expression("f"["unconverged"])),
             col = c("black", "blue", "green", "red", "orange"), pch = rep(19, 4))
      par(las = 0)
      dev.off()
      
      colors <- rep("red", N)
      for (i in 1:N) {
        if (state[i,] == "unconverged")
          colors[i] <- "orange"
      }
      
      open.dev(plot.type, "ei_convergence", N)
      par(las = 1)
      plot(c(0, seq(1, (N + 2))),
              c(min(f.temp), max(f.temp), min(Obj.opt.value), max(Obj.opt.value),
                rep(min(f.temp), N - 1)),
              xlab = "evaluations", ylab = "objective function",
              main = "EI optimization convergence",
              xaxt = "n", yaxt = "n", pch = 19, col = "white")
      lines(f.temp[1:N, , drop = FALSE], type = "p", pch = 19, col = colors)
      plotCI((N - iter.number + 1):N + 0.5, Obj.opt.value, 2 * sqrt(Obj.opt.var),
             pch = 19, col = "black", add = TRUE)
      legend("bottomleft", xpd = TRUE, inset = c(0, -0.5), bty = "n",
             legend = c("f", expression('f'['unconverged']),
                        expression('f'['min'])),
             col = c("red", "orange", "black"), pch = c(19, 19, 19))

      axis(1, at = NULL, tcl = 0.25, lwd = 2, las = 1)
      axis(2, at = NULL, tcl = 0.25, lwd = 2, las = 2)
      par(las = 0)
      custom.grid()
      
      dev.off()
    }
    
  } else if (imp.type == "REI") {
    
    if (L == 1) {
      
      x.min <- min(x.min, x0.min, x1.min)
      x.max <- max(x.max, x0.max, x1.max)
      
      N.p <- 200
      X.obj.fun <- matrix(seq(x.min, x.max, len = N.p), N.p, 1)
      f.obj.fun <- ObjectFunction(X.obj.fun)$val
      
      r.km.est <- km.est
      r.km.var <- km.var
      km.est <- KrigingModel(X[-(N + 1), , drop = FALSE], f.temp, theta, v, alpha)$EST
      
      f.est <- calc(km.est, X.obj.fun)
      mr.est <- calc(r.km.est, X.obj.fun)
      mr.var.est <- calc(r.km.var, X.obj.fun)
      
      mr.var.est <- NumZeroSubstitute(mr.var.est)
      # mr.var.est <- ifelse(abs(mr.var.est) < .Machine$double.eps, rep(0, N.p) , mr.var.est)
      
      mu.est <- calc(km.mu, X.obj.fun)
      
      rei.aug <- RelativeExpectedImprovementAug(r.km.est, InducedVariance,
                                               X[-(N + 1), , drop = FALSE],
                                               X.obj.fun,
                                               matrix(tail(Rei.x1.coord, 1), N.p, 1),
                                               M.inv, theta, v, alpha,
                                               tail(Obj.opt.value, 1))
      rei <- rei.aug$rei
      
      ## do sprawdzenia:
      ## Induced Variance
      mr.var.est <- rei.aug$var
      
      ## Mean value of the objective function
      N.conv <- 50
      x.conv <- matrix(seq(-0.5, 1, len = N.conv), N.conv, L)
      mr.obj.fun <- calc(RobustKrigingModel(x.conv, ObjectFunction(x.conv)$val,
                                            theta, v, alpha)$EST, x.conv)
      ##
      
      y.min <- min(rbind(f.obj.fun, f.est, mr.est, mr.est - 2 * sqrt(mr.var.est)))
      y.max <- max(rbind(f.obj.fun, f.est, mr.est, mr.est + 2 * sqrt(mr.var.est)))
      rei.max <- max(rei)
      rei <- rei / rei.max * (y.max - y.min) * 0.05
      
      colors <- rep("black", N)
      for (i in 1:N) {
        if (state[i,] == "unconverged")
          colors[i] <- "red"
      }
      
      open.dev(plot.type, "kriging", N)
      par(las = 1)
      matplot(X.obj.fun,
              cbind(f.obj.fun, f.est, mr.est, mr.est + 2 * sqrt(mr.var.est), rei, mu.est),
              type = "l", col = c(1, 2, 4, 3, 6, 5),
              lty = c(1, 1, 1, 2, 1, 1), lwd = c(2, 1.5, 1.5, 1.5, 1.5, 1.5),
              xlab = "", ylab = "", xaxt = "n", yaxt = "n",
              ylim = c(y.min, y.max))
      matplot(X.obj.fun, mr.est - 2 * sqrt(mr.var.est),
              type = "l", col = 3,
              lty = 2, lwd = 1.5,
              add = TRUE)
      
      ## do sprawdzenia
      matplot(x.conv, mr.obj.fun, type = "l", col = "black", lty = 2, lwd = 1.5, add = TRUE)
      ##
      
      polygon(c(X.obj.fun, rev(X.obj.fun)),
              c(mr.est + 2 * sqrt(mr.var.est), rev(mr.est - 2 * sqrt(mr.var.est))),
              col = rgb(0, 1, 0, 0.1), border = FALSE)
      abline(v = c(X[N + 1, ], tail(Rei.x1.coord, 1)), col = c("red", "blue"))
      custom.grid(lwd = 1, col = "gray")
      legend("bottomleft",
             c(expression('F'(x)), expression(widehat('F')(x)),
               expression(widehat('MR')(x)),
               expression(widehat('MR')(x)%+-%2*sqrt(widehat(sigma[MR])(x))),
               expression(REI(x)),
               expression(mu(x))),
             col = c(1, 2, 4, 3, 6, 5), lty = c(1, 1, 1, 2, 1, 1),
             lwd = c(2, 1.5, 1.5, 1.5, 1.5, 1.5), bty = "n")
      legend("bottomright",
             c(expression(x^'N+1'), expression(x['REI,max'])),
             col = c("red", "blue"), lty = c(1, 1), bty = "n")
      points(X[1:N, ], f.temp[1:N, ], pch = 19, col = colors)
      axis(1, at = NULL, tcl = 0.25, lwd = 2, las = 1)
      axis(2, at = NULL, tcl = 0.25, lwd = 2, las = 2)
      par(las = 0)
      dev.off()
      
    } else {
      
      colors <- c(rep("black", N - iter.number + 1),
                  rep("green", iter.number - 1),
                  "red", "blue", "yellow")
      for (i in 1:N) {
        if (state[i, ] == "unconverged")
          colors[i] <- "orange"
      }
      
      open.dev(plot.type, "kriging", N)
      par(las = 1)
      pairs(rbind(X, tail(Rei.x1.coord, 1), tail(Obj.opt.coord, 1)),
            labels = sapply(1:L, function(i) paste0("x_", i)), 
            main = "Points locations",
            col = colors,
            pch = 19, oma = c(5, 3, 5, 7))
      legend("topright", xpd = TRUE, inset = c(-0.04, 0.05), bty = "n",
             legend = c(expression("x"["LHS"]),
                        expression("x"["REI"]),
                        expression("x"["unconverged"]),
                        expression("x"["max(REI)"]),
                        expression("x"^"N+1"),
                        expression("x"["opt(ME)"])),
             col = c("black", "green", "orange", "blue", "red", "yellow"),
             pch = 19)
      par(las = 0)
      dev.off()
      
      colors <- rep("red", N)
      for (i in 1:N) {
        if (state[i,] == "unconverged")
          colors[i] <- "orange"
      }
      
      open.dev(plot.type, "rei_convergence", N)
      par(las = 1)
      plot(c(0, seq(1, (N + 2))),
           c(min(f.temp), max(f.temp), min(Obj.opt.value), max(Obj.opt.value),
             rep(min(f.temp), N - 1)),
           xlab = "evaluations", ylab = "objective function, mean response",
           main = "REI optimization convergence",
           xaxt = "n", yaxt = "n", pch = 19, col = "white")
      lines(f.temp[1:N, , drop = FALSE], type = "p", pch = 19, col = colors)
      plotCI((N - iter.number + 1):N + 0.5, Obj.opt.value, 2 * sqrt(Obj.opt.var),
             pch = 19, col = "black", add = TRUE)
      legend("bottomleft", xpd = TRUE, inset = c(0, -0.5), bty = "n",
             legend = c("f",
                        expression('f'['unconverged']),
                        expression('ME'['opt'])),
             col = c("red", "orange", "black"), pch = c(19, 19, 19))
      
      axis(1, at = NULL, tcl = 0.25, lwd = 2, las = 1)
      axis(2, at = NULL, tcl = 0.25, lwd = 2, las = 2)
      par(las = 0)
      custom.grid()
      
      dev.off()
    }

  } else {
    stop("Wrong \"imp.type\" value.")
  }
}

krig.settings <- list(CovarianceType = CovarianceGauss, poly.deg = 1,
                      par.fun.type = "LogLikelihood", eigen.rat = 1e6,
                      opt.search.type = "genetic", gen.enhance = TRUE,
                      gen.size = 140, pop.size = 300,
                      lhs.point.number = 100, lhs.add.points = NULL,
                      lhs.prev.points = TRUE,
                      lhs.enhance = TRUE, lhs.enhance.number = 0.05,
                      theta.min = 0.1, theta.max = 1.0,
                      alpha.min = 0.90, alpha.max = 0.9999)

optim.settings <- list(find.minimum = TRUE,
                       opt.search.type = "genetic", gen.enhance = TRUE,
                       gen.size = 140, pop.size = 300,
                       lhs.point.number = 200, lhs.add.points = NULL,
                       lhs.enhance = TRUE, lhs.enhance.number = 0.05)

ei.settings <- list(imp.type = "Expected Improvement",
                    find.minimum = FALSE,
                    opt.search.type = "genetic", gen.enhance = TRUE,
                    gen.size = 140, pop.size = 300,
                    lhs.point.number = 100, lhs.add.points = NULL,
                    lhs.enhance = TRUE, lhs.enhance.number = 0.05,
                    threshold = 0.005)

rei.settings <- list(imp.type = "Relative Expected Improvement",
                     find.minimum = FALSE,
                     x0.min = NULL, x0.max = NULL,
                     x1.min = NULL, x1.max = NULL,
                     opt.search.type = "genetic", gen.enhance = TRUE,
                     gen.size = 140, pop.size = 300,
                     lhs.point.number = 100, lhs.add.points = NULL,
                     lhs.enhance = TRUE, lhs.enhance.number = 0.05,
                     threshold = 0.005)

uncertainty.settings <- list(ErrorType = ErrorGauss, st.dev = 1, delta = 1)

