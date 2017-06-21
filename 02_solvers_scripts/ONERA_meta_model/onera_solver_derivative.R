
library(gvector)
load(paste0("~/0_cloud_storage/01_projects/03_R/00_share/",
            "02_solvers_scripts/ONERA_meta_model/Metamodel.Rdata"))

ObjectDerivative <- function (x, vectorized = FALSE) {
  eps <- 1e-6
  
  FunDer <- function (X, Eps) {
    eps <- Eps[1,1]
    der <- -(t(calc(Fun, X + Eps)) - t(calc(Fun, X - Eps))) / (2 * eps)
    return (der)
  }
  
  if (vectorized == TRUE) {
    L <- ncol(x)
    N <- nrow(x)
    X <- do.call(rbind, lapply(1:L, function(i) x))
    Eps <- diag(eps, nrow = L)
    Eps <- sapply(1:L, function (i) t(sapply(1:N, function(n) Eps[i, , drop = FALSE])) )
    
    val <- matrix(FunDer(X, Eps), N, L)
    state <- matrix("converged", N, L)
  } else {
    L <- length(x)
    N <- 1
    X <- do.call(rbind, lapply(1:L, function(i) x))
    Eps <- diag(eps, nrow = L)
    
    val <- matrix(FunDer(X, Eps), 1, L)
    state <- matrix("converged", 1, L)
  }
  return (list(val = val, state = state))
}

ObjectFunction <- function (x, vectorized = FALSE) {
  if (vectorized == TRUE) {
    val <- -t(calc(Fun, x))
    state <- matrix("converged", nrow(x), 1)
  } else {
    val <- -calc(Fun, t(x))
    state <- "converged"
  }

  return (list(val = val, state = state))
}
