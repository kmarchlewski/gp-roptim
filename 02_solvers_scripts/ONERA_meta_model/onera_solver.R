
library(gvector)
load(paste0("~/0_cloud_storage/01_projects/03_R/00_share/",
            "02_solvers_scripts/ONERA_meta_model/Metamodel.Rdata"))

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
