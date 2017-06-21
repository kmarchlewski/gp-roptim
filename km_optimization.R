#!/home/kmarchlewski/Programs/bin/Rscript

# source("km_optimization.R")

library(gaussAlgebra)
library(fields)
library(nloptr)
library(mco)
library(lhs)
library(rgl)
library(plotrix)

source("02_solvers_scripts/Hartman_6d_test_fun/hartmann_6d_solver.R")
source("01_kriging_model_optimization/km_optim_data.R")
source("03_design_of_experiment/doe.R")
source("configuration/settings.R")

assign("last.warning", NULL, envir = baseenv())
options(nwarnings = 1000, warn = 1)

# X.aug <- DoE(N, L, x.min, x.max, add.points = NULL, method = "geneticLHS")

X.aug <- read.table(file = "doe_data.csv", header = TRUE,
                    sep = ";", dec = ".")
X.aug <- list(coord = as.matrix(X.aug[, 1:L]),
              val = as.matrix(X.aug[, L + 1]),
              state = as.matrix(X.aug[, L + 2]))

rei.optim <- REIoptim(X.aug, x.min, x.max, iter.number, ObjectFunction,
                      krig.settings, optim.settings, rei.settings,
                      enable.cv = FALSE, write.to.file = TRUE, continue = FALSE)

# ei.optim <- EIoptim(X.aug, x.min, x.max, iter.number, ObjectFunction,
#                     krig.settings, optim.settings, ei.settings,
#                     enable.cv = FALSE, write.to.file = TRUE, continue = FALSE)

X <- ei.optim$X
f <- ei.optim$f
state <- ei.optim$state

Theta <- ei.optim$Theta
Alpha <- ei.optim$Alpha
V <- ei.optim$V

Obj.opt.value <- ei.optim$Obj.opt.value
Obj.opt.var <- ei.optim$Obj.opt.var
Obj.opt.coord <- ei.optim$Obj.opt.coord

kriging.model <- KrigingModel(X, f, Theta[iter.number,], V[iter.number],
                              Alpha[iter.number])
km.est <- kriging.model$EST
km.var <- kriging.model$VAR

# for (i in 1:N) {
# FunctionProjection(ObjectFunction, km.est, X, f, x.min, x.max,
#                    X[which.min(f),], c(2, 3))
#   Sys.sleep(2)
# }
