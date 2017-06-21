
# Settings: Hartmann 6d model 25.05.2016
# population size for nsga2 algorithm must be multiple of 4

L <- 6
N <- 30
x.min <- 0
x.max <- 1
x.min <- rep(x.min, L)
x.max <- rep(x.max, L)

krig.settings$CovarianceType <- CovarianceGauss # CovarianceMatern
krig.settings$par.fun.type <- "LogLikelihood" # "LOOCV"
krig.settings$poly.deg <- 1
krig.settings$theta.min <- rep(0.05, L)
krig.settings$theta.max <- rep(0.5, L)
krig.settings$alpha.max <- 0.99999
krig.settings$alpha.min <- 0.90
krig.settings$gen.size <- 100 + 40 * floor(2 * log(L))
krig.settings$pop.size <- 160 + 40 * floor(2 * log(L))
krig.settings$opt.search.type <- "disabled" # "gradient" 
krig.settings$lhs.point.number <- 150 * (L + 1) + 1
krig.settings$lhs.prev.points <- TRUE

optim.settings$find.minimum <- TRUE
optim.settings$gen.size <- 100 + 40 * floor(2 * log(L))
optim.settings$pop.size <- 160 + 40 * floor(2 * log(L))
optim.settings$opt.search.type <- "gradient"
optim.settings$lhs.point.number <- 100 * (L + 1) + 1

ei.settings$gen.size <- 240 + 40 * floor(2 * log(L))
ei.settings$pop.size <- 320 + 40 * floor(2 * log(L))
ei.settings$opt.search.type <- "gen-EI-2" # "gen-EI-1" "gen-EI-2" "gradient" "genetic"
ei.settings$lhs.point.number <- 150 * (L + 1) + 1

rei.settings$x0.min <- x.min
rei.settings$x0.max <- x.max
rei.settings$x1.min <- x.min
rei.settings$x1.max <- x.max
rei.settings$gen.size <- 240 + 40 * floor(2 * log(L))
rei.settings$pop.size <- 320 + 40 * floor(2 * log(L))
rei.settings$opt.search.type <- "genetic" # "gen-EI-1" "gen-EI-2" "gradient" "genetic"
rei.settings$lhs.point.number <- 400 * (L + 1) + 1

uncertainty.settings$ErrorType <- ErrorGauss # ErrorGauss ErrorBeta22
uncertainty.settings$st.dev <- rep(0.04, L)
uncertainty.settings$delta <- rep(0.04, L)

iter.number <- 40
plot.type <- "png"
