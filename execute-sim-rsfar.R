library(secure)
library(stringr)
source("R/rvsmanopt.R")
source("R/generate-dataset.R")

# generate dataset ------------------------------------------------------------
n <- 400; p <- 120; q <- 60; nrank <- 12
s2n <- 0.5; rho <- 0.3
dataset <- generate.data(n, p , q, nrank, s2n, rho)

# load data from dataset --------------------------------------------------
rx <- qr(dataset$X)$rank
max.rank.est <- min(rx, dataset$q)

parameters <- rvs.manopt.parameters()
fit.rvs.manopt <- rvs.manopt.path(dataset$X, dataset$Y, nrank=max.rank.est, parameters)
