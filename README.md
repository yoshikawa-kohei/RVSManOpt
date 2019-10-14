# RVSManOpt
Authors: Kohei Yoshikawa and Shuichi Kawano

# Introduction 
This is an implementation of RVSManOpt in "Sparse Reduced-Rank Regression for Simultaneous Rank and Variable Selection via Manifold Optimization"

# Usage example: Demonstraition of synthetic dataset 
This is a simple example which shows how to use these codes.
## 1. Import some libraries and our source codes.
```R
library(secure)
source("R/rvsmanopt.R")
source("R/generate-dataset.R")
```

## 2. Generate a dataset
Set the following parameters, and generate dataset.
  - n: Sample size 
  - p: The dimension of the predictor vectors
  - q: The dimension of the response vectors
  - nrank: The true rank of the model
  - s2n: Signal to noise ratio
  - rho: The correlation of error vectors

```R
n <- 400; p <- 120; q <- 60; nrank <- 12
s2n <- 0.5; rho <- 0.3
dataset <- generate.data(n, p , q, nrank, s2n, rho)
```

## 3. Perform RVSManOpt
Now, you can perform the RVSManOpt.
```R
rx <- qr(dataset$X)$rank
max.rank.est <- min(rx, dataset$q)

parameters <- rvs.manopt.parameters()
fit.rvs.manopt <- rvs.manopt.path(dataset$X, dataset$Y, nrank=max.rank.est, parameters)
```

# Licence
These code are free, non-commercial and open source.
