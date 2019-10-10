## Library for optimization on manifold
## Written by Kohei Yoshikawa
library(MASS)

sym <- function(X){
  return(0.5*(t(X) + X))
}

## Stiefel Manifold : St(p,n) = {X in R^(n x p) : X^T X = I_p}
## generate matrix on Stiefel Manifold
Stiefel.rnorm <- function(p, n, mu=numeric(p), sigma=diag(p)){
  X <- mvrnorm(n, mu, sigma)
  qd <- qr(X)
  return(qr.Q(qd))
}

# innner product
Stiefel.inner <- function(xi, eta){
  return(tr(t(xi)%*%eta))
}

# projection onto Stifel Manifold
Stiefel.proj <- function(X, xi){
  return(xi - X%*%sym( (t(X)%*%xi) ) )
}

# Euclid gradient to Riemannian gradient
Stiefel.eugrad2rgrad <- function(X, eta){
  Stiefel.proj(X, eta)
}

# Retraction to Stiefel Manifold
Stiefel.retract <- function(X, xi, parameters){
# Stiefel.retract <- function(X, xi, parameters){
  #qd <- spQR(X + xi)
  qd <- qr(X + xi)
  # qd <- qr(rslre.Round2zero(X + xi, parameters))
  #return(qd$Q)
  return(qr.Q(qd))
}
