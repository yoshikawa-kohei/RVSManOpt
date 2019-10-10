## Library for optimization on Generalized Stiefel manifold
## Written by Kohei Yoshikawa
library(MASS)

sym <- function(X){
  return(0.5*(t(X) + X))
}

tr <- function(X){
  return(sum(diag(X)))
}

# matrix square root
sqrtm <- function(X){
  U <- svd(X)$u
  V <- svd(X)$v
  D <- diag(sqrt(svd(X)$d))
  return(U %*% D %*% t(V))
}

## Generalized Stiefel Manifold : GSt(p,n) = {X in R^(n x p) : X^T G X = I_p} for some G > 0
## generate matrix on Generalized Stiefel Manifold
GStiefel.rnorm <- function(q, p, G, mu=numeric(q), sigma=diag(q)){
  X = mvrnorm(p, mu, sigma)
  U <- svd(X)$u
  V <- svd(X)$v
  L <- diag(eigen(t(U) %*% G %*% U)$values)
  P <- eigen(t(U) %*% G %*% U)$vectors
  return(U%*%P%*%solve(sqrtm(L)))
}

# innner product on tangent space
GStiefel.inner <- function(xi, eta, G){
  return(tr(t(xi) %*% G %*% eta))
}

# projection onto Generalized Stifel Manifold
GStiefel.proj <- function(X, xi, G){
  return(xi - X %*% sym(t(X) %*% G %*% xi))
}

# Euclid gradient to Riemannian gradient
GStiefel.eugrad2rgrad <- function(X, eta, G){
  return(GStiefel.proj(X, eta, G))
}

# Retraction to Generalized Stiefel Manifold
# GStiefel.retract <- function(X, xi, G, sqrt_G, inv_sqrt_G, parameters){
#   qd <- qr(sqrt_G %*% (X + xi))
#   return(inv_sqrt_G %*% qr.Q(qd))
# }

GStiefel.retract <- function(X, xi, G, sqrt_G, inv_sqrt_G, parameters){
  Z <- t(X + xi) %*% G %*% (X + xi)
  R <- chol(Z)
  return((X+xi) %*% backsolve(r = R, x = diag(ncol(R))))
}
