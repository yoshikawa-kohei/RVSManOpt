######################################################################
# The factor extraction algorithm with rank and variable selection via sparse regularization and manifold optimization (RVSManOpt).
# 
######################################################################
library(MASS)
source("func/manifold/Stiefel.R")
source("func/manifold/GeneralizedStiefel.R")
source("R/m-admm.R")

rvs.manopt.parameters <- function(lambda.max=1, lambda.min=1e-15, lambda.split=100, lambda1=1, lambda2=1,
                              step=1, contraction.factor=0.5, max.ls.steps=10, sigma=0.5, epoch=2000,
                              gamma1=1, gamma2=1, gamma3=1, alpha=0.05, U.gamma=1, V.gamma=1, D.gamma=10){
  return(list(lambda.max=lambda.max, lambda.min=lambda.min, lambda.split=lambda.split, lambda1=lambda1, lambda2=lambda2,
              step=step, contraction.factor=contraction.factor, max.ls.steps=max.ls.steps, sigma=sigma, epoch=epoch,
              gamma1=gamma1, gamma2=gamma2, gamma3=gamma3, alpha=alpha, U.gamma=U.gamma, V.gamma=V.gamma, D.gamma=D.gamma))
}

rvs.manopt <- function(X, Y, nrank, parameters, init.U=NULL, init.D=NULL, init.V=NULL){
  n <- nrow(Y); p <- ncol(X); q <- ncol(Y);
  G <- (t(X) %*% X) / n

  # Generate weight matrixes W_U, W_V, W_D ----------------------------------
  fit<-NULL
  fit.eig <- eigen((t(Y) %*% X %*% ginv(t(X) %*% X) %*% t(X) %*% Y)/n)
  fit$V <- fit.eig$vectors[, 1:nrank]
  fit$D <- diag(sqrt(fit.eig$values[1:nrank]))
  fit$U <- (ginv(t(X) %*% X) %*% t(X) %*% Y %*% fit.eig$vectors%*% ginv(diag(sqrt(fit.eig$values))))[, 1:nrank]
  
  W <- list()
  W$U <- 1 / abs(fit$U)
  W$V <- 1 / abs(fit$V)
  W$D <- diag(1 / diag(fit$D))
  
  W$U <- W$U^parameters$U.gamma
  W$V <- W$V^parameters$V.gamma
  W$D <- W$D^parameters$D.gamma


  # Initialize estimators ---------------------------------------------------
  if(is.null(init.U)){
    U <- fit$U
  }else{
    zero.mat <- matrix(0, nrow=p, ncol=nrank-dim(init.U)[2])
    U <- cbind(init.U, zero.mat)
  }
  if(is.null(init.D)){
    D <- fit$D
  }else{
    D <- diag(c(diag(init.D),rep(0, nrank-dim(init.D)[2])))
  }
  if(is.null(init.V)){
    V <- fit$V
  }else{
    zero.mat <- matrix(0, nrow=q, ncol=nrank-dim(init.V)[2])
    V <- cbind(init.V, zero.mat)
  }
  
  # Initialize dual variables
  U.aster <- U
  V.aster <- V
  V.aster2 <- V
  Omega <- matrix(0, nrow=dim(U)[1], ncol=dim(U)[2])
  Phi <- matrix(0, nrow=dim(V)[1], ncol=dim(V)[2])
  Psi <- matrix(0, nrow=dim(V)[1], ncol=dim(V)[2])
  
  # Solve problem with Manifold-ADMM ----------------------------------------
  fit.final <- m_admm(X, Y, W, G, U, D, V, U.aster, V.aster, V.aster2, Omega, Phi, Psi, parameters)
  if(!is.null(fit.final$C.est)){
    if(length(fit.final$D) == 1){
      if(fit.final$D<0){
        fit.final$D <- -fit.final$D
        fit.final$U <- -fit.final$U
      }
    }else{
      idx <- sort(abs(diag(fit.final$D)),T, index.return=TRUE)$ix
      fit.final$U <- fit.final$U[,idx]
      fit.final$V <- fit.final$V[,idx]
      fit.final$D <- diag(diag(fit.final$D)[idx])
      for(i in 1:dim(fit.final$D)[1]){
        if(fit.final$D[i, i]<0){
          fit.final$D[i, i] <- -fit.final$D[i, i]
          fit.final$U[, i] <- -fit.final$U[, i]
        }
      }
    }
  }
  return(fit.final)
}

rvs.manopt.path <- function(X, Y, nrank, parameters){
  cat("[info] Start the estimation... \n")
  lambda1.seq <- seq(parameters$lambda.max, parameters$lambda.min, length=parameters$lambda.split)
  lambda2.seq <- seq(parameters$lambda.max, parameters$lambda.min, length=parameters$lambda.split)
  
  conditions <- expand.grid(lambda1=lambda1.seq, lambda2=lambda2.seq)
  list.fit.rvs.manopt <- list()
  array.IC <- numeric(nrow(conditions))
  
  # BIC
  for(itr in 1:nrow(conditions)){
    parameters$lambda1 <- conditions$lambda1[itr]
    parameters$lambda2 <- conditions$lambda2[itr]
    list.fit.rvs.manopt[[itr]] <- rvs.manopt(X, Y, nrank, parameters)
    array.IC[itr] <- list.fit.rvs.manopt[[itr]]$BIC
  }
  
  # Find the index of the fitted model that minimizes information criteria value
  min.idx <- which.min(array.IC)
  cat("[info] Finish the estimation! \n")
  return(list.fit.rvs.manopt[[min.idx]])
}