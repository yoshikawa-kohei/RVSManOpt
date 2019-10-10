m_admm <- function(X, Y, W, G, U, D, V, U.aster, V.aster, V.aster2, Omega, Phi, Psi, parameters){
  
  hist <- NULL
  hist$Cost <- numeric(parameters$epoch)
  n <- dim(X)[1]
  
  p <- ncol(X)
  q <- ncol(Y)
  
  sqrt_G <- sqrtm(G)
  inv_sqrt_G <- solve(sqrt_G)
  
  for (itr in 1:parameters$epoch){
    ######################################################################
    ## U-Step ------------------------------------------------------------
    # Calculate U gradient
    U.eugrad <- calcGradientU(X, Y, U, D, V, U.aster, Omega, parameters)
    U.rgrad <- GStiefel.eugrad2rgrad(U, U.eugrad, G)
    # calculate step on Armijo conditon
    step <- calcStep(X, Y, U, D, V, U.rgrad, G, sqrt_G, inv_sqrt_G, parameters, "U")
    # Update U
    U <- GStiefel.retract(U, -step * U.rgrad, G, sqrt_G, inv_sqrt_G, parameters)
    
    ######################################################################
    ## V-Step ------------------------------------------------------------
    # Calculate V gradient
    V.eugrad <- calcGradientV(X, Y, U, D, V, V.aster, V.aster2, Phi, Psi, parameters)
    V.rgrad <- Stiefel.eugrad2rgrad(V, V.eugrad)
    # calculate step on Armijo conditon
    step <- calcStep(X, Y, U, D, V, V.rgrad, G, sqrt_G, inv_sqrt_G, parameters, "V")
    # Update V
    V <- Stiefel.retract(V, -step * V.rgrad, parameters)
    
    ######################################################################
    ## D-Step ------------------------------------------------------------
    D.old <- D
    D <- updateD(X, Y, U, D, V)
    
    ######################################################################
    ## U*-Step -----------------------------------------------------------
    U.aster.old <- U.aster
    U.aster <- matrix_soft_thresold(U + Omega, W$U, n*parameters$lambda2/parameters$gamma2)

    ######################################################################
    ## V*-Step -----------------------------------------------------------
    V.aster.old <- V.aster
    V.aster <- matrix_soft_thresold(V + Phi, W$V, n*parameters$lambda1*parameters$alpha/parameters$gamma3)

    ######################################################################
    ## V**-Step (rank selection) -----------------------------------------
    V.aster2.old <- V.aster2
    V.aster2 <- matrix_group_hard_thresold(V + Psi, W$D, sqrt(2*n*parameters$lambda1*(1-parameters$alpha)*sqrt(q)/parameters$gamma1))
    
    ######################################################################
    ## Omega, Phi, Pshi-Step ------------------------------------------------------------
    Omega <- Omega + U - U.aster
    Phi <- Phi + V - V.aster
    Psi <- Psi + V - V.aster2
    
    ######################################################################
    ## Calculate Cost and Check stopping rule ----------------------------
    hist$Cost[itr] <- calc.cost(X, Y, U, D, V.aster2)
    
    if ((norm(U.aster.old, "F")!=0) & (norm(V.aster.old, "F")!=0) & (norm(V.aster2.old, "F")!=0)){
      if (norm(U.aster - U.aster.old, "F")/norm(U.aster.old, "F")
          + norm(V.aster - V.aster.old, "F")/norm(V.aster.old, "F")
          + norm(V.aster2 - V.aster2.old, "F")/norm(V.aster2.old, "F")
          < 1e-5){
        cat("Iteration = ", itr, "\n")
        break
      }
    }
    
  }
  
  # Select estimated factors
  selected.factor.indexes <- which(apply(V.aster2, 2, function(x){norm(x, type="2")}) != 0)
  
  if(length(selected.factor.indexes) == 0){
    U.est <- NULL
    V.est <- NULL
    D.est <- NULL
    C.est <- NULL
    SSE <- tr(Y %*% t(Y))/(n*q)
    df <- -1
    GIC <- log(SSE) + (log(q*n)/(n*q)) * df
    BIC <- log(SSE) + (df*log(q*n))/(q*n)
    BICP <- log(SSE) + 2*df*log(p*q)/(q*n)
    AIC <- log(SSE) + 2/q/n*(df)
  }else{
    U.est <- matrix(0, nrow=p,ncol=ncol(U))
    V.est <- matrix(0, nrow=q,ncol=ncol(V))
    U.est[U.aster!=0] <- U[U.aster!=0]
    V.est[V.aster!=0] <- V[V.aster!=0]
    U.est <- as.matrix(U.est[, selected.factor.indexes])
    V.est <- as.matrix(V.est[, selected.factor.indexes])
    D.est <- as.matrix(D[selected.factor.indexes,selected.factor.indexes])
    C.est <- U.est %*% D.est %*% t(V.est)
    # Calclate information criterion(BIC)
    n <- nrow(Y); p <- ncol(X); q <- ncol(Y);
    SSE <- (calc.cost(X, Y, U.est, D.est, V.est)^2)/(n*q)
    df <- (length(which(U.est != 0)) + length(which(V.est != 0)) -1)
    if(p < n){
      GIC <- log(SSE) + ((log(log(n*q))*log(p*q))/(n*q)) * df
      BIC <- log(SSE) + (df*log(q*n))/(q*n)
      BICP <-  log(SSE) + 2*df*log(p*q)/(q*n)
      AIC <- log(SSE) + 2/q/n*(df)
    }else if(p >= n){
      GIC <- log(SSE) + ((log(log(n*q))*log(p*q))/(n*q)) * df
      BIC <- log(SSE) + (df*log(q*n))/(q*n)
      BICP <- log(SSE) + 2*df*log(p*q)/(q*n)
      AIC <- log(SSE) + 2/q/n*(df)
    }
  }
  
  return(list(U = U.est, D = D.est, V = V.est, GIC = GIC, BIC = BIC, BICP=BICP, AIC=AIC, C.est = C.est, parameters=parameters))
}

######################################################################
## Functions
######################################################################
# U-step U gradient
calcGradientU <-function(X, Y, U, D, V, U.aster, Omega, parameters){
  grad <- - t(X) %*% Y %*% V %*% D + parameters$gamma1 * (U - U.aster + Omega)
  return(grad/norm(grad))
}

# V-step V gradient
calcGradientV <-function(X, Y, U, D, V, V.aster, V.aster2, Phi, Psi, parameters){
    grad <- -(t(Y) %*% X %*% U %*% D) + parameters$gamma1 * (V - V.aster2 + Psi) + parameters$gamma3 * (V- V.aster + Phi)
    return(grad/(norm(grad)))
}

updateD <- function(X, Y, U, D, V){
  VYXU = t(V) %*% t(Y) %*% X %*% U
  n = dim(X)[1];
  D = diag(diag(VYXU/n))
  return(D)
}

matrix_soft_thresold <- function(X, W, lambda){
  softX <- sign(X)*pmax((abs(X) - W*(lambda)),0)
  return(softX)
}

matrix_group_hard_thresold <- function(X, W, lambda){
  matrix_dim <- dim(X)
  weight <- diag(W)
  for(i in 1:matrix_dim[2]){
    if( 1 - (lambda*sqrt(weight[i]))/norm(X[, i], type="2") < 0 ){
      X[, i] <- 0
    }
  }
  return(X)
}

calcStep <- function(X, Y, U, D, V, rgrad, G, sqrt_G, inv_sqrt_G, parameters, type){
  m <- 0
  step <- parameters$step*(parameters$contraction.factor^m)
  hist_m_cost <- numeric(parameters$max.ls.steps)
  
  if(type == "U"){
    nextU <- GStiefel.retract(U, -step*rgrad, G, sqrt_G, inv_sqrt_G, parameters)
    nextCost <- calc.cost(X, Y, nextU, D ,V)
    currentCost <- calc.cost(X, Y, U, D, V)
    
    while(nextCost > currentCost - parameters$sigma*step*tr(t(rgrad)%*%G%*%(rgrad))){
      # Reduce the step size
      m <- m + 1
      step <- parameters$step*(parameters$contraction.factor^m)
      nextU <- GStiefel.retract(U, -step*rgrad, G, sqrt_G, inv_sqrt_G, parameters)
      nextCost <- calc.cost(X, Y, nextU, D, V)
      currentCost <- calc.cost(X, Y, U, D, V)
      
      if(m > parameters$max.ls.steps){
        break
      }
    }
  }else if(type == "V"){
    nextV <- Stiefel.retract(V, -step*rgrad, parameters)
    nextCost <- calc.cost(X, Y, U, D, nextV)
    currentCost <- calc.cost(X, Y, U, D, V)
    
    while(nextCost > currentCost - parameters$sigma*step*tr(t(rgrad)%*%(rgrad))){
      # Reduce the step size
      m <- m + 1
      step <- parameters$step*(parameters$contraction.factor^m)
      nextV <- Stiefel.retract(V, -step*rgrad, parameters)
      nextCost <- calc.cost(X, Y, U, D, nextV)
      currentCost <- calc.cost(X, Y, U, D, V)
      
      if(m > parameters$max.ls.steps){
        break;
      }
    }
  }
  
  return(step)
}

calc.cost <- function(X, Y, U, D, V){
  return(norm(Y - X %*% U %*% D %*% t(V), type="F"))
}
