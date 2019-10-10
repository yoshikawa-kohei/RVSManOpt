# Generate synthentic testing data --------------------------------------------
# Simulation setting: 
# n: sample size 
# p: the number of dimension of the predictor vector
# q: the number of dimension of the response vector 
# nrank: true rank of the model
# s2n: signal to noise ratio
# rho: error correlation

generate.data <- function(n, p , q, nrank, s2n, rho){
  # Generate data
  U <- matrix(0, ncol = nrank, nrow = p)
  V <- matrix(0, ncol = nrank, nrow = q)
  D <- rep(0, nrank)
  
  u <- c(1, -1, 1, -1, 0.5, -0.5, 0.5, -0.5)
  v <- c(1, -1, 0.5, -0.5)
  
  U[, 1] <- c(u, rep(0, p-8))
  if(nrank != 1){
    for(i in 2:nrank){
      U[, i] <- c(rep(0,5*(i-1)), u, rep(0, p-(5*i+3) ) )
    }
  }
  
  V[, 1] <- c(v,rep(0,q-4))
  if(nrank != 1){
    for(i in 2:nrank){
      V[, i] <- c(rep(0,4*(i-1)), v, rep(0, q-4*i))
    }
  }
  if(nrank != 1){
  U[,1:nrank] <- apply(U[,1:nrank],2,function(x)x/sqrt(sum(x^2)))
  V[,1:nrank] <- apply(V[,1:nrank],2,function(x)x/sqrt(sum(x^2)))
  D[1:nrank] <- seq((5.0+(nrank-1)*0.1), 5.0, by=-0.1); D <- diag(D,nrow=length(D))
  }else{
    U <- U/sum(U^2)
    V <- V/sum(V^2)
    D <- as.matrix(5)
  }
  
  C <- U%*%D%*%t(V)
  
  Xsigma <- 0.5^abs(outer(1:p, 1:p,FUN="-"))
  
  sim.sample <- secure.sim(U, D, V, n, s2n, Xsigma, rho)
  Y <- sim.sample$Y[1:n,]
  X <- sim.sample$X[1:n,]  
  
  return(dataset = list(X=X, Y=Y, U=U, D=D, V=V, C=C, Xsigma=Xsigma,
                        nrank=nrank, n=n, p=p, q=q,
                        s2n=s2n, rho=rho))
}

