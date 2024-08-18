library(R2jags)
load.module("lecuyer")
load.module("glm")
load.module("dic")

IRT_LSM.string <- "model{
  for (i in 2:N){
   theta[i-1] <- beta0 + inprod(edgeX[i, 1, ], beta.edge) + Z[i, 1]*Z[1, 1] + Z[i, 2]*Z[1, 2]
  }

  for (j in 2:(N-1)) {
    for (i in 1:(j-1)) {
      theta[i + (j-1)*(N-1)] <- beta0 + inprod(edgeX[i, j, ], beta.edge) + Z[i, 1]*Z[j, 1] + Z[i, 2]*Z[j, 2]   
    }
    for (i in (j+1):N){
      theta[(i-1) + (j-1)*(N-1)] <- beta0 + inprod(edgeX[i, j, ], beta.edge) + Z[i, 1]*Z[j, 1] + Z[i, 2]*Z[j, 2]    
    }
  }
  
  for (i in 1:(N-1)){
    theta[(N-1)*(N-1) + i] <- beta0 + inprod(edgeX[i, N, ], beta.edge) +  Z[i, 1]*Z[N, 1] + Z[i, 2]*Z[N, 2]
  }
  
  for (i in 1:(N*(N-1))){
    for (j in 1:J){
      logit(p[i, j]) <- int[j] + slp[j] * theta[i]
      Y[i, j] ~ dbern(p[i, j])
    }
  }

  # prior
  slp[1] = 1
    for (i in 2:J){
      slp[i] ~ dt(0, pow(1.25, -2), 1)T(0, )
    }
  
  for (i in 1:J){
      int[i] ~ dt(0, pow(10, -2), 1)
  }
  
  beta0 = 0
  for (i in 1:n_cov){
    beta.edge[i] ~ dnorm(0, 0.1)
  }
  
  for (i in 1:N){
      Z[i, 1:2] ~ dmnorm.vcov(Z_mu, Z_vcov[1:2, 1:2])
  }

  Z_mu[1] <- 0
  Z_mu[2] <- 0
  Z_vcov[1, 1] <- sigsqz
  Z_vcov[2, 2] <- sigsqz
  Z_vcov[2, 1] <- 0
  Z_vcov[1, 2] <- 0
  
  sigsqz ~ dt(0, 1, 4)T(0,) # half t prior
  
}"
