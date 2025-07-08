##-----------------------------------------#
## [item] 2PL [rater] 1PL model 
##-----------------------------------------#
library(nimble)

modelCode <- nimbleCode({
  # Likelihood
  for(n in 1:tot){
    
    y[n] ~ dcat(ppi[n,1:K])
    
    pic[n,1]              <- 0
    pi[n,1]               <- 0
    

    for(k in 2:K){
      pic[n,k]          <-  trans_alpha_r[RRi[n]]*alpha_i[II[n]] * (eta[PPi[n]] - beta_ir[II[n],RRi[n]] - category_prm[II[n], RRi[n], k])  
      ## cumumlative sum
      pi[n,k]           <-  sum(pic[n,1:k])     
    }
    ppi[n,1:K]          <- exp(pi[n,1:K])/sum(exp(pi[n,1:K]))
  }
  
  ## Individual effects - ability
  
  for(p in 1:P) {
    eta[p] ~ dnorm(0, var = 1)  
  }

  # Item discrimination
  for(i in 1:I){
    alpha_i[i] ~ dlnorm(0, 0.5)
  }
  
  # Rater discrimination (with constraints)
  trans_alpha_r[1] <- 1 / prod(alpha_r[1:(R-1)])
  for(r in 2:R){
    trans_alpha_r[r] <- alpha_r[r-1]
  }
  
  for(r in 1:(R-1)){
    alpha_r[r] ~ dlnorm(0, 0.5)
  }
  
  # Item-rater difficulty
  for(i in 1:I){
    for(r in 1:R){
      beta_ir[i, r] ~ dnorm(0, var = 1)
    }
  }

  # Rater category thresholds
  for(r in 1:R){
    for(k in 1:(K-1)){
      beta_rk[r, k] ~ dnorm(0, var = 1)
    }
    category_est_r[r, 1:(K-1)] <- beta_rk[r, 1:(K-1)]
    category_est_r[r, K-1]     <- -sum(beta_rk[r, 1:(K-1)])
  }
  
  # Item category thresholds
  for(i in 1:I) {
    for(k in 1:(K-1)) {
      beta_ik[i, k] ~ dnorm(0, var = 1)
    }
  
  category_est_i[i, 1:(K-1)] <- beta_ik[i, 1:(K-1)]
  category_est_i[i, K-1]     <- -sum(beta_ik[i, 1:(K-1)])
    
    for(r in 1:R){
      category_prm[i, r, 1] <- 0
      for(k in 2:K){
        category_prm[i, r, k] <- category_prm[i, r, k-1] + 
        category_est_r[r, k-1] + category_est_i[i, k-1]
      }
    }
  }
})


constants <- list(I = max(Data$II), P = max(Data$PPi),R = max(Data$RRi), K=max(Data$y), 
                  tot = length(Data$y), II=Data$II , PPi=Data$PPi, RRi=Data$RRi)
data <- list(y = Data$y) #check

inits <- list(eta    = rnorm(constants$P, 0, 3),
              alpha_i = rep(1, constants$I),
              alpha_r = rep(1, constants$R), 
              beta_ir = matrix(rnorm(constants$I * constants$R), nrow = constants$I ), 
              beta_rk = matrix(rnorm(constants$R * (constants$K - 1)), nrow = constants$R ),
              beta_ik = matrix(rnorm(constants$I * (constants$K -1)), nrow = constants$I))

monitors = c("eta","alpha_i", "alpha_r", "beta_ir")


