##-----------------------------------------#
## Semiparametric - [item] 2PL [rater] 2PL model 
##-----------------------------------------#
modelCode <- nimbleCode({
  # Likelihood
  for(n in 1:tot){
    
    y[n,1] ~ dcat(ppi[n,1:K])
    
    pic[n,1]           <- 0
    pi[n,1]            <- 0
    
    for(k in 2:K){
      pic[n,k]          <-  exp(l_lambda[II[n]] + l_phi[RRi[n]]) * (eta[PPi[n]] - beta[II[n]] + delta[k-1] - tau[RRi[n]] )   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
      pi[n,k]           <-  sum(pic[n,1:k])     
    }
    ppi[n,1:K]          <- exp(pi[n,1:K])/sum(exp(pi[n,1:K]))
  }
  
  # Items parameters
  for(k in 1:(K-1)){
    delta[k]  ~ dnorm(0,  var = 25)   
  }
  
  for(i in 1:(I-1)) {
    beta[i]     ~ dnorm(0,  var = 25)  
    l_lambda[i] ~ dnorm(0,  var = 3)
  }
  
  beta[I]   <-  -sum(beta[1:(I-1)])
  
  l_lambda[I] <-  -sum(l_lambda[1:(I-1)])
  

  # Raters parameters
  
  for(i in 1:(R-1)) {
    tau[i]   ~ dnorm(0,  var = 25) 
    l_phi[i] ~ dnorm(0,  var = 3)
  }
  
  tau[R]    <-  -sum(tau[1:(R-1)])
  l_phi[R]  <-  -sum(l_phi[1:(R-1)])
  

  ## Individual effects
  
  ## CRP for clustering individual effects
  zi[1:P] ~ dCRP(alpha, size = P)
  alpha ~ dgamma(a, b)  
  ## Mixture component parameter drawn from the base measure
  for(j in 1:P) {
    eta[j] ~ dnorm(mu[j], var = s2[j])  
    mu[j] <- muTilde[zi[j]]                 
    s2[j] <- s2Tilde[zi[j]]   
  }
  
  for(m in 1:M) {
    muTilde[m] ~ dnorm(0, var = s2_mu)
    s2Tilde[m] ~ dinvgamma(nu1, nu2)
  }
  
  
})


constants <- list(I = max(Data$II), P = max(Data$PPi),R = max(Data$RRi), K=max(Data$y), 
                  tot = length(Data$y), II=Data$II , PPi=Data$PPi, RRi=Data$RRi,
                  M=50)

data <- list(y = Data$y) #check



inits <- list(beta     = rnorm(constants$I, 0, 3),
              delta    = rnorm(constants$K-1, 0, 3),
              eta      = rnorm(constants$P, 0, 3),
              l_lambda = rep(0, constants$I),
              tau      = rnorm(constants$R, 0, 3),
              l_phi    = rep(0, constants$R),
              zi     = sample(1:constants$M,constants$P,replace = TRUE),
              alpha  = 1,
              muTilde= rep(0,constants$M),
              s2Tilde= rep(1,constants$M),
              mu     = rep(0,constants$P),
              s2     = rep(1,constants$P),
              nu1    = 2.01, 
              nu2    = 1.01,
              s2_mu  = 2,
              a      = 1,
              b      = 3
)

monitors = c("beta","delta", "l_lambda","tau","l_phi", "eta","zi", "muTilde", "s2Tilde", "alpha")


# #----

# model2PL2PL         <- nimbleModel(code2PL2PL, constants, data, inits)

# cmodel2PL2PL        <- compileNimble(model2PL2PL)

# conf2PL2PL          <- configureMCMC(model2PL2PL, monitors = monitors)

# modelMCMC           <- buildMCMC(conf2PL2PL)
# cModelMCMC          <- compileNimble(modelMCMC, project = model2PL2PL)

# system.time(samples <- runMCMC(cModelMCMC, niter=55000, nburnin = 5000, thin=10 ))

# ################################################################################

# betaCols <- grep("beta", colnames(samples))
# deltaCols <- grep("delta", colnames(samples))
# lambdaCols <- grep("lambda", colnames(samples))
# tauCols <- grep("tau", colnames(samples))
# phiCols <- grep("phi", colnames(samples))

# samplesSummary(samples[, c(betaCols)])
# samplesSummary(samples[, c(deltaCols)])
# samplesSummary(samples[, c(lambdaCols)])

# samplesSummary(samples[, c(tauCols)])
# samplesSummary(samples[, c(phiCols)])


# #--- Trace
# K = max(Data$y)
# R = max(Data$RRi)
# I = max(Data$II)

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:I)
#   ts.plot(samples[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ betaCols[i]])

# par(mfrow = c(1, 2), cex = 1.1)
# for(i in 1:(K-1))
#   ts.plot(samples[ , deltaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ deltaCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:I)
#   ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ lambdaCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:R)
#   ts.plot(samples[ , tauCols[i]], xlab = 'iteration', ylab = colnames(samples)[ tauCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:R)
#   ts.plot(samples[ , phiCols[i]], xlab = 'iteration', ylab = colnames(samples)[ phiCols[i]])

# ################################################################################
# calculateWAIC(samples, model2PL2PL)

# ################################################################################
