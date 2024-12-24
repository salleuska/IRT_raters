rm(list=ls())
library(nimble)
setwd("C:/Users/39388/Dropbox/Il mio PC (LAPTOP-NO4UO9GH)/Desktop/Bocconi/Sally")

Data  <- read.csv("Data_2PL1PL.csv")


# ----
code2PL1PL <- nimbleCode({
  # Likelihood
  for(n in 1:tot){
    y[n,1] ~ dcat(ppi[n,1:K])
    
    for(k in 1:K){
      pic[n,k]          <-  exp(  lambda[II[n]] * (eta[PPi[n]] + beta[II[n]] + delta[k] + tau[RRi[n]] ) )   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
      pi[n,k]           <-  sum(pic[n,1:k])     
    }
    ppi[n,1:K]      <- pi[n,1:K]/sum(pi[n,1:K])
  }
  
  # Items parameters
  delta[1]  <- 0
  for(k in 2:K){
    delta[k]  ~ dnorm(0,  var = 25)   
  }
  
  for(i in 1:(I-1)) {
    beta[i]   ~ dnorm(0,  var = 25)  
    lambda[i] ~ dlnorm(meanlog = 0, sdlog = 3)
  }
  
  beta[I]   <-  -sum(beta[1:(I-1)])
  lambda[I] <-  -sum(log(lambda[1:(I-1)]))
  
  # Raters parameters
  
  for(i in 1:(R-1)) {
    tau[i] ~ dnorm(0,  var = 100) 
    }
  
  tau[R] <- -sum(tau[1:(R-1)])
  
  ## Individual effects
  
  for(p in 1:P) {
    eta[p] ~ dnorm(0, var = 1)  
  }
  
})

constants <- list(I = max(II), P = max(PPi),R = max(RRi), K=max(Data$y), 
                  tot = length(Data$y), II=II , PPi=PPi, RRi=RRi,
                  M=50)

data <- list(y = Data$y) #check

set.seed(2)



inits <- list(beta   = rnorm(constants$I, 0, 3),
              delta  = rnorm(constants$K, 0, 3),
              eta    = rnorm(constants$P, 0, 3),
              lambda = rep(1, constants$I),
              tau    = rnorm(constants$R, 0, 3)
              )

monitors = c("beta","delta", "lambda","tau", "eta")



#----

model2PL1PL   <- nimbleModel(code2PL1PL, constants, data, inits)

cmodel2PL1PL  <- compileNimble(model2PL1PL)

conf2PL1PL    <- configureMCMC(model2PL1PL, monitors = monitors)

modelMCMC     <- buildMCMC(conf2PL1PL)
cModelMCMC    <- compileNimble(modelMCMC, project = model2PL1PL)

system.time(samples <- runMCMC(cModelMCMC, niter=55000, nburnin = 5000, thin=10 ))

################################################################################

betaCols <- grep("beta", colnames(samples))
deltaCols <- grep("delta", colnames(samples))
lambdaCols <- grep("lambda", colnames(samples))

tauCols <- grep("tau", colnames(samples))

samplesSummary(samples[, c(betaCols)])
samplesSummary(samples[, c(deltaCols)])

samplesSummary(samples[, c(lambdaCols)])
samplesSummary(samples[, c(tauCols)])


#--- Trace
par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:4)
  ts.plot(samples[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ betaCols[i]])

par(mfrow = c(1, 2), cex = 1.1)
for(i in 2:K)
  ts.plot(samples[ , deltaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ deltaCols[i]])

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:4)
  ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ lambdaCols[i]])

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:R)
  ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ tauCols[i]])
################################################################################
calculateWAIC(samples, model2PL1PL)

################################################################################
