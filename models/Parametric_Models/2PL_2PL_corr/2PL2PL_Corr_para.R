rm(list=ls())
library(nimble)
setwd("C:/Users/39388/Dropbox/Il mio PC (LAPTOP-NO4UO9GH)/Desktop/Bocconi/Sally")

Data  <- read.csv("Data_2PL2PL.csv")

# ---- 
uppertri_mult_diag <- nimbleFunction(
  run = function(mat = double(2), vec = double(1)) {
    returnType(double(2))
    p <- length(vec)
    out <- matrix(nrow = p, ncol = p, init = FALSE)
    for(i in 1:p)
      out[ , i] <- mat[ , i] * vec[i]
    return(out)
  })

# ----
code2PL2PL <- nimbleCode({
  # Likelihood
  for(n in 1:tot){
    
    y[n] ~ dcat(ppi[n,1:K])
    
    pic[n,1]           <- 0
    pi[n,1]            <- 0
    
    for(k in 2:K){
      pic[n,k]          <-  exp(l_lambda[II[n]] + r_features[RRi[n],2]) * (eta[PPi[n]] - beta[II[n]] + delta[k-1] - r_features[RRi[n],1] )   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
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
  Ustar[1:2,1:2] ~ dlkj_corr_cholesky(1, 2)
  
  U[1:2,1:2]     <- uppertri_mult_diag(Ustar[1:2, 1:2], sds[1:2])
 
   for(i in 1:(R-1)) {
    r_features[i,1:2]   ~ dmnorm(mu_raters[1:2], cholesky = U[1:2, 1:2], prec_param = 0)
    
  }
  
  r_features[R,1]    <-  -sum(r_features[1:(R-1),1])
  r_features[R,2]    <-  -sum(r_features[1:(R-1),2])
  
  
  ## Individual effects
  
  for(p in 1:P) {
    eta[p] ~ dnorm(0, var = 1)  
  }
  
})


constants <- list(I = max(Data$II), P = max(Data$PPi),R = max(Data$RRi), K=max(Data$y), 
                  tot = length(Data$y), II=Data$II , PPi=Data$PPi, RRi=Data$RRi,
                  mu_raters=c(0,0))

data <- list(y = Data$y) #check

set.seed(2)


r_features = matrix(NA,constants$R,2)
r_features[,1]      = rnorm(constants$R, 0, 3)
r_features[,2]      = rep(0, constants$R)
inits <- list(beta                = rnorm(constants$I, 0, 3),
              delta               = rnorm(constants$K, 0, 3),
              eta                 = rnorm(constants$P, 0, 3),
              l_lambda            = rep(0, constants$I),
              r_features          = r_features, 
              Ustar               = matrix(1,2,2),
              U                   = matrix(1,2,2),
              sds                 = rep(1,2)
)

monitors = c("beta","delta", "l_lambda","r_features", "Ustar", "eta")



#----

model2PL2PL         <- nimbleModel(code2PL2PL, constants, data, inits)

cmodel2PL2PL        <- compileNimble(model2PL2PL)

conf2PL2PL          <- configureMCMC(model2PL2PL, monitors = monitors)

modelMCMC           <- buildMCMC(conf2PL2PL)
cModelMCMC          <- compileNimble(modelMCMC, project = model2PL2PL)

system.time(samples <- runMCMC(cModelMCMC, niter=55000, nburnin = 5000, thin=10 ))

################################################################################

betaCols <- grep("beta", colnames(samples))
deltaCols <- grep("delta", colnames(samples))
lambdaCols <- grep("lambda", colnames(samples))
tauCols <- grep("tau", colnames(samples))
phiCols <- grep("phi", colnames(samples))

samplesSummary(samples[, c(betaCols)])
samplesSummary(samples[, c(deltaCols)])
samplesSummary(samples[, c(lambdaCols)])

samplesSummary(samples[, c(tauCols)])
samplesSummary(samples[, c(phiCols)])


#--- Trace
K = max(Data$y)
R = max(Data$RRi)
I = max(Data$II)

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:I)
  ts.plot(samples[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ betaCols[i]])

par(mfrow = c(1, 2), cex = 1.1)
for(i in 1:(K-1))
  ts.plot(samples[ , deltaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ deltaCols[i]])

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:I)
  ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ lambdaCols[i]])

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:R)
  ts.plot(samples[ , tauCols[i]], xlab = 'iteration', ylab = colnames(samples)[ tauCols[i]])

par(mfrow = c(2, 2), cex = 1.1)
for(i in 1:R)
  ts.plot(samples[ , phiCols[i]], xlab = 'iteration', ylab = colnames(samples)[ phiCols[i]])

################################################################################
calculateWAIC(samples, model2PL2PL)

################################################################################
