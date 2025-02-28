##-----------------------------------------#
## Semiparametric - AR - [item] 2PL [rater] 2PL model 
## Correlated raters' features
##-----------------------------------------#

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
modelCode <- nimbleCode({
  # Likelihood
  for(n in 1:tot){
    
    y[n] ~ dcat(ppi[n,1:K])
    
    pic[n,1]           <- 0
    pi[n,1]            <- 0
    
    for(k in 2:K){
      pic[n,k]          <-  exp(l_lambda[II[n]] + r_features[RRi[n],2]) * (eta[PPi[n]+1] - beta[II[n]] + delta[k-1]
                                                                           - r_features[RRi[n],1] - r_features[RRi[n],3]*eta[ARi[n]+1])   # Adjacent categories logits (Agresti, 2013, pp.309-310; Masters, 1982, p.158)
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
  Ustar[1:3,1:3] ~ dlkj_corr_cholesky(1.3, 3)
  
  U[1:3,1:3]     <- uppertri_mult_diag(Ustar[1:3, 1:3], sds[1:3])
  
  for(i in 1:(R-1)) {
    r_features[i,1:3]   ~ dmnorm(mu_raters[1:3], cholesky = U[1:3, 1:3], prec_param = 0)
    
  }
  r_features[R,3]    ~ dnorm(0, var = 3) 
  
  r_features[R,1]    <-  -sum(r_features[1:(R-1),1])
  r_features[R,2]    <-  -sum(r_features[1:(R-1),2])
  
  
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
                  mu_raters=c(0,0,0), ARi=Data$ARi, M=50)

data <- list(y = Data$y) #check




r_features          = matrix(NA,constants$R,3)
r_features[,1]      = rnorm(constants$R, 0, 3)
r_features[,2]      = rep(0, constants$R)
r_features[,3]      = rnorm(constants$R, 0, 0.2)

inits <- list(beta                = rnorm(constants$I, 0, 3),
              delta               = rnorm(constants$K-1, 0, 3),
              eta                 = c(0, rnorm(constants$P, 0, 3)),
              l_lambda            = rep(0, constants$I),
              r_features          = r_features, 
              Ustar               = matrix(1,3,3),
              U                   = matrix(1,3,3),
              sds                 = rep(1,3),
#              zi     = sample(1:constants$M,constants$P,replace = TRUE),
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

monitors = c("beta","delta", "l_lambda","r_features", "U", "eta","zi", "muTilde", "s2Tilde", "alpha")


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
# tauCols <- grep("r_features", colnames(samples))
# UstarCols <- grep("U", colnames(samples))
# etaCols   <- grep("eta", colnames(samples))

# samplesSummary(samples[, c(betaCols)])
# samplesSummary(samples[, c(deltaCols)])
# samplesSummary(samples[, c(lambdaCols)])


# samplesSummary(samples[, c(tauCols)])


# #--- Trace
# K = max(Data$y)
# R = max(Data$RRi)
# I = max(Data$II)

# # Item 
# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:I)
#   ts.plot(samples[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ betaCols[i]])

# par(mfrow = c(1, 2), cex = 1.1)
# for(i in 1:(K-1))
#   ts.plot(samples[ , deltaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ deltaCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:I)
#   ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ lambdaCols[i]])

# # Rater
# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:(2*R))
#   ts.plot(samples[ , tauCols[i]], xlab = 'iteration', ylab = colnames(samples)[ tauCols[i]])


# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:9)
#   ts.plot(samples[ , UstarCols[i]], xlab = 'iteration', ylab = colnames(samples)[ UstarCols[i]])

# ################################################################################
# calculateWAIC(samples, model2PL2PL)

# ################################################################################
# # ACCURACY CHECK
# par(mfrow = c(1, 1), cex = 1.1)
# # --- ETA
# eta_hat <- apply(samples[ , 20:119],2,mean) #check the indeces
# plot(Data$eta,eta_hat,pch=19,ylim=c(-10,10),xlim = c(-10,10),main="Eta - Accuracy",
#      xlab="TRUE ETA", ylab="ESTIMATES")
# abline(c(0,1),xpd=FALSE,col="red")

# # --- BETA 
# beta_hat <- apply(samples[ , betaCols],2,mean) #check the indeces
# plot(Data$beta,beta_hat,pch=19,ylim=c(-0.8,0.8),xlim = c(-0.8,0.8),main="Beta - Accuracy",
#      xlab="TRUE BETA", ylab="ESTIMATES")
# abline(c(0,1),xpd=FALSE,col="red")
# # --- DELTA 
# delta_hat <- apply(samples[ , deltaCols[1:3]],2,mean) #check the indeces
# plot(Data$delta,delta_hat,pch=19,ylim=c(-10,10),xlim = c(-10,10),main="Delta - Accuracy",
#      xlab="TRUE DELTA", ylab="ESTIMATES")
# abline(c(0,1),xpd=FALSE,col="red")
# # --- l_LAMBDA 
# l_lambda_hat <- apply(samples[ , lambdaCols],2,mean) #check the indeces
# plot(Data$lambda,l_lambda_hat,pch=19,ylim=c(-2,2),xlim = c(-2,2),main="Delta - Accuracy",
#      xlab="TRUE l_LAMBDA", ylab="ESTIMATES")
# abline(c(0,1),xpd=FALSE,col="red")



