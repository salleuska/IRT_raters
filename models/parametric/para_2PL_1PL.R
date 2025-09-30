##-----------------------------------------#
## [item] 2PL [rater] 1PL model 
## - Adjacent-categories (GPCM-style)
## - Item slope:  l_lambda[i]
## - Item diff:   beta[i]
## - Rater bias:  tau[r]
## - Thresholds:  delta[k]  
## - Subjects' eta  N(0,1)
##-----------------------------------------#
modelCode <- nimbleCode({

  ##---------------- Likelihood ----------------##
  for(n in 1:tot){

    # Observation n belongs to subject j, item i, rater r
    # indices provided as data: PPi[n] (j), II[n] (i), RRi[n] (r)

    y[n] ~ dcat(ppi[n, 1:K])

    # Baseline category (k = 1) has zero adjacent/cumulative logit
    pic[n, 1] <- 0
    pi[n, 1]  <- 0

    # Common log-slope term: exp(l_lambda[item] + l_phi[rater])
    slope[n]<- exp( l_lambda[ II[n] ])

    for(k in 2:K){
      # Adjacent-category logit piece for category k:
      # slope * ( eta - beta - delta - tau )
      # NOTE: delta (thresholds) ENTERS WITH A MINUS SIGN.
      pic[n, k] <- slope[n] * ( eta[ PPi[n] ] 
                               -  beta[ II[n] ] 
                               -  delta[k - 1] 
                               -  tau[ RRi[n] ] )
      # Cumulative logit up to k
      pi[n, k]  <- sum( pic[n, 1:k] )
    }

    # Softmax over cumulative logits to get probabilities
    ppi[n, 1:K] <- exp( pi[n, 1:K] ) / sum( exp( pi[n, 1:K] ) )
  }

  ##---------------- Item parameters ----------------##
  # Threshold steps (shared across raters)
  for(k in 1:(K-1)){
    delta[k] ~ dnorm(0, var = 3)
  }

  # Item locations and log-slopes with sum-to-zero constraints
  for(i in 1:(I-1)) {
    beta[i]     ~ dnorm(0, var = 3)
    l_lambda[i] ~ dnorm(0, var = 3)
  }
  beta[I]      <- -sum(beta[1:(I-1)])
  l_lambda[I]  <- -sum(l_lambda[1:(I-1)])

  ##---------------- Rater parameters ----------------##
  # Rater biases (tau) and log-slopes (l_phi) with sum-to-zero constraints
  for(r in 1:(R-1)) {
    tau[r]   ~ dnorm(0, var = 3)
  }
  tau[R]   <- -sum(tau[1:(R-1)])

  ##---------------- Subject latent traits ----------------##
  for(p in 1:P) {
    eta[p] ~ dnorm(0, var = 1)  
  }
})



constants <- list(I = max(Data$II), P = max(Data$PPi),R = max(Data$RRi), K=max(Data$y), 
                  tot = length(Data$y), II=Data$II , PPi=Data$PPi, RRi=Data$RRi)
data <- list(y = Data$y) #check

inits <- list(beta   = rnorm(constants$I, 0, 3),
              delta  = rnorm(constants$K, 0, 3),
              eta    = rnorm(constants$P, 0, 3),
              l_lambda = rep(0, constants$I),
              tau    = rnorm(constants$R, 0, 3)
              )

monitors = c("beta","delta", "l_lambda","tau", "eta")




# model2PL1PL         <- nimbleModel(code2PL1PL, constants, data, inits)

# cmodel2PL1PL        <- compileNimble(model2PL1PL)

# conf2PL1PL          <- configureMCMC(model2PL1PL, monitors = monitors)

# modelMCMC           <- buildMCMC(conf2PL1PL)
# cModelMCMC          <- compileNimble(modelMCMC, project = model2PL1PL)

# system.time(samples <- runMCMC(cModelMCMC, niter=55000, nburnin = 5000, thin=10 ))

# ################################################################################

# betaCols <- grep("beta", colnames(samples))
# deltaCols <- grep("delta", colnames(samples))
# lambdaCols <- grep("lambda", colnames(samples))

# tauCols <- grep("tau", colnames(samples))

# samplesSummary(samples[, c(betaCols)])
# samplesSummary(samples[, c(deltaCols)])

# samplesSummary(samples[, c(lambdaCols)])
# samplesSummary(samples[, c(tauCols)])


# #--- Trace
# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:4)
#   ts.plot(samples[ , betaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ betaCols[i]])

# par(mfrow = c(1, 2), cex = 1.1)
# for(i in 2:max(Data$y))
#   ts.plot(samples[ , deltaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ deltaCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:4)
#   ts.plot(samples[ , lambdaCols[i]], xlab = 'iteration', ylab = colnames(samples)[ lambdaCols[i]])

# par(mfrow = c(2, 2), cex = 1.1)
# for(i in 1:max(Data$RRi))
#   ts.plot(samples[ , tauCols[i]], xlab = 'iteration', ylab = colnames(samples)[ tauCols[i]])
# ################################################################################
# calculateWAIC(samples, model2PL1PL)

# ################################################################################
