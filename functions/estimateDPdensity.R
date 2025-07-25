## The following function calculates the DP predictive on a grid of points
## Assume a DP mixture of normals
## The predictive is calculated for each MCMC samples
##
## @ samples. MCMC samples matrix with R samples and P parameters
## @ paraNames. Named list where characteres refer to the parameters name in the sample matrix
##              "alpha" refers to the concentration parameter, "muTilde" to the cluster means, 
##              "s2Tilde" to the cluster variances, "zi" to the cluster memberships
## @ nIndividuals. Number of individuals  
## @ grid. Grid of points where the predictive density will be calculated
##
## The function returns a matrix with R rows and the predictive density for each point on the grid


estimateDPdensity <- function(samples, 
                              paramNames = list(alpha = "alpha", 
                                                muTilde = "muTilde", 
                                                s2Tilde = "s2Tilde", 
                                                zi = "zi"), 
                              nIndividuals, 
                              grid = seq(-10, 10, len = 200)
                              ) {

  # posterior samples of the concentration parameter
  alphaSamples <- samples[ , paramNames$alpha] 
  # posterior samples of the cluster means
  muTildeSamples <- samples[ , grep(paramNames$muTilde, colnames(samples))] 
  # posterior samples of the cluster variances
  s2TildeSamples <- samples[ , grep(paramNames$s2Tilde, colnames(samples))] 
  # posterior samples of the cluster memberships
  ziSamples <- samples [ , grep(paramNames$zi, colnames(samples))] 


  densitySamples <- matrix(0, ncol = length(grid), nrow = nrow(samples)) 

  ## Compute the predictive density for every point on the grid

  for(i in 1:nrow(samples)){
    k <- unique(ziSamples[i, ])
    kNew <- max(k) + 1
    mk <- c()
    li <- 1
    for(l in 1:length(k)) {
      mk[li] <- sum(ziSamples[i, ] == k[li])
      li <- li + 1
    }
    alpha <- alphaSamples[i]
    
    muK <-  muTildeSamples[i, k]
    s2K <-  s2TildeSamples[i, k]
    muKnew <-  muTildeSamples[i, kNew]
    s2Knew <-  s2TildeSamples[i, kNew]
    
    densitySamples[i, ] <- sapply(grid, 
                  function(x)(sum(mk * dnorm(x, muK, sqrt(s2K))) +
                  alpha * dnorm(x, muKnew, sqrt(s2Knew)) )/(alpha+nIndividuals))
  }

  ## Return values of the function
  out <- list(grid = grid, densitySamples= densitySamples)
  return(out)
  
}