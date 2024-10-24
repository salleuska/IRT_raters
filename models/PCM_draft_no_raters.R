########################################################
## PCM  wide format for the data -- no ratesr
## Y[i,j]  indicating student i , item j
## N - number of students
## nItems - number of items
## K - number of possible scores; Note that scores are defined from 1 to K, not from 0!
########################################################
library(nimble)
## nimble model

modelCode =  nimbleCode({

	for(i in 1:N) 
	{
		for(j in 1:nItems)
		{
	         ## data generative mechanism
		 Y[i,j] ~ dcat(prob[i,j,1:K])
		}
	 ## students' ability
	 theta[i] ~ dnorm(0, 1)  
	}


	for(i in 1:N)
	{
         for(j in 1:nItems)
	 {
          for(k in 1:K)
	  {
           eta[i,j, k] <- alpha[j] * (theta[i] - beta[j, k])
           psum[i,j, k] <- sum(eta[i, j, 1:k]) ## cumulative sums
           exp.psum[i, j, k] <- exp(psum[i, j, k])
	  }
	 }
	}

	## normalize probabilities 
	for(i in 1:N)
	{
         for(j in 1:nItems)
	 {

          for(k in 1:K)
	  {
           prob[i, j, k] <- exp.psum[i, j, k] /(sum(exp.psum[i, j, 1:K]))
	  }
	 }
	}
	## prior specification
        for(j in 1:nItems)
	{
     	 log(alpha[j]) ~ dnorm(0.5, 0.5)
          beta[j, 1] <- 0
          for (k in 2:K)
	  {
           beta[j,  k] ~ dnorm(0, var = 2) 
          }
         }
})
