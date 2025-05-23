##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
## These script containts functions to compute the expectation and variance
## of the a priori number of clusters via Monte Carlo when using a prior 
## distribution for the DP concentration parameter alpha
#################################################################

expectedNumberOfClusters <- function(N, alpha){
	sum(alpha/(alpha + seq(1:N) - 1))
}

varianceNumberOfClusters <- function(N, alpha){
	num <- alpha*(seq(1:N) -1)
	den <- (alpha - 1 + seq(1:N))^2
	sum(num/den)
}

#########################
## Example 
#########################
## number of observations in the data - simulation
N <- 1000
## parameters for the gamma distribution
a <- 2
b <- 4

## number of Monte carlo replicates
R <- 10^5

## generate alpha ~ Ga(a, b)
set.seed(1)
alphaVals <- rgamma(R, a, b)

kVals <- sapply(alphaVals, function(x) expectedNumberOfClusters(N, x))
kValsVar <- sapply(alphaVals, function(x) varianceNumberOfClusters(N, x))

## Estimate of the a priori expected number of clusters
mean(kVals)

## Estimate of the a priori variance for thenumber of clusters
## uses law of total variance -  Mean(var(x | alpha)) + var(mean(x|alpha))
mean(kValsVar) + var(kVals)

#########################

## number of observations in the data - OCSE
N <- 30
## parameters for the gamma distribution
a <- 1
b <- 3

## number of Monte carlo replicates
R <- 10^5

## generate alpha ~ Ga(a, b)
set.seed(1)
alphaVals <- rgamma(R, a, b)

kVals <- sapply(alphaVals, function(x) expectedNumberOfClusters(N, x))
kValsVar <- sapply(alphaVals, function(x) varianceNumberOfClusters(N, x))

## Estimate of the a priori expected number of clusters
mean(kVals)

## Estimate of the a priori variance for thenumber of clusters
## uses law of total variance -  Mean(var(x | alpha)) + var(mean(x|alpha))
mean(kValsVar) + var(kVals)
