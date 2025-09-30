##-----------------------------------------#
## This script process posterior samples
## and estimate quantites for inference
##-----------------------------------------#
library(here)
source("functions/estimateDPdensity.R")
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --resFileName
## --simFileName
##-----------------------------------------#
## TMP for testing
args <- list(resFileName = "output/data_bimodal_noAR/semi/semi_2PL_2PL.rds", 
			 simFileName = "data/simulated/data_bimodal_noAR.rds")
##-----------------------------------------#
# if(is.null(args$outDir)) outDir <- "output/posterior_samples_elaborated/" else dir <- args$outDir

listLength <- length(strsplit(args$resFileName, "\\/|.rds")[[1]])
data <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength -2]
fileName <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength]

## modelType - parametric or semi
modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]
# ## 
# modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]

## read objects
resObj <- readRDS(args$resFileName)

simData <- readRDS(args$simFileName)
##-------------------------------------------------------##
## TO DO: likely some rescaling for comparison - parametric and semiparametric
## SP: Need to be sure on the constraints?
##-------------------------------------------------------##
samples <- resObj$samples
colnames(samples)

etaCols      <- grep("^eta(\\[|$)",    colnames(samples))
betaCols     <- grep("^beta(\\[|$)",   colnames(samples))
deltaCols    <- grep("^delta(\\[|$)",  colnames(samples))   # thresholds (K-1)
tauCols      <- grep("^tau(\\[|$)",    colnames(samples))   # rater biases
llambaCols   <- grep("^l_lambda(\\[|$)", colnames(samples)) # item log-slopes
lphiCols     <- grep("^l_phi(\\[|$)",    colnames(samples)) # rater log-slopes

P <- length(etaCols)
##-------------------------------------------------------##
## !!!!!!!!!! START FROM HERE
##-------------------------------------------------------##

betaMeans <- colMeans(samples[, betaCols])

plot(simData$beta, betaMeans)
abline(0,1)

deltaMeans <- colMeans(samples[, deltaCols])

plot(simData$delta, deltaMeans)
abline(0,1)

## check constraints
out <- estimateDPdensity(samples, nIndividuals = P, 
	grid = seq(-8, 8, length = 200))
res <- data.frame(grid = out$grid, mean = apply(out$densitySamples, 2, mean))

library(ggplot2)

ggplot(res, aes(x = grid, y = mean)) +
  geom_line(aes(color = "Posterior (mean)"), linewidth = 1) +
  geom_density(
    data = data.frame(eta = simData$eta),       # <- or simData$eta
    aes(x = eta, y = after_stat(density), color = "Simulated (kernel)"),
    inherit.aes = FALSE, linewidth = 0.8, alpha = 0.7
  ) +
  scale_color_manual(values = c(
    "Posterior (mean)"      = "black",
    "Simulated (kernel)"    = "steelblue"
  )) +
  labs(x = "ability", y = "density", color = "") +
  theme_minimal()


