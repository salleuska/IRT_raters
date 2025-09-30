##-----------------------------------------#
## This script process posterior samples
## and estimate quantites for inference
##-----------------------------------------#
library(here)
source("functions/estimateDPdensity.R")
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --resFileName
##-----------------------------------------#
## TMP for testing
## args <- list(resFileName = "output/OSCE_Long/para/para_uto.rds")
args <- list(resFileName = "output/data_bimodal_noAR/semi/semi_2PL_2PL.rds")
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

##-------------------------------------------------------##
## !!!!!!!!!! START FROM HERE
##-------------------------------------------------------##



hist(samples[, grep("^eta", colnames(samples))], 
	breaks = 100)


tauCols <- grep("tau", colnames(samples))
nimble::samplesSummary(samples[, c(betaCols)])
nimble::samplesSummary(samples[, c(deltaCols)])

nimble::samplesSummary(samples[, c(lambdaCols)])
nimble::samplesSummary(samples[, c(tauCols)])

out <- estimateDPdensity(samples, nIndividuals = 30, 
	grid = seq(-8, 8, length = 200))
res <- data.frame(grid = out$grid, mean = apply(out$densitySamples, 2, mean))

library(ggplot2)

ggplot(res, aes(x=grid, y = mean)) + geom_line()
