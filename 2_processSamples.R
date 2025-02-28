##-----------------------------------------#
## This script post-process the posterior samples 
## meaning? 
##-----------------------------------------#
library(here)
source("functions/estimateDPdensity.R")
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --resFileName
##-----------------------------------------#
## TMP for testing
args <- list(resFileName = "output/OSCE_Long/semi/semi_2PL_2PL.rds")
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
## TO DO: likely some rescaling for comparison - parametric and semi
##-------------------------------------------------------##
samples <- resObj$samples

##-------------------------------------------------------##
## !!!!!!!!!! START FROM HERE
##-------------------------------------------------------##
betaCols <- grep("beta", colnames(samples))
deltaCols <- grep("delta", colnames(samples))
lambdaCols <- grep("lambda", colnames(samples))

hist(samples[, grep("^eta", colnames(samples))])


tauCols <- grep("tau", colnames(samples))
nimble::samplesSummary(samples[, c(betaCols)])
nimble::samplesSummary(samples[, c(deltaCols)])

nimble::samplesSummary(samples[, c(lambdaCols)])
nimble::samplesSummary(samples[, c(tauCols)])

out <- estimateDPdensity(samples, nIndividuals = 30)
res <- data.frame(grid = out$grid, mean = apply(out$densitySamples, 2, mean))

library(ggplot2)

ggplot(res, aes(x=grid, y = mean)) + geom_line()
