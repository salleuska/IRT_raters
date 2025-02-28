#!/bin/bash
#############################################
## Bash script to replicate analysis in the paper
#############################################
## create all directories needed for output if not present 

mkdir -p figures
mkdir -p output



#############################################
## 1) RUN MODELS 
#############################################
## Run models

Rscript 1_runNimbleModel.R  \
--model=models/semiparametric/semi_2PL_2PL.R \
--data=data/osce/OSCE_Long.rds \
--niter=1000 \
--nburnin=500 \
--nthin=1 


#############################################
## 2) get and process MCMC samples
#############################################



#############################################
## 4) Make plots
############################################


