#!/bin/bash
#############################################
## Bash script to replicate analysis in the paper
#############################################
## create all directories needed for outpu 

mkdir -p figures
mkdir -p output


#############################################
## 0) Simulate data
#############################################

## The script called below simulate data 

#############################################
## 1) RUN ALL MODELS 
#############################################
## Run all models for simulated data 
#############################################
## UNIMODAL SIMULATION 
#############################################

## Parametric models

Rscript 1_runNimbleModel.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulated/data_2PL2PL.rds \
--niter=1000 \
--nburnin=500 \
--nthin=1 


#############################################
## 2) get and process MCMC samples
#############################################



#############################################
## 4) Make plots
############################################


