##-----------------------------------------#
## last update: 
## R version 
## nimble version
##-----------------------------------------#
## Same usage as main.sh but to run models from R
##-----------------------------------------#

args <- list()
args$model 		= "models/parametric/para_2PL_2PL.R" 
args$data 		= "data/simulated/data_2PL2PL.rds"
args$niter 		= 1000 
args$nburnin 	= 500 
args$nthin 		= 1

source("1_runNimbleModel.R")
