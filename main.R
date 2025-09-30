##-----------------------------------------#
## last update: 
## R version 
## nimble version
##-----------------------------------------#
## Same usage as main.sh but to run models from R
##-----------------------------------------#

## Here exaple
args <- list(
  model   = "models/parametric/para_2PL_2PL.R",
  data    = "data/simulated/data_2PL2PL.rds",
  niter   = 1000,
  nburnin = 500,
  nthin   = 1
)

source("1_runNimbleModel.R")