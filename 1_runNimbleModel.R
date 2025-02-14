##-----------------------------------------#
## last update: 
## R version 
## nimble version
##-----------------------------------------#
## This scripts runs a nimble model
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)

## Accepted arguments
## model 				path to R script with model code
## data 				path to data file
## niter        number of iterations for MCMC
## nburnin 			number of burnin iterations for MCMC
## nthin 				thinning (deafault is 1)
## dirResults	  (optional) path directory for results
##-----------------------------------------##
## TMP - for testing the script
##-----------------------------------------##
args <- list()
args$model 		= "models/parametric/para_2PL_2PL.R" 
args$data 		= "data/simulated/data_2PL2PL.rds"
args$niter 		= 1000 
args$nburnin 	= 500 
args$nthin 		= 1

##-----------------------------------------##
## Load libraries and functions
##-----------------------------------------##
library(nimble)
library(here)

##-----------------------------------------##
## Set variables from args list
##-----------------------------------------##
calcWAIC <- TRUE

## results directory
if(is.null(args$dirResults)) dir <- "output" else dir <- args$dirResults

## filename used for output
filename <- unlist(strsplit(basename(args$model), "[\\.]"))[1]

## MCMC settings
MCMCcontrol 		<- list()
MCMCcontrol$niter 	<- as.numeric(args$niter)
MCMCcontrol$nburnin <- as.numeric(args$nburnin)

## set seed based on slurm task id if available
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if(task_id == "") seed <- 1 else seed <- 1 + as.numeric(task_id)

MCMCcontrol$seed <- seed

cat("##--------------------------------##\n")
cat("Model ", filename, "\n")
cat("Data ", args$data, "\n")
cat("##--------------------------------##\n")

##---------------------------------------##
## Read data 
## Note: constants and inits are handles within each model script

Data 	<- readRDS(args$data)

##---------------------------------------##
## Source model code  
source(args$model)

##---------------------------------------##
## START HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!
##---------------------------------------##
## init random effects using standardized raw score

# data.frame(Data$y, Data$PPi)
# by(Data$y, as.factor(Data$PPi),function(x) sum(x$y), simplify = T)

# scores 		<- as.vector(by(Data$y, as.factor(Data$PPi),function(x) sum(x$y), simplify = T))
# Sscores 	<- (scores - mean(scores))/sd(scores)
# inits$eta 	<- Sscores



## BNP inits for data application
if(grepl("bnp", args$model)) {
	if(grepl("health", args$data)) {
		inits$zi 	<- kmeans(Sscores, 3)$cluster
		inits$a 	<- 1        # gamma prior parameters
		inits$b 	<- 3        # gamma prior parameters
	}

	## hyperamenters for simulated data
	if(grepl("simulation", args$data)) {
		inits$zi 	<- kmeans(Sscores, 4)$cluster
		inits$a 	<- 2
		inits$b 	<- 4
	}
	## hyperamenters for timss data
	if(grepl("timss", args$data)) {
		inits$zi 	<- kmeans(Sscores, 3)$cluster
		constants$M <- 30   	# number of clusters
		inits$a 	<- 1        # gamma prior parameters
		inits$b 	<- 3        # gamma prior parameters
	}
}

##---------------------------------------------------##
## Create model and MCMC configuration
##---------------------------------------------------##


model <- nimbleModel(code 			= modelCode,
										 data 			= data,  
										 constants	= constants,
										 inits 			= inits, 
										 calculate 	= FALSE)


## update monitors
monitors <- c(monitors, "myLogProbAll", "myLogProbSome", "myLogLik")

## Flag for WAIC
if(calcWAIC) {
## conditional WAIC - grouped students
  if(grepl("timss", args$data)) {
		indList <- split(seq_along(alldata$id), alldata$id)
		groups <- sapply(indList, function(x)  paste0('y[', x, ']'))  
  }
  else {
		groups <- paste0('y[', 1:constants$N, ', 1:', constants$I, ']')
  }


	mcmcConf <- configureMCMC(model, monitors = monitors, 
		enableWAIC = TRUE, 
		waicControl = list(dataGroups = groups))

} else {
	mcmcConf <- configureMCMC(model, monitors = monitors)

}

## Add samplers to monitor log posterior probability and likelihood
nodeList = c("beta", "lambda", "eta")
if("gamma" %in% monitors) nodeList[1] <- "gamma"

mcmcConf$removeSampler("myLogLik")
mcmcConf$addSampler("myLogLik", type =  "logProb_summer", 
  control = list(nodeList = c("y")))

mcmcConf$removeSampler("myLogProbAll")
mcmcConf$addSampler("myLogProbAll", type =  "logProb_summer")

mcmcConf$removeSampler("myLogProbSome")
mcmcConf$addSampler("myLogProbSome", type =  "logProb_summer", 
  control = list(nodeList = nodeList))

mcmcConf


## sampler configuration changes according to mode
if(args$mode == "centered" ) {

	if(("gamma" %in% monitors) & grepl("constrainedAbilities|unconstrained", filename)){ 

		  mcmcConf$removeSamplers("log_lambda")
		  # mcmcConf$removeSamplers("gamma")
  
 		  mcmcConf$addSampler(type = 'centered',
       			target = c('log_lambda', 'gamma'),
       			control = list(nodesToCenter = 'eta', scale = 0.1, adaptive = TRUE))

	} else {
	  q(save = 'no')
	}
}


mcmcConf$addMonitors2("eta")
mcmcConf$setThin2(MCMCcontrol$thin2)
mcmc <- buildMCMC(mcmcConf)	

## Add thinning for all variables when running on data
if(grepl("timss", args$data) | grepl("health", args$data)) {
	mcmcConf$setThin(MCMCcontrol$thin2)
}



##---------------------------------------------------##
## Compile model & MCMC 
##---------------------------------------------------##
compilationTime <- system.time({
    Cmodel <- try(compileNimble(model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble model.")
    }
    Cmcmc <- try(compileNimble(mcmc, project = model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble MCMC.")
    }
})


##---------------------------------------------------##
## Run MCMC 
##---------------------------------------------------##
runningTime <- system.time({try(
	res <- runMCMC(Cmcmc, 
				   niter 	= MCMCcontrol$niter, 
				   nburnin  = MCMCcontrol$nburnin,
				   setSeed  = seed))
	if(inherits(res, 'try-error')) {
  		warning(paste0("There was a problem running nimble MCMC."))
  }
})

##---------------------------------------------------##
## Save results, times, settings
##---------------------------------------------------##
results <- list(samples  = res,
				compilationTime  = compilationTime,
				samplingTime     = runningTime*(1 - MCMCcontrol$niter/MCMCcontrol$nburnin),
				runningTime      = runningTime,
				MCMCcontrol      = MCMCcontrol)

if(calcWAIC) results$modelWAIC <- Cmcmc$getWAIC()$WAIC


##---------------------------------------------------##
## directory for output
##---------------------------------------------------##
modelType       <- unlist(strsplit(basename(args$model), "[\\_\\.]"))[1]
dataName        <- unlist(strsplit(basename(args$data), "[.]"))[1]

outDir <- paste0(dir, "/", dataName, "/", modelType, "/")

dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)


if(grepl("centered", args$mode)) {
	filenameOutput <- paste0(outDir, filename, "_centered.rds")
} else {
	filenameOutput <- paste0(outDir, filename, ".rds")
}

saveRDS(results, file = filenameOutput )
