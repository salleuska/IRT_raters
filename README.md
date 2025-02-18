[Draft]

Supplementary code for the paper "TBD" (2025+)

This repository contains code to estimate the semiparametric IRT raters model. 

### Requirements

Make sure to install the following R packages and relative dependencies. 
(*Insert here all packages use in this repo or make an install.R file*) 


```{}
install.packages("nimble")  ## to run models
install.packages("here")	## to use relative paths within the folder

install.packages("MASS")	## simulate from a mult. normal distribution
```

#### Folders 

- `data` 
	- `simulated` : contains simulated data and code for generate the data
	- `ocse` : contains OCSE data from Uto et. al
- `models`  folder containing R scripts with nimble models definition, constants and inits. Filename structure `[modelType]_[subjectModel]_[itemModel]_[ratersModel].R`. E.g., `para_2PL_1PL` stands for `[subjectModel]` is parametric, `[itemModel]` is 2PL (two parameters per item), `[ratersModel]` is 1PL (one parameter per rater)
	- `parametric` 
	- `semiparametric` 
- `output` folder that collects posterior samples (and other info) 
- `functions` helper functions to process posterior samples
- `uto_el_al` folder containing README and codes of (Uto et al. paper)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309887#references]

### Run a model 

The script `1_runNimbleModel` is a parametrized R script that can be run from /terminal. Conceptually similar to R functions, as arguments are passed to the script. Run from the terminal is going to be convenient for simulations (we will likely run those on a cluster). 


`Rscript 1_runNimbleModel.R --model= --dirResults= --data= --niter= --nburnin=`


```{}
 --model=           # path to the model code to to run  
 --dirResults=      # directory to results  
 --data=            # directory to data   
 --niter=           # number of iterations  
 --nburnin=         # number of burnin iteration  
 --nthin=			# thinning (optional) default 1
 ```

 Some conventions:

 * The model code script uses `modelCode` as variable name for the model code, and contains definitions for `data, constants, inits`
 * Results of the MCMC are saved in the `output` folder using the convention `output/[dataName]/[modelType]/[modelName.rds]`
 * Results are saved as `list` containing `samples`, timings (compilation, sampling excluding burning, running time) and MCMC settings `MCMCcontrols`