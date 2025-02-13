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
- `models`  folder containing R scripts with nimble models definition, constants and inits. Filename structure `[subjectModel]_[itemModel]_[ratersModel].R`. E.g., `para_2PL_1PL` stands for `[subjectModel]` is parametric, `[itemModel]` is 2PL (two parameters per item), `[ratersModel]` is 1PL (one parameter per rater)
	- `parametric` 
	- `semiparametric` 
- `output` folder that collects posterior samples (and other info) 
- `functions` helper functions to process posterior samples
- `uto_el_al` folder containing README and codes of (Uto et al. paper)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309887#references]
