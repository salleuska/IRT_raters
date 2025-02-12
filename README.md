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
- `models`  folder containing R script with nimble models definition (parametric and semiparametric)
- `uto_el_al` folder containing README and codes of (Uto et al. paper)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0309887#references]
