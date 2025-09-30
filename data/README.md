- `data` 
	- `simulated` : contains simulated data and code for generate the data
	- `ocse` : contains OCSE data from Uto et. al



`0_generate_ARdata`

Simulates raters scores according to our model. 
We assume an autoregressive process where the rater's scoring is influences by the ability of the previous scored subject


### Simulation dimensions

N               = 100   # Subjects sample size
I               = 5     # Number of items
R               = 5     # Number of raters

K               = 4     # Number of categories
C               = 1     # Number of students' clusters

sp              = 3     # number of subjects per rater

tot             = sp * N * I


### Item parameters description


### Subject latent ability distribution


### Raters' Features Distribution

Raters' features: systematic bias, autoregressive coefficient, consistency. 

## `data`

- **`simulated/`**  
  Contains simulated datasets and R scripts to generate them.  
  - `0_generate_ARdata.R`: simulates rater–subject–item responses under the Extended Generalized Many-Facet Rasch Model (E-GMFRM) and its dynamic extension with autoregressive (anchoring) effects.  
  - Output `.rds` files include multiple scenarios (see below).  

- **`ocse/`**  
  Contains the OSCE data from Uto et al. (2024), used as the real-data application.

---

## `0_generate_ARdata.R`

This script simulates rater scores according to the **E-GMFRM** and **Dynamic E-GMFRM**.  
We allow for both **autoregressive effects** — where a rater’s severity for the current subject depends on the latent trait of the previously scored subject — and scenarios without autoregression.

### Simulation Scenarios

The script generates four datasets by default:

1. **Unimodal latent trait, with AR effect** → `data_unimodal_AR.rds`  
2. **Unimodal latent trait, without AR effect** → `data_unimodal_noAR.rds`  
3. **Bimodal latent trait, with AR effect** → `data_bimodal_AR.rds`  
4. **Bimodal latent trait, without AR effect** → `data_bimodal_noAR.rds`

### Simulation Dimensions (defaults)

- **N** = 100 subjects  
- **I** = 5 items  
- **R** = 5 raters  
- **K** = 4 score categories (0,…,K)  
- **sp** = 3 raters per subject  
- **tot** = `sp * N * I` total observations  

### Parameters

#### Item Parameters
- **lambda[i]**: item discrimination parameters (constrained to sum to zero).  
- **beta[i]**: item location/easiness parameters (constrained to sum to zero).  
- **delta[k]**: category step thresholds shared across raters.  

#### Subject Latent Traits
- **eta[j]**: latent ability of subject *j*.  
  - *Unimodal*: eta[j] ~ Normal(0,1)  
  - *Bimodal*: eta[j] ~ 0.5·Normal(-2,1) + 0.5·Normal(2,1)  

#### Rater Features
For each rater *r*:
- **tau[r,1]**: global bias (severity).  
- **tau[r,2]**: autoregressive coefficient (anchoring effect of previous subject’s ability).  
- **phi[r]**: discrimination/consistency parameter.  

Rater features are sampled jointly from a multivariate normal distribution with identifiability constraints (sum-to-zero).

---

### Output

Each `.rds` file contains a list with:
- `y`: simulated categorical responses  
- `eta`: subject latent traits  
- `PPi`, `II`, `RRi`, `ARi`: indices linking responses to subjects, items, raters, and previous subjects  
- item parameters (`lambda`, `beta`, `delta`)  
- rater parameters (`tau`, `phi`)  


