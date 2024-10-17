# Medical Interview Score Data from PostCC-OSCE and Programs for an Extended Many-Facet IRT Model

[https://doi.org/10.5061/dryad.tmpg4f56q](https://doi.org/10.5061/dryad.tmpg4f56q)

This dataset includes the actual score data collected from a medical interview test conducted at Tokyo Medical and Dental University as part of a PostCC OSCE, as well as the program for estimating parameters of the new item response theory model.

## **Structure of data and related programs:**

This dataset includes the following data:

* **dat.csv:** The actual OSCE score data.

It also includes the following related programs:

* **main.R:** The main program to estimate the proposed IRT model parameters.
* **model.stan:** The Stan code for the parameter estimation of the proposed model.

## Description of the data

This is the actual score data collected from a medical interview test. The data consist of rating scores assigned by five raters to 30 videos recorded as material for reconfirmation of medical interview tests for sixth-year dental students. Of the five raters, two each evaluated all 30 examinees, while the remaining three each evaluated 10 different examinees. The raters performed their evaluations independently, without communicating or consulting with one another. The 30 evaluation items on a scoring rubric were scored on a 4-point scale (4: excellent; 3: good; 2: acceptable; 1: unsatisfactory).

* The first column of the **dat.csv** file represents the index of the rubric's evaluation items, starting from 1.
* The second column represents the index of raters, starting from 1.
* The third column to the last column represents the examinees' responses. Each response must be a continuous integer value starting from 0. **The empty cells represent missing data**.

## **Related programs and u**sage notes

The following two main codes are included. Note that to run these codes, R and RStan environments are required. Furthermore, they must be located at the same folder level as the **dat.csv** file.

* **model.stan:** The Stan code for the parameter estimation of the proposed model.
* **main.R:** The main program to estimate the proposed IRT model parameters through the Markov chain Monte Carlo algorithm.

To estimate the parameters from **dat.csv**, simply run the entire code in the **main.R** file. The parameter estimates will be stored in the variables "alpha_r", "alpha_i", "beta_ir", "beta_ik", "d_rk", and "theta". The variable "alpha_r" represents the consistency of the r-th rater, "alpha_i" represents the discrimination power of the i-th item, "beta_ir" represents the rater-item interaction parameter that indicates the severity of rater r on item i, "beta_ik" represents the item-specific step-difficulty parameter that represents the difficulty of transition between scores k-1 and k for item i, "d_rk" represents the rater-specific step-difficulty parameter denoting the severity of rater r of transition from score k-1 to k, and "theta" represents the latent ability of examinee j.

To apply our model to another dataset, please create a dataset with the same structure as **dat.csv** and update the file name specified in line 6 of **main.R**, as well as the number of examinees, raters, and rubric items specified in line 8. "K" in line 8 refers to the number of score categories.
