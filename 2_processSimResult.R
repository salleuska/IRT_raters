##-----------------------------------------#
## This script process posterior samples
## and estimate quantites for inference
##-----------------------------------------#
library(here)
source(here("functions/estimateDPdensity.R"))
source(here("functions/utils.R"))
##-----------------------------------------#
args <- R.utils::commandArgs(asValue=TRUE)
## --resFileName
## --simFileName
##-----------------------------------------#
## TMP for testing
args <- list(resFileName = "output/uto_sim_bimodal/semi/semi_uto.rds", 
			 simFileName = "data/simulated/uto_sim_bimodal.rds")
##-----------------------------------------#
# if(is.null(args$outDir)) outDir <- "output/posterior_samples_elaborated/" else dir <- args$outDir

listLength <- length(strsplit(args$resFileName, "\\/|.rds")[[1]])
data <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength -2]
fileName <- strsplit(args$resFileName, "\\/|.rds")[[1]][listLength]

## modelType - parametric or semi
modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]
# ## 
# modelType <- strsplit(basename(fileName), "\\_|.rds")[[1]][1]

## read objects
resObj <- readRDS(args$resFileName)

simData <- readRDS(args$simFileName)
##-------------------------------------------------------##
## TO DO: likely some rescaling for comparison - parametric and semiparametric
## SP: Need to be sure on the constraints?
##-------------------------------------------------------##
samples <- resObj$samples
cols <- extract_param_cols(colnames(samples),
                           param_names = c("eta", "beta_ir", "alpha_i",
                                           "trans_alpha_r", "category_est_r",
                                           "category_est_i"))

cols
##-------------------------------------------------------##
## TMP - Chatgpt stuff
##-------------------------------------------------------##
## ------------------------------------------------------------ ##
## 2. Extract columns for your parameters
## ------------------------------------------------------------ ##
samples <- resObj$samples
cn <- colnames(samples)

param_cols <- extract_param_cols(
    cn,
    c("eta", "beta_ir", "category_est_i", "category_est_r",
      "alpha_i", "trans_alpha_r")
)

## ------------------------------------------------------------ ##
## 3. Load simulated true parameters and infer dims
## ------------------------------------------------------------ ##

I <- nrow(simData$beta_ir_true)
R <- ncol(simData$beta_ir_true)
K <- ncol(simData$category_est_i_true)
P <- length(simData$eta_true)

## ------------------------------------------------------------ ##
## 4. Posterior means (vectors)
## ------------------------------------------------------------ ##
eta_mean   <- colMeans(samples[, param_cols$eta, drop = FALSE])
beta_mean  <- colMeans(samples[, param_cols$beta_ir, drop = FALSE])
cat_i_mean <- colMeans(samples[, param_cols$category_est_i, drop = FALSE])
cat_r_mean <- colMeans(samples[, param_cols$category_est_r, drop = FALSE])
alpha_i_mean   <- colMeans(samples[, param_cols$alpha_i, drop = FALSE])
alpha_r_mean   <- colMeans(samples[, param_cols$trans_alpha_r, drop = FALSE])

## ------------------------------------------------------------ ##
## 5. Reshape matrices to (I,R,K)
## ------------------------------------------------------------ ##

beta_mean_mat  <- matrix(beta_mean,  nrow = I, ncol = R, byrow = FALSE)
cat_i_mean_mat <- matrix(cat_i_mean, nrow = I, ncol = K, byrow = FALSE)
cat_r_mean_mat <- matrix(cat_r_mean, nrow = R, ncol = K, byrow = FALSE)

## ------------------------------------------------------------ ##
## 6. Comparison plots  (true vs posterior mean)
## ------------------------------------------------------------ ##

compare_scatter <- function(true, est, main = "") {
    plot(true, est, pch = 19, cex = 0.6,
         xlab = "True", ylab = "Posterior mean", main = main)
    abline(0, 1, col = "red", lwd = 2)
}

# β_ir
compare_scatter(as.vector(simData$beta_ir_true),
                as.vector(beta_mean_mat),
                main = "beta_ir (item–rater)")

# Item thresholds
compare_scatter(as.vector(simData$category_est_i_true),
                as.vector(cat_i_mean_mat),
                main = "Item thresholds")

# Rater thresholds
compare_scatter(as.vector(simData$category_est_r_true),
                as.vector(cat_r_mean_mat),
                main = "Rater thresholds")

# α_i
compare_scatter(simData$alpha_i_true,
                alpha_i_mean,
                main = "Item discrimination alpha_i")

# α_r (trans_alpha_r)
compare_scatter(simData$trans_alpha_r_true,
                alpha_r_mean,
                main = "Rater discrimination trans_alpha_r")

# η (abilities)
eta_centered <- (simData$eta_true - mean(simData$eta_true))/sd(simData$eta_true)
compare_scatter(eta_centered,
                eta_mean,
                main = "Eta (abilities)")

## ------------------------------------------------------------ ##
## 7. DP density check (your original code, adapted)
## ------------------------------------------------------------ ##

out <- estimateDPdensity(
  samples      = samples,
  paramNames   = list(
    alpha   = "alpha_dp",   # <- concentration parameter (NOT eta)
    muTilde = "muTilde",     # cluster means
    s2Tilde = "s2Tilde",     # cluster variances
    zi      = "zi"           # cluster labels
  ),
  nIndividuals = P,
  grid         = seq(-8, 8, length = 200)
)


res <- data.frame(
  grid = out$grid,
  mean = apply(out$densitySamples, 2, mean)
)



library(ggplot2)

ggplot(res, aes(x = grid, y = mean)) +
  geom_line(aes(color = "Posterior (mean)"), linewidth = 1) +
  geom_density(
    data = data.frame(eta = eta_centered),
    aes(x = eta, y = after_stat(density), color = "Simulated (kernel)"),
    inherit.aes = FALSE, linewidth = 0.8, alpha = 0.7
  ) +
  scale_color_manual(values = c(
    "Posterior (mean)"   = "black",
    "Simulated (kernel)" = "steelblue"
  )) +
  labs(x = "Ability", y = "Density", color = "") +
  theme_minimal()