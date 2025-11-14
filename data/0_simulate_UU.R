#####################################################
library(here)
#####################################################################
# function to draw subject abilities (eta) under a chosen latent-trait scenario
#####################################################################
draw_eta <- function(N, latent_scenario = c("unimodal","bimodal"),
                     mu_uni = 0, sd_uni = 1,
                     mix_mu = c(-2, 2), mix_sd = c(1, 1), mix_w = c(0.5, 0.5)) {

  latent_scenario <- match.arg(latent_scenario)

  if (latent_scenario == "unimodal") {
    eta <- rnorm(N, mu_uni, sd_uni)
    list(eta = eta, mu_s = mu_uni, sigma_s = sd_uni^2)
  } else {
    # bimodal (finite mixture for simulation)
    comp <- sample.int(2, N, replace = TRUE, prob = mix_w)
    eta <- rnorm(N, mix_mu[comp], mix_sd[comp])
    list(eta = eta, mu_s = mix_mu, sigma_s = mix_sd^2)
  }
}

#####################################################################
# build design: assign each subject to exactly 'sp' raters
# and create indices vector PPi, RRi, ARi

## @N  : number of subjects
## @R  : number of raetehrs
## @I  : number of items
## @sp : subjects per raters
#####################################################################

build_design <- function(N, R, I, sp) {
  Q_mat <- matrix(0, nrow = N, ncol = R)  # subject-by-rater assignment
  # (balanced-ish) each subject gets 'sp' raters
  for (j in 1:N) {
    ind <- if (sp >= R) 1:R else sample(1:R, sp, replace = FALSE)
    Q_mat[j, ind] <- 1
  }

  # For each rater, get the ordered list of subjects they rate, then set previous (AR)
  PP <- integer(0)  # subject id per rater "block"
  RR <- integer(0)  # rater id per rater "block"
  AR <- integer(0)  # previous subject (same rater), or 0 if none

  for (r in 1:R) {
    subj_r <- which(Q_mat[, r] == 1)
    # if you have a known scoring order per rater, sort/permute here
    if (length(subj_r) > 0) {
      PP <- c(PP, subj_r)
      RR <- c(RR, rep(r, length(subj_r)))
      AR <- c(AR, c(0, subj_r[-length(subj_r)]))
    }
  }

  # expand by items
  PPi <- rep(PP, each = I)
  RRi <- rep(RR, each = I)
  ARi <- rep(AR, each = I)
  II  <- rep(rep(1:I, times = 1), times = length(PP))

  stopifnot(length(PPi) == length(RRi), length(PPi) == length(ARi), length(PPi) == length(II))
  list(Q_mat = Q_mat, PPi = PPi, RRi = RRi, ARi = ARi, II = II)
}

## -----------------------------------------------------------
## 3. True parameters under the SAME constraints as NIMBLE
## -----------------------------------------------------------
gen_true_params_uto <- function(I, R, K,
                                meanlog_alpha_i = 0, sdlog_alpha_i = 0.5,
                                meanlog_alpha_r = 0, sdlog_alpha_r = 0.5,
                                sd_beta_ir = 1, sd_steps = 1) {
  
  alpha_i <- rlnorm(I, meanlog_alpha_i, sdlog_alpha_i)
  
  alpha_r_free <- rlnorm(R - 1, meanlog_alpha_r, sdlog_alpha_r)
  trans_alpha_r <- numeric(R)
  trans_alpha_r[2:R] <- alpha_r_free
  trans_alpha_r[1]   <- 1 / prod(alpha_r_free)
  
  beta_ir <- matrix(rnorm(I * R, 0, sd_beta_ir), nrow = I, ncol = R)
  
  d_rk <- matrix(0, nrow = R, ncol = K)
  category_est_r <- matrix(0, nrow = R, ncol = K)
  for (r in 1:R) {
    d_rk[r, 1] <- 0
    if (K > 2) d_rk[r, 2:(K - 1)] <- rnorm(K - 2, 0, sd_steps)
    category_est_r[r, 1:(K - 1)] <- d_rk[r, 1:(K - 1)]
    category_est_r[r, K]         <- -sum(d_rk[r, 2:(K - 1)])
  }
  
  d_ik <- matrix(0, nrow = I, ncol = K)
  category_est_i <- matrix(0, nrow = I, ncol = K)
  for (i in 1:I) {
    d_ik[i, 1] <- 0
    if (K > 2) d_ik[i, 2:(K - 1)] <- rnorm(K - 2, 0, sd_steps)
    category_est_i[i, 1:(K - 1)] <- d_ik[i, 1:(K - 1)]
    category_est_i[i, K]         <- -sum(d_ik[i, 2:(K - 1)])
  }
  
  list(
    alpha_i        = alpha_i,
    trans_alpha_r  = trans_alpha_r,
    beta_ir        = beta_ir,
    d_rk           = d_rk,
    d_ik           = d_ik,
    category_est_r = category_est_r,
    category_est_i = category_est_i
  )
}


#####################################################################
## Function that generates data (raters' scores for each subject/item)
## given latent abilities (unimodal, the desing, and item parameters
## Can include the autoregressive scoring or not

## The output is a dataset in "long" format, meaning that for each row we have
## score, rater id, subject id, item id
#####################################################################

generateData_uto <- function(eta, N, I, R, K,
                             PPi, II, RRi,
                             alpha_i, trans_alpha_r,
                             beta_ir, category_est_r, category_est_i) {

  tot <- length(PPi)
  stopifnot(tot == length(II), tot == length(RRi))

  y    <- integer(tot)
  pic  <- matrix(0, nrow = tot, ncol = K)  # adjacent–category logits
  pi   <- matrix(0, nrow = tot, ncol = K)  # cumulative logits
  ppi  <- matrix(0, nrow = tot, ncol = K)  # probabilities

  for (n in 1:tot) {
    j <- PPi[n]
    i <- II[n]
    r <- RRi[n]

    # category 1 reference
    pic[n, 1] <- 0
    pi[n, 1]  <- 0

    # φ_ir = α_i * α_r  (no D)
    phi_ir <- alpha_i[i] * trans_alpha_r[r]

    for (k in 2:K) {
      # same linear predictor as in your NIMBLE code:
      # pic[n,k] <- φ_ir * (η_j - β_ir - cat_r - cat_i)
      lin <- phi_ir * (eta[j] - beta_ir[i, r] -
                         category_est_r[r, k] -
                         category_est_i[i, k])

      pic[n, k] <- lin
      pi[n, k]  <- pi[n, k - 1] + lin
    }

    # stable softmax over cumulative logits
    z  <- pi[n, 1:K] - max(pi[n, 1:K])
    ez <- exp(z)
    ppi[n, 1:K] <- ez / sum(ez)

    y[n] <- sample.int(K, size = 1, prob = ppi[n, 1:K])
  }

  list(
    y = y,
    eta = eta,
    pic = pic,
    pi = pi,
    ppi = ppi
  )
}

## -----------------------------------------------------------
## 5. Unimodal ability scenario
## -----------------------------------------------------------
set.seed(99)

N <- 100
I <- 5
R <- 5
K <- 4
sp <- 3

design <- build_design(N, R, I, sp)
PPi <- design$PPi
RRi <- design$RRi
II  <- design$II

eta_uni <- draw_eta(N, "unimodal", 0, 1)$eta

true_par <- gen_true_params_uto(
  I = I, R = R, K = K,
  sd_beta_ir = 0.5,
  sd_steps   = 0.5
)

Data_uni <- generateData_uto(
  eta = eta_uni,
  N = N, I = I, R = R, K = K,
  PPi = PPi, RRi = RRi, II = II,
  alpha_i        = true_par$alpha_i,
  trans_alpha_r  = true_par$trans_alpha_r,
  beta_ir        = true_par$beta_ir,
  category_est_r = true_par$category_est_r,
  category_est_i = true_par$category_est_i
)

## ---- Save everything in one object ----
out <- list(
  y              = Data_uni$y,
  eta_true       = eta_uni,
  PPi = PPi, RRi = RRi, II = II,
  alpha_i_true        = true_par$alpha_i,
  trans_alpha_r_true  = true_par$trans_alpha_r,
  beta_ir_true        = true_par$beta_ir,
  category_est_r_true = true_par$category_est_r,
  category_est_i_true = true_par$category_est_i
)

dir.create(here("data","simulated"), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, here("data","simulated","uto_sim_unimodal.rds"))

## -----------------------------------------------------------
## 6. Bimodal ability scenario 
## -----------------------------------------------------------

## bimodal abilities (finite Normal mixture for simulation)
eta_bi <- draw_eta(
  N,
  latent_scenario = "bimodal",
  mix_mu = c(-2, 2),
  mix_sd = c(1, 1),
  mix_w  = c(0.5, 0.5)
)$eta

## Center the generated abiltities
eta_bi <- (eta_bi - mean(eta_bi))/sd(eta_bi)

## generate data with the SAME true item/rater parameters
Data_bi <- generateData_uto(
  eta = eta_bi,
  N = N, I = I, R = R, K = K,
  PPi = PPi, RRi = RRi, II = II,
  alpha_i        = true_par$alpha_i,
  trans_alpha_r  = true_par$trans_alpha_r,
  beta_ir        = true_par$beta_ir,
  category_est_r = true_par$category_est_r,
  category_est_i = true_par$category_est_i
)

## pack and save
out_bi <- list(
  y              = Data_bi$y,
  eta_true       = eta_bi,
  PPi = PPi, RRi = RRi, II = II,
  alpha_i_true        = true_par$alpha_i,
  trans_alpha_r_true  = true_par$trans_alpha_r,
  beta_ir_true        = true_par$beta_ir,
  category_est_r_true = true_par$category_est_r,
  category_est_i_true = true_par$category_est_i
)

saveRDS(out_bi, here("data","simulated","uto_sim_bimodal.rds"))