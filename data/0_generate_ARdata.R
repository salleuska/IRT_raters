#####################################################
library(here)
library(MASS)

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


#####################################################################
## Function that generates data (raters' scores for each subject/item)
## given latent abilities (unimodal, the desing, and item parameters
## Can include the autoregressive scoring or not

## The output is a dataset in "long" format, meaning that for each row we have
## score, rater id, subject id, item id
#####################################################################

generateData <- function(eta, N, R, tot, PPi, II, RRi, ARi,
                         lambda, beta_ir, rho_rk, delta_ik, R_tau, R_phi,
                         use_AR = FALSE, K) {

  y    <- numeric(tot)
  pic  <- matrix(0, nrow = tot, ncol = K)  # adjacent–category logits
  pi   <- matrix(0, nrow = tot, ncol = K)  # cumulative logits
  ppi  <- matrix(0, nrow = tot, ncol = K)  # probabilities


  # MAIN LOOP
  for (n in 1:tot) {
    pic[n, 1] <- 0
    pi[n, 1]  <- 0

    ar_prev_idx <- ARi[n]

    for (k in 2:K) {
      step <- k -1

      # δ_{jirk} = β_{ir} + b_{ik} + ρ_{rk}
      delta_step <- beta_ir[II[n], RRi[n]] +  delta_ik[II[n], step] + rho_rk[RRi[n], step]
      base <- exp(R_phi[RRi[n]] + lambda[II[n]]) * (eta[PPi[n]] - delta_step)

      if (use_AR) {
        stop("AR contribution not implemented in this generator. Set use_AR=FALSE.")
        # If re-enabling later:
        # if (ar_prev_idx != 0) {
        #   base <- base - exp(R_phi[RRi[n]] + lambda[II[n]]) *
        #                   (R_tau[RRi[n], 2] * eta[ar_prev_idx])
        # }
      }

      pic[n, k] <- base
      pi[n, k]  <- pi[n, k - 1] + base
    }

    # stable softmax over cumulative logits
    z <- pi[n, 1:K] - max(pi[n, 1:K])
    ez <- exp(z)
    ppi[n, 1:K] <- ez / sum(ez)

    y[n] <- sample(1:K, 1, prob = ppi[n, 1:K])
  }

  list(
    y = y, eta = eta, ppi = ppi,
    PPi = PPi, II = II, RRi = RRi, ARi = ARi,
    lambda = lambda, beta_ir = beta_ir, delta_ik = delta_ik,
    tau = R_tau, phi = R_phi, rho_rk = rho_rk, K = K
  )
}

## ---------------- CONFIG (all-in-one) ----------------

set.seed(1)

N <- 100     # subjects
I <- 5       # items
R <- 5       # raters
K <- 4       # categories (scored 0..K-1; generator uses 1..K)
sp <- 3      # raters per subject
tot <- sp * N * I

## ---- RATER FEATURES (severity, AR anchor, consistency) ----
Omega  <- diag(1, 3, 3)
sigma3 <- c(1, 0.3, 0.2)
Sigma  <- diag(sigma3) %*% Omega %*% diag(sigma3)
mu3    <- c(0, 0, 0)

R_features <- MASS::mvrnorm(R, mu3, Sigma)
# 1: severity, 2: anchoring, 3: consistency (log-ξ)
R_tau <- rbind(
  R_features[1:(R-1), 1:2],
  c(-sum(R_features[, 1]), R_features[R, 2])
)
l_phi <- c(R_features[1:(R-1), 3], -sum(R_features[1:(R-1), 3]))
R_phi <- l_phi
if (any(R_tau[,2] > 1) || any(R_tau[,2] < -1)) message("WARNING: divergent AR(1) path (|rho|>1)")

## ---- ITEM PARAMETERS ----
l_lambda <- rnorm(I - 1, 0, 0.3); l_lambda[I] <- -sum(l_lambda[1:(I - 1)])  # sum-to-zero on log-λ
lambda   <- l_lambda
beta_i   <- rnorm(I - 1, 0, 0.5); beta_i[I] <- -sum(beta_i[1:(I - 1)])      # item main effect
delta    <- rnorm(K - 1, 1, 1)                                              # base step vector

## ---- BUILD MODEL MATRICES EXPECTED BY generateData ----
# β_ir (I x R): item–rater base difficulty (start from item + rater severity)
beta_ir <- outer(beta_i, rep(1, R)) + matrix(R_tau[, 1], nrow = I, ncol = R, byrow = TRUE)

# ρ_rk (R x (K-1)): rater step adjustments (start at 0)
rho_rk <- matrix(0, nrow = R, ncol = K - 1)

# Δ_ik (I x (K-1)): item-specific steps (replicate delta across items)
delta_ik <- matrix(delta, nrow = I, ncol = K - 1, byrow = TRUE)

## ---- LIGHT IDENTIFIABILITY (center steps and β_ir interaction) ----
# center steps within item
delta_ik <- sweep(delta_ik, 1, rowMeans(delta_ik), FUN = "-")
# center rater step effects within rater
rho_rk   <- sweep(rho_rk, 1, rowMeans(rho_rk), FUN = "-")
# make β_ir a pure interaction (remove row/col means)
beta_ir  <- beta_ir -
            matrix(rowMeans(beta_ir), nrow = I, ncol = R, byrow = FALSE) -
            matrix(colMeans(beta_ir), nrow = I, ncol = R, byrow = TRUE) +
            mean(beta_ir)

## ---- DESIGN (subjects × raters, then expand by items) ----
design <- build_design(N, R, I, sp)   # uses your existing helper
PPi <- design$PPi; RRi <- design$RRi; ARi <- design$ARi; II <- design$II

## ---- ABILITIES ----
eta_uni <- rnorm(N, 0, 1)

## ---- SIMULATE (NO AR) ----
Data_uni_noAR <- generateData(
  eta = eta_uni, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda,                # log-scale λ_i (φ_ir = exp(λ_i + R_φ_r))
  beta_ir = beta_ir,              # I x R
  rho_rk  = rho_rk,               # R x (K-1)
  delta_ik = delta_ik,            # I x (K-1)
  R_tau = R_tau, R_phi = R_phi,
  use_AR = FALSE, K = K
)

eta_bi <- draw_eta(
  N, "bimodal",
  mix_mu = c(-2, 2), mix_sd = c(1, 1), mix_w = c(0.5, 0.5)
)$eta
# ---- Scenario: BIMODAL + NO AR ----
Data_bi_noAR <- generateData(
  eta = eta_uni, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda,                # log-scale λ_i (φ_ir = exp(λ_i + R_φ_r))
  beta_ir = beta_ir,              # I x R
  rho_rk  = rho_rk,               # R x (K-1)
  delta_ik = delta_ik,            # I x (K-1)
  R_tau = R_tau, R_phi = R_phi,
  use_AR = FALSE, K = K)

# ---- Save (clear, scenario-specific filenames) ----
dir.create(here::here("data", "simulated"), showWarnings = FALSE, recursive = TRUE)
saveRDS(Data_uni_noAR, here("data","simulated","data_unimodal_noAR.rds"))
saveRDS(Data_bi_noAR,  here("data","simulated","data_bimodal_noAR.rds"))

