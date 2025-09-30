#
library(here)
library(MASS)

# --- Helpers --------------------------------------------------------------

# function to draw subject abilities (eta) under a chosen latent-trait scenario
draw_eta <- function(N, latent_scenario = c("unimodal","bimodal"),
                     mu_uni = 0, sd_uni = 1,
                     mix_mu = c(-2, 2), mix_sd = c(1, 1), mix_w = c(0.5, 0.5)) {

  latent_scenario <- match.arg(latent_scenario)

  if (latent_scenario == "unimodal") {
    eta <- rnorm(N, mu_uni, sd_uni)
    list(eta = eta, mu_s = mu_uni, sigma_s = sd_uni^2, c_p_s_s = 1)
  } else {
    # bimodal (finite mixture for simulation)
    comp <- sample.int(2, N, replace = TRUE, prob = mix_w)
    eta <- rnorm(N, mix_mu[comp], mix_sd[comp])
    list(eta = eta, mu_s = mix_mu, sigma_s = mix_sd^2, c_p_s_s = mix_w)
  }
}

# build design: assign each subject to exactly 'sp' raters
# and create indices vector PPi, RRi, ARi

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

generateData <- function(eta, N, R, tot, PPi, II, RRi, ARi,
                         lambda, beta, delta, R_tau, R_phi,
                         use_AR = TRUE, K) {

  y    <- numeric(tot)
  # pic: stores adjacent–category logits for each observation and category
  # pi:  stores cumulative logits (sums of adjacent logits) for each obs/category
  # ppi: stores normalized probabilities (softmax over cumulative logits)
  
  pic  <- matrix(0, nrow = tot, ncol = K)
  pi   <- matrix(0, nrow = tot, ncol = K)
  ppi  <- matrix(0, nrow = tot, ncol = K)

  # MAIN LOOP
  for (n in 1:tot) {
    # Category 1 (score 0 in the model) is the *baseline*.
    # By definition, its adjacent logit and cumulative logit are zero.
    pic[n, 1] <- 0
    pi[n, 1]  <- 0

    ar_prev_idx <- ARi[n]

    for (k in 2:K) {
      # delta now enters with a minus sign (treated as thresholds)
      base <- exp(R_phi[RRi[n]] + lambda[II[n]]) *
              (eta[PPi[n]] - beta[II[n]] - delta[k - 1] - R_tau[RRi[n], 1])

      if (use_AR && ar_prev_idx != 0) {
        # AR contribution: subtract rho_r * eta_prev (matches paper’s definition)
        base <- base - exp(R_phi[RRi[n]] + lambda[II[n]]) * (R_tau[RRi[n], 2] * eta[ar_prev_idx])
      }

      pic[n, k] <- base
      pi[n, k]  <- sum(pic[n, 1:k])
    }

    ppi[n, 1:K] <- exp(pi[n, 1:K]) / sum(exp(pi[n, 1:K]))
    y[n] <- sample(1:K, 1, prob = ppi[n, 1:K])
  }

  list(
    y = y, eta = eta, ppi = ppi,
    PPi = PPi, II = II, RRi = RRi, ARi = ARi,
    lambda = lambda, beta = beta, delta = delta,
    tau = R_tau, phi = R_phi
  )
}

# ---------------- CONFIG ----------------

N               = 100  # Subjects sample size
I               = 5    # Number of items
R               = 5    # Number of raters

K               = 4    # Number of categories
sp              = 3    # number of subjects per rater
tot             = sp * N * I

# ---- RATER FEATURES  ----
# raters features are correlated

Omega  <- diag(1, 3, 3)              #
sigma3 <- c(1, 0.3, 0.2)
Sigma  <- diag(sigma3) %*% Omega %*% diag(sigma3)
mu3    <- c(0, 0, 0)

R_features <- MASS::mvrnorm(R, mu3, Sigma)

# 1. severity 
# 2. anchoring
# 3. consistency/discrimination

R_tau <- rbind(
  R_features[1:(R-1), 1:2],
  c(-sum(R_features[, 1]), R_features[R, 2])
)
l_phi <- c(R_features[1:(R-1), 3], -sum(R_features[1:(R-1), 3]))
R_phi <- l_phi

if (any(R_tau[,2] > 1) || any(R_tau[,2] < -1)) {
  message("WARNING: divergent AR(1) path (|rho|>1)")
}

# ---- ITEM PARAMETERS (same constraints as yours) ----
l_lambda <- rnorm(I - 1, 0, 0.3); l_lambda[I] <- -sum(l_lambda[1:(I - 1)])
beta     <- rnorm(I - 1, 0, 0.5); beta[I]     <- -sum(beta[1:(I - 1)])
lambda   <- l_lambda
delta    <- rnorm(K - 1, 1, 1)

# ---- DESIGN (subjects × raters, then expand by items) ----
design <- build_design(N, R, I, sp)
PPi <- design$PPi; RRi <- design$RRi; ARi <- design$ARi; II <- design$II

# ---- A small convenience for filenames ----
write_scenario <- function(latent, ar_on) paste0(latent, if (ar_on) "_AR" else "_noAR")

# ---- Scenario A: UNIMODAL + AR ----
eta_uni <- draw_eta(N, "unimodal", mu_uni = 0, sd_uni = 1)$eta
Data_uni_AR <- generateData(
  eta = eta_uni, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda, beta = beta, delta = delta,
  R_tau = R_tau, R_phi = R_phi,
  use_AR = TRUE, K = K
)

# ---- Scenario B: UNIMODAL + NO AR ----
Data_uni_noAR <- generateData(
  eta = eta_uni, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda, beta = beta, delta = delta,
  R_tau = R_tau, R_phi = R_phi,
  use_AR = FALSE, K = K
)

# ---- Scenario C: BIMODAL + AR ----
eta_bi <- draw_eta(
  N, "bimodal",
  mix_mu = c(-2, 2), mix_sd = c(1, 1), mix_w = c(0.5, 0.5)
)$eta
Data_bi_AR <- generateData(
  eta = eta_bi, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda, beta = beta, delta = delta,
  R_tau = R_tau, R_phi = R_phi,
  use_AR = TRUE, K = K
)

# ---- Scenario D: BIMODAL + NO AR ----
Data_bi_noAR <- generateData(
  eta = eta_bi, N = N, R = R, tot = length(PPi),
  PPi = PPi, II = II, RRi = RRi, ARi = ARi,
  lambda = lambda, beta = beta, delta = delta,
  R_tau = R_tau, R_phi = R_phi,
  use_AR = FALSE, K = K
)

# ---- Save (clear, scenario-specific filenames) ----
dir.create(here::here("data", "simulated"), showWarnings = FALSE, recursive = TRUE)
saveRDS(Data_uni_AR,   here("data","simulated","data_unimodal_AR.rds"))
saveRDS(Data_uni_noAR, here("data","simulated","data_unimodal_noAR.rds"))
saveRDS(Data_bi_AR,    here("data","simulated","data_bimodal_AR.rds"))
saveRDS(Data_bi_noAR,  here("data","simulated","data_bimodal_noAR.rds"))

