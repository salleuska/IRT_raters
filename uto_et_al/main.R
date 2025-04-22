library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dat_org <- read.table("dat.csv", header=TRUE, sep=",")
  
setting <- list(K = 4, n_person = 30, n_item = 30, n_rater = 5)
dat <- c()
for(i in 1:length(dat_org[,1])){
  for(j in 3:length(dat_org[1,])){
    if(!is.na(dat_org[i, j])){
      dat <- rbind(dat, cbind(dat_org[i,1:2], j-2, dat_org[i, j]+1))
    }
  }
}
colnames(dat) <- c("Items","Raters","Examinees","Score")

data_irt=list(K=setting$K, J=setting$n_person, R=setting$n_rater, I=setting$n_item, N=nrow(dat), 
              Examinees=dat$Examinees, Raters=dat$Raters, Items=dat$Items, X=dat$Score)

stan <- stan_model(file="model.stan")
fit <- sampling(stan, data=data_irt, 
  iter=5000, warmup=2000, chains=3)

get_prod_restricted_prm <- function(d){
  return(append(1.0/prod(d), d))
}

get_beta_ir <- function(prm_org, I, R){
  for(i in 1:I){
    prm = prm_org[((i-1)*R+1):((i-1)*R+R)]
    if(i == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}

get_category_prm <- function(category_prm, N, K){
  for(n in 1:N){
    prm = category_prm[((n-1)*(K-2)+1):((n-1)*(K-2)+(K-2))]
    prm = append(prm, -1*sum(prm))
    if(n == 1){
      mat = t(data.frame(prm))
    } else {
      mat = rbind(mat, t(data.frame(prm)))
    }
  }  
  return(mat)
}


saveRDS("fit", file="fit.rds")
alpha_r <- get_prod_restricted_prm(summary(fit, par="alpha_r")$summary[,"mean"])
d_rk <- get_category_prm(summary(fit, par="beta_rk")$summary[,"mean"], setting$n_rater, setting$K)
alpha_i <- summary(fit, par="alpha_i")$summary[,"mean"]
d_ik <- get_category_prm(summary(fit, par="beta_ik")$summary[,"mean"], setting$n_item, setting$K)
beta_ir <- get_beta_ir(summary(fit, par="beta_ir")$summary[,"mean"], setting$n_item, setting$n_rater)
theta <- summary(fit, par="theta")$summary[,"mean"]

fit_ss <- extract(fit, permuted = TRUE) # fit_ss is a list 
theta <- fit_ss$theta

hist(theta, breaks = 100)