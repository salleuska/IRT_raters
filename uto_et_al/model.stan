data{
  int <lower=0> J;//n_examinee
  int <lower=0> I;//n_rubric_item
  int <lower=0> R;//n_rater
  int <lower=2> K;//n_score_category
  int <lower=0> N;//n_samples
  int <lower=1, upper=J> Examinees [N];
  int <lower=1, upper=I> Items [N];
  int <lower=1, upper=R> Raters [N];
  int <lower=1, upper=K> X [N];
}
transformed data{
  vector[K] c = cumulative_sum(rep_vector(1, K)) - 1;
}
parameters {
  real theta[J];
  real<lower=0> alpha_r [R-1];
  real<lower=0> alpha_i [I];
  vector[R] beta_ir[I];
  vector[K-2] beta_rk [R];
  vector[K-2] beta_ik [I];
}
transformed parameters{
  real<lower=0> trans_alpha_r[R];
  vector[K-1] category_est_r [R];
  vector[K-1] category_est_i [I];
  vector[K] category_prm[I, R];
  trans_alpha_r[1] = 1.0 / prod(alpha_r);
  trans_alpha_r[2:R] = alpha_r;  
  for(r in 1:R){
    category_est_r[r, 1:(K-2)] = beta_rk [r];
    category_est_r[r, K-1] = -1*sum(beta_rk [r]);  
  }
  for(z in 1:I){
    category_est_i[z, 1:(K-2)] = beta_ik [z];
    category_est_i[z, K-1] = -1*sum(beta_ik [z]);  
    for(r in 1:R){
      category_prm[z, r] = cumulative_sum(append_row(0, (category_est_r[r] + category_est_i[z])));    
    }
  }
}
model{
  theta ~ normal(0, 1);
  alpha_i ~ lognormal(0, 0.5);
  trans_alpha_r ~ lognormal(0, 0.5);
  for (z in 1:I) beta_ir[z,] ~ normal(0, 1);
  for (r in 1:R) category_est_r [r,] ~ normal(0, 1);
  for (z in 1:I) category_est_i [z,] ~ normal(0, 1);
  for (n in 1:N){
    X[n] ~ categorical_logit(1.7 * trans_alpha_r[Raters[n]] * alpha_i[Items[n]] 
    * ( c * (theta[Examinees[n]] - beta_ir[Items[n],Raters[n]]) -category_prm[Items[n], Raters[n]]));
  }
}
generated quantities {
  vector[N] log_lik;
  int<lower=1, upper=K> X_tilde[N];
  for (n in 1:N){
    log_lik[n] = categorical_logit_log(X[n], 1.7 * trans_alpha_r[Raters[n]] * alpha_i[Items[n]]
    * ( c * (theta[Examinees[n]] - beta_ir[Items[n],Raters[n]]) -category_prm[Items[n], Raters[n]]));
    X_tilde[n] = categorical_logit_rng(1.7 * trans_alpha_r[Raters[n]] * alpha_i[Items[n]] 
    * ( c * (theta[Examinees[n]] - beta_ir[Items[n],Raters[n]]) -category_prm[Items[n], Raters[n]]));
  }
}
