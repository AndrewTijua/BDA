data {
  int<lower=0> N;
  int<lower=0> P_f;
  int<lower=0> P_r;
  
  int<lower=0> success[N];
  int<lower=1> trials[N];
  
  matrix[N, P_f] X_f;
  matrix[N, P_r] X_r;
  
  vector[2] n_params;
}

parameters {
  vector[P_f] beta_f_raw;
  vector[P_r] sn_vec;
  real<lower=0,upper=10> reff_sdv;
}

transformed parameters{
  vector[P_f] beta_f = n_params[2] * beta_f_raw + n_params[1]; 
  vector[P_r] beta_r = reff_sdv * sn_vec;
  vector[N] lg_p = X_f * beta_f + X_r * beta_r;
}

model {
  reff_sdv ~ uniform(0, 10);
  sn_vec ~ std_normal(); //hence beta_r ~ normal(0, reff_sdv)
  beta_f_raw ~ std_normal(); //hence beta_f ~ normal(n_params[1], n_params[2])
  //beta_f ~ normal(n_params[1], n_params[2]);
  success ~ binomial(trials, inv_logit(lg_p));
}
generated quantities{
  int data_ppred[N] = binomial_rng(trials, inv_logit(lg_p));
  vector[N] data_prop = inv_logit(lg_p);
}

