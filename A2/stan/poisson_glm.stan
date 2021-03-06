data {
  int<lower=0> N;
  int<lower=0> P;
  
  int<lower=0> y[N];
  
  matrix[N, P] X;
  
  int<lower=0> N_new;
  matrix[N_new, P] X_new;
  
  vector[2] n_params;
  
  vector[N] offset;
  vector[N_new] offset_new;
}
transformed data{
}

parameters {
  vector[P] lambda;
}

transformed parameters{
  vector[N] log_rate = X * lambda + offset;
  vector[N_new] log_rate_new = X_new * lambda + offset_new;
  vector<lower=0>[N] rate = exp(log_rate);
}

model {
  lambda ~ normal(n_params[1], n_params[2]);
  y ~ poisson_log(log_rate);
}

generated quantities{
  int<lower=0> y_new[N_new] = poisson_log_rng(log_rate_new);
  int<lower=0> data_ppred[N] = poisson_log_rng(log_rate);
}
