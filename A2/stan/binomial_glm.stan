data {
  int<lower=0> N;
  int<lower=0> P;
  
  int<lower=0> y[N];
  
  matrix[N, P] X;
  
  vector[2] n_params;
}

parameters {
  vector[P] beta;
}

transformed parameters{
  vector[N] lg_p = X * beta;
}

model {
  beta ~ normal(n_params[1], n_params[2]);
  y ~ binomial(1, inv_logit(lg_p));
}
generated quantities{
  int data_ppred[N] = binomial_rng(1, inv_logit(lg_p));
}

