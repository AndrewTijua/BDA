data {
  int<lower=1> K;
  int<lower=0> N;
  matrix[N, K] X;
  vector[N] y;
  int<lower=0> N_new;
  matrix[N_new, K] x_new;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  beta ~ normal(0, 10000);
  sigma ~ inv_gamma(0.1, 0.1);
  y ~ normal(X * beta, sigma);
}
generated quantities {
  vector[N_new] y_new;
  for (n in 1:N_new)
    y_new[n] = normal_rng(x_new[n] * beta, sigma);
}
