data {
  int<lower=1> K;
  int<lower=0> N;
  matrix[N, K] X;
  vector[N] y;
  int<lower=0> N_new;
  matrix[N_new, K] x_new;
  vector[4] p_params;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
  real<lower=0> nu;
}
model {
  nu ~ exponential(p_params[4]);
  beta ~ normal(0, p_params[1]);
  sigma ~ inv_gamma(p_params[2], p_params[3]);
  y ~ student_t(nu+1, X * beta, sigma);
}
generated quantities {
  vector[N_new] y_new;
  for (n in 1:N_new)
    y_new[n] = student_t_rng(nu+1, x_new[n] * beta, sigma);
}
