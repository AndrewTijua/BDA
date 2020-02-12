data {          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  simplex[2] theta;          // mixing proportions
  positive_ordered[2] alpha; // locations of mixture components
  vector<lower=0>[2] beta;  // scales of mixture components
  
}
parameters {
  real<lower=0> lambda;
}
model {
  y ~ exponential(lambda);
  target += log_mix(theta[1],
  gamma_lpdf(lambda | alpha[1], beta[1]),
  gamma_lpdf(lambda | alpha[2], beta[2]));
}
generated quantities {
  real<lower=0> mwt = 1 / lambda;
  real postdraw = exponential_rng(lambda);
}
