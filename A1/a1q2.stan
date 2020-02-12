data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  simplex[K] theta;          // mixing proportions
  positive_ordered[K] alpha; // locations of mixture components
  vector<lower=0>[K] beta;  // scales of mixture components
  
}
parameters {
  real<lower=0> lambda;
}
model {
  vector[K] log_theta = log(theta);  // cache log calculation
  y ~ exponential(lambda);
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K)
      lps[k] += gamma_lpdf(lambda | alpha[k], beta[k]);
    target += log_sum_exp(lps);
  }
}
generated quantities {
  real<lower=0> mwt = 1 / lambda;
  real postdraw = exponential_rng(lambda);
}
