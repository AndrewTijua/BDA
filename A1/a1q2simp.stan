data {
  int<lower=0> N;
  vector[N] y;
  vector[2] gprior;
}

parameters {
  real<lower=0> lambda;
}


model {
  lambda ~ gamma(gprior[1], gprior[2]);
  y ~ exponential(lambda);
}

generated quantities {
  real postdraw = exponential_rng(lambda);
  real ewt = 1/lambda;
}
