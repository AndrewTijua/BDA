functions{
  real lindley_lpdf(real y, real lambda){
    real lpdf = 2*log(lambda) - log(lambda+1) + log(1+y) - lambda * y;
    return lpdf;
  }
  real lindley_rng(real lambda){
    real u = uniform_rng(0,1);
    real p = lambda/(lambda+1);
    if (u<=p){
      real v = exponential_rng(lambda);
      return v;
    }
    else{
    real w = gamma_rng(2, lambda);
    return w;
    }
  }
}
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
  for (n in 1:N){
    y[n] ~ lindley(lambda);
  }
  target += log_mix(theta[1],
  gamma_lpdf(lambda | alpha[1], beta[1]),
  gamma_lpdf(lambda | alpha[2], beta[2]));
}
generated quantities {
  real<lower=0> mwt =  ((lambda + 2)/(lambda*(lambda + 1)));
  real postdraw = lindley_rng(lambda);
}
