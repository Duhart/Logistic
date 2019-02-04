// Bayesian logistic regresion
// file saved as bayes_log_reg.stan
data {
  int<lower=0> N;             // number of data points 
  int<lower=0,upper=1> y[N];  // value 0 or 1
  vector[N] x1;                 // First covariate
  vector[N] x2;                 // Second covariate
  real b0_init;
  real b1_init;
  real b2_init;
  real<lower=0> sig;
}
parameters {
  real b0;                    // intercept
  real b1;                    // slope of x1
  real b2;                    // slope of x2
}
transformed parameters {
  vector[N] p = 1 ./ (1 + exp(-b0 - x1 * b1 - x2 * b2)); // probabilities
}
model {
  b0 ~ normal(b0_init,sig);
  b1 ~ normal(b1_init,sig);
  b2 ~ normal(b2_init,sig);
  for(k in 1:N)
      y[k] ~ bernoulli(p[k]);
}
