data {
  int<lower=0> N_C;
  int<lower=0> Y_C;
  int<lower=0> N_T;
  int<lower=0> Y_T;
  real<lower=0> alpha_C_prior;
  real<lower=0> beta_C_prior;
  real mu_d_prior;
  real<lower=0> sigma_d_prior;
}

parameters {
  real<lower=0, upper=1> p_C;
  real d;
}

transformed parameters {
  real logit_pC = logit(p_C);
  real p_T = inv_logit(logit_pC - d);
}


model {
  // Priors
  p_C ~ beta(alpha_C_prior, beta_C_prior);
  d ~ normal(mu_d_prior, sigma_d_prior);

  // Likelihoods
  Y_C ~ binomial(N_C, p_C);
  Y_T ~ binomial(N_T, p_T);
}
