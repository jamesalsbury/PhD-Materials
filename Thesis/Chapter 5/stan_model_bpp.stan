data {
  // ---- Interim data ----
  int<lower=0> N_C;           // interim control n
  int<lower=0> Y_C;           // interim control events
  int<lower=0> N_T;           // interim treatment n
  int<lower=0> Y_T;           // interim treatment events

  // ---- Future (remaining) data ----
  int<lower=0> N_future_C;    // remaining control sample
  int<lower=0> N_future_T;    // remaining treatment sample

  // ---- Prior parameters ----
  real<lower=0> alpha_C_prior;
  real<lower=0> beta_C_prior;
  real mu_d_prior;
  real<lower=0> sigma_d_prior;
}

parameters {
  real<lower=0, upper=1> p_C;   // control event rate
  real d;                       // treatment effect (log-odds difference)
}

transformed parameters {
  // Treatment probability (positive d => beneficial treatment)
  real logit_pC = logit(p_C);
  real p_T = inv_logit(logit_pC - d);
}

model {
  // ---- Priors ----
  p_C ~ beta(alpha_C_prior, beta_C_prior);
  d ~ normal(mu_d_prior, sigma_d_prior);

  // ---- Likelihood for interim data ----
  Y_C ~ binomial(N_C, p_C);
  Y_T ~ binomial(N_T, p_T);
}

generated quantities {
  int Y_future_C = binomial_rng(N_future_C, p_C);
  int Y_future_T = binomial_rng(N_future_T, p_T);

  int total_Y_C = Y_C + Y_future_C;
  int total_Y_T = Y_T + Y_future_T;
  int total_N_C = N_C + N_future_C;
  int total_N_T = N_T + N_future_T;

  real p_c_final = total_Y_C / (1.0 * total_N_C);
real p_t_final = total_Y_T / (1.0 * total_N_T);
real p_pool = (total_Y_C + total_Y_T) / (1.0 * (total_N_C + total_N_T));


  real se;
  real z_stat;
  int success;
  real rho_final;

  // Compute SE safely
  se = sqrt(p_pool * (1 - p_pool) * (1.0 / total_N_C + 1.0 / total_N_T));

  // If SE == 0 (degenerate case), assign z_stat = 0 and success = 0
  if (se > 0) {
    z_stat = (p_c_final - p_t_final) / se;
    success = (z_stat > 1.96);
  } else {
    z_stat = 0;
    success = 0;
  }

  rho_final = p_c_final - p_t_final;
}
