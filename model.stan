data {
  int<lower=0> N;               // Number of observations
  vector[N] time;               // Time to event or censoring
  int<lower=0, upper=1> status[N];  // Event indicator
  int<lower=0, upper=1> group[N];   // Group (0=Control, 1=Treatment)
  real<lower=0> a_delay;
  real<lower=0> b_delay;
  real<lower=0> alpha_lambda;
  real<lower=0> beta_lambda;
  real mu_HR;
  real<lower=0> sigma_HR;
}
parameters {
  real<lower=0> lambda;
  real<lower=0> HR1;
  real<lower=0> HR2;
  real<lower=a_delay, upper=b_delay> delayT;
}
model {
  // Priors
  lambda ~ gamma(alpha_lambda, beta_lambda);
  HR1 ~ lognormal(mu_HR, sigma_HR);
  HR2 ~ lognormal(mu_HR, sigma_HR);
  delayT ~ uniform(a_delay, b_delay);
  
  // Likelihood
  for (i in 1:N) {
    real hazard;
    real cumulative_hazard;
    if (group[i] == 0) { // Control group
      hazard = lambda;
      cumulative_hazard = lambda * time[i];
    } else {             // Treatment group
      if (time[i] < delayT) {
        hazard = lambda * HR1;
        cumulative_hazard = lambda * HR1 * time[i];
      } else {
        hazard = lambda * HR2;
        cumulative_hazard = lambda * HR1 * delayT + lambda * HR2 * (time[i] - delayT);
      }
    }
    if (status[i] == 1) {
      target += log(hazard) - cumulative_hazard; // Event
    } else {
      target += -cumulative_hazard;             // Censor
    }
  }
}
