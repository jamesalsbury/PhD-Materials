// stan_model_binary_outcome.stan
data {
  int<lower=0> N_C; // Number of patients in control group
  int<lower=0> Y_C; // Number of events in control group
  int<lower=0> N_T; // Number of patients in treatment group
  int<lower=0> Y_T; // Number of events in treatment group

  // Prior parameters for p_C (Beta distribution)
  real<lower=0> alpha_C_prior;
  real<lower=0> beta_C_prior;

  // Prior parameters for d (Normal distribution)
  real mu_d_prior;
  real<lower=0> sigma_d_prior;
}

parameters {
  real<lower=0, upper=1> p_C; // Event rate in control group
  real d_raw; // Unconstrained difference for Normal prior (will be transformed to d)
}

transformed parameters {
  real d = d_raw; // The actual difference. Could add transformations here if needed.
  real p_T = p_C + d; // Proposed treatment rate
  
  // No need for explicit 'target +=' here.
  // The binomial likelihood in the model block will implicitly
  // assign -infinity to the log-likelihood if p_T is out of [0,1].
  // However, for safety and clarity, it's often good practice to
  // ensure p_T is within bounds *before* it's used in the likelihood.
  // One common way is to "clip" it for deterministic parameters if not directly used
  // in a distribution that handles the bounds gracefully.
  // For parameters used in likelihoods, Stan's distributions handle this.
  // But let's add a soft constraint for robustness, if necessary, or rely on likelihood.
  // Given the direct use in binomial, the likelihood handles it.
}

model {
  // Priors
  p_C ~ beta(alpha_C_prior, beta_C_prior);
  d_raw ~ normal(mu_d_prior, sigma_d_prior);

  // Likelihood
  // Stan's binomial distribution will automatically handle cases where p_T
  // falls outside [0,1] by assigning a log-probability of -infinity,
  // effectively rejecting those parameter combinations.
  Y_C ~ binomial(N_C, p_C);
  Y_T ~ binomial(N_T, p_T);
}

generated quantities {
  // You can generate any quantities of interest here, e.g., relative risk, odds ratio
  real RR; // Relative Risk
  real OR; // Odds Ratio

  // It's good practice to ensure p_T and p_C are within bounds for these calculations
  // before division to prevent NaNs or Infs, though the model should mostly prevent this.
  if (p_C > 0 && p_T > 0 && p_T < 1) { // Added p_T < 1 for OR denominator
    RR = p_T / p_C;
    OR = (p_T / (1 - p_T)) / (p_C / (1 - p_C));
  } else {
    // Assign a sensible default or indicator for undefined cases
    RR = -1; // Indicate a problematic value
    OR = -1; // Indicate a problematic value
  }
}