
######################################
#Calculate the conditional power
######################################

# Interim data
n_t_interim <- 78
n_c_interim <- 63
events_t_interim <- 31
events_c_interim <- 23

# Remaining patients
n_t_total <- 162
n_c_total <- 162
n_t_remaining <- n_t_total - n_t_interim  
n_c_remaining <- n_c_total - n_c_interim  

# Observed event rates
p_t_obs <- events_t_interim / n_t_interim  
p_c_obs <- events_c_interim / n_c_interim  

# Simulation settings
n_sim <- 10000000
alpha <- 0.025

CP_Func <- function(control_event_rate, treatment_event_rate){
  
  
  # Simulate future outcomes
  treatment_sim <- rbinom(n_sim, n_t_remaining, treatment_event_rate)
  control_sim <- rbinom(n_sim, n_c_remaining, control_event_rate)
  
  # Total events
  treatment_total <- treatment_sim + events_t_interim
  control_total <- control_sim + events_c_interim
  
  # Final sample sizes
  n_t_final <- n_t_interim + n_t_remaining
  n_c_final <- n_c_interim + n_c_remaining
  
  # Proportions
  p_t_final <- treatment_total / n_t_final
  p_c_final <- control_total / n_c_final
  
  p_pool <- (control_total + treatment_total) / (2 * n_t_total)
  se <- sqrt(2 * p_pool * (1 - p_pool) / n_t_total)
  
  
  # Z-statistics
  z <- (p_c_final - p_t_final)  / se
  
  conditional_power <- mean(z>qnorm(1 - 0.025))
  
  conditional_power
  # Output
  #cat("Conditional Power (simulated):", round(conditional_power * 100, 4), "%\n")
}



CP_Func(0.35, 0.35)

######################################
#Calculate the predictive power
######################################


alpha_C_prior <- 10.7
beta_C_prior <- 13.1


mu_d_prior <- 0.15
sigma_d_prior <- 1000 

N_C_interim <- 63 # Number of patients in control group at interim
Y_C_interim <- 23  # Number of events in control group at interim (e.g., 6/50 = 12%)

N_T_interim <- 78 # Number of patients in treatment group at interim
Y_T_interim <- 31  # Number of events in treatment group at interim (e.g., 3/50 = 6%)


# --- 3. Prepare data for Stan ---
stan_data <- list(
  N_C = N_C_interim,
  Y_C = Y_C_interim,
  N_T = N_T_interim,
  Y_T = Y_T_interim,
  alpha_C_prior = alpha_C_prior,
  beta_C_prior = beta_C_prior,
  mu_d_prior = mu_d_prior,
  sigma_d_prior = sigma_d_prior
)


stan_model <- stan_model("Thesis/Chapter 4/stan_model_binary_outcome.stan")


fit <- sampling(stan_model,
                data = stan_data,
                chains = 4,      # Number of MCMC chains
                iter = 2000,     # Total iterations per chain (includes warmup)
                warmup = 1000,   # Warmup iterations per chain
                thin = 1,        # Thinning rate (1 means no thinning)
                seed = 1234      # For reproducibility
)



print(fit, pars = c("p_C", "d"),
      probs = c(0.025, 0.5, 0.975)) # 2.5%, 50% (median), 97.5% quantiles for 95% CI


mcmc_dens(fit, pars = c("p_C", "d"),
          prob_outer = 0.95) + # Show 95% credible interval
  vline_at(0, linetype = 2, color = "red") # Add a vertical line at d=0 for the difference


# Extract posterior samples
posterior_samples <- as.data.frame(fit)



