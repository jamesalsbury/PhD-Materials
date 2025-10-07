
# Calculate the conditional power ----------------------------------------------------------------

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
n_sim <- 1000000
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
  
}


theta_f <- seq(-10, 30, by=1)/100
CP_Vec <- rep(NA, length(theta_f))

for (i in 1:length(CP_Vec)){
  CP_Vec[i] <- CP_Func(0.45, 0.45 - theta_f[i])
}

png("CP_Moxo.png", units="in", width=10, height=6, res=700)
plot(theta_f, CP_Vec, type = "l", lwd = 2, col = "steelblue",
     xlab = expression(theta[f]), 
     ylab = "Conditional Power", ylim = c(0,1),
     cex.axis=1.5, cex.lab=1.5, cex.main=2)

theta_obs <- p_c_obs - p_t_obs

# Observed effect
abline(v = theta_obs, col = "#8A2BE2", lty = 3, lwd = 2)  # a brighter purple
text(theta_obs, 0.6, expression(theta[f] == hat(theta)[t]), col = "#8A2BE2", pos =2)

# Null effect
abline(v = 0, col = "#D62728", lty = 2, lwd = 2)  # a softer red (less harsh)
text(0, 0.55, expression(theta[f] == 0), col = "#D62728", pos = 4)

# Clinically relevant effect
abline(v = 0.15, col = "#FF7F0E", lty = 1, lwd = 2)  # warm orange
text(0.15, 0.5, expression(theta[f] == theta[delta]), col = "#FF7F0E", pos = 4)


dev.off()


# Do conjugate updating for theta_c and theta_t ----------------------------------------------------------------

# Priors for theta_c and theta_t
theta_c_alpha_prior <- 10.7
theta_c_beta_prior  <- 13.1

theta_t_alpha_prior <- 6
theta_t_beta_prior  <- 14

#Interim data
c_x <- 23
c_n <- 63
t_x <- 31
t_n <- 78

# Posteriors for theta_c and theta_t
theta_c_alpha_post <- theta_c_alpha_prior + c_x
theta_c_beta_post  <- theta_c_beta_prior + c_n - c_x
theta_t_alpha_post <- theta_t_alpha_prior + t_x
theta_t_beta_post  <- theta_t_beta_prior + t_n - t_x


# Sequence of theta values
theta <- seq(0, 1, length.out = 1000)

# Prior densities
theta_c_prior <- dbeta(theta, theta_c_alpha_prior, theta_c_beta_prior)
theta_t_prior <- dbeta(theta, theta_t_alpha_prior, theta_t_beta_prior)

# Likelihoods (Binomial likelihood scaled for plotting)
# theta_c_likelihood <- dbinom(c_x, c_n, theta)
# theta_c_likelihood <- theta_c_likelihood / max(theta_c_likelihood) * max(theta_c_prior)
# theta_t_likelihood <- dbinom(t_x, t_n, theta)
# theta_t_likelihood <- theta_t_likelihood / max(theta_t_likelihood) * max(theta_t_prior) 

# Posterior densities
theta_c_posterior <- dbeta(theta, theta_c_alpha_post, theta_c_beta_post)
theta_t_posterior <- dbeta(theta, theta_t_alpha_post, theta_t_beta_post)


# Start high-quality PNG
png("Moxo_conj_dists.png", units="in", width=10, height=6, res=700)

# Set plotting range
ylim_max <- max(theta_c_prior, theta_t_prior, theta_c_posterior, theta_t_posterior)

# Plot priors and posteriors
plot(theta, theta_c_prior, type="l", lwd=2, col="blue",
     ylim=c(0, ylim_max),
     ylab="Density",
     xlab=expression(theta),
     cex.axis=1.5, cex.lab=1.5, cex.main=1.5)

lines(theta, theta_c_posterior, lwd=2, col="blue", lty=2)
lines(theta, theta_t_prior, lwd=2, col="red")
lines(theta, theta_t_posterior, lwd=2, col="red", lty=2)

legend("topright",
       legend=c(expression("Prior " * theta[c]),
                expression("Posterior " * theta[c]),
                expression("Prior " * theta[t]),
                expression("Posterior " * theta[t])),
       col=c("blue", "blue", "red", "red"),
       lwd=2, lty=c(1,2,1,2),
       bty="n", cex=1.2)

dev.off()

# Do conjugate updating for theta_t ----------------------------------------------------------------


# Parameters
alpha_prior_scen1<- 5.93
beta_prior_scen1  <- 13.87

alpha_prior_scen2 <- 2.95
beta_prior_scen2  <- 6.86

alpha_prior_scen3 <- 3.68
beta_prior_scen3  <- 6.83


x <- 31
n <- 78

# Posterior parameters
alpha_post_scen1 <- alpha_prior_scen1 + x
beta_post_scen1  <- beta_prior_scen1 + n - x

alpha_post_scen2 <- alpha_prior_scen2 + x
beta_post_scen2 <- beta_prior_scen2 + n - x

alpha_post_scen3 <- alpha_prior_scen3 + x
beta_post_scen3  <- beta_prior_scen3 + n - x

# Sequence of theta values
theta <- seq(0, 1, length.out = 1000)

# Prior density
prior_scen1 <- dbeta(theta, alpha_prior_scen1, beta_prior_scen1)
prior_scen2 <- dbeta(theta, alpha_prior_scen2, beta_prior_scen2)
prior_scen3 <- dbeta(theta, alpha_prior_scen3, beta_prior_scen3)



# Likelihood (Binomial likelihood scaled for plotting)
likelihood <- dbinom(x, n, theta)
likelihood <- likelihood / max(likelihood) * max(prior)  # scale to match prior for visualization

# Posterior density
posterior_scen1 <- dbeta(theta, alpha_post_scen1, beta_post_scen1)
posterior_scen2 <- dbeta(theta, alpha_post_scen2, beta_post_scen2)
posterior_scen3 <- dbeta(theta, alpha_post_scen3, beta_post_scen3)


png("Theta_t_Dens.png", units="in", width=10, height=6, res=700)

# Plot
plot(theta, prior_scen1, type="l", lwd=2, col="blue", ylim=c(0, max(prior_scen1, posterior_scen1)),
     ylab="Density", xlab=expression(theta[t]), cex.axis=1.5, cex.lab=1.5, cex.main=2)

lines(theta, posterior_scen1, lwd=2, col="blue", lty = 2)


lines(theta, prior_scen2, lwd=2, col="red", lty = 1)
lines(theta, posterior_scen2, lwd=2, col="red", lty = 2)

lines(theta, prior_scen3, lwd=2, col="green", lty = 1)
lines(theta, posterior_scen3, lwd=2, col="green", lty = 2)


legend("topright",
       legend=c("Scenario 1", "Scenario 2", "Scenario 3"),
       col=c("blue", "red", "green"),
       lwd=2,
       lty=c(1,1,1),
       title="Scenario",
       cex=1)

legend("right",
       legend=c("Prior", "Posterior"),
       col="black",
       lwd=2,
       lty=c(1,2),
       title="Line type",
       cex=1,
       inset=c(0,0))   # pushes to the side


dev.off()


# Calculate the predictive power - Conjugate updating ----------------------------------------------------------------

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

# Parameters
theta_c_alpha_prior <- 10.7
theta_c_beta_prior  <- 13.1
theta_t_alpha_prior <- 2.95
theta_t_beta_prior  <- 6.86


# Posterior parameters
theta_c_alpha_post <- theta_c_alpha_prior + events_c_interim
theta_c_beta_post  <- theta_c_beta_prior + n_c_interim - events_c_interim
theta_t_alpha_post <- theta_t_alpha_prior + events_t_interim
theta_t_beta_post  <- theta_t_beta_prior + n_t_interim - events_t_interim


n_sims <- 10000
ass_vec <- rep(NA, n_sims)

for (i in 1:n_sims){
  control_event_rate <- rbeta(1, theta_c_alpha_post, theta_c_beta_post)
  treatment_event_rate <- rbeta(1, theta_t_alpha_post, theta_t_beta_post)
  ass_vec[i] <- BPP_Func(control_event_rate = control_event_rate,
                         treatment_event_rate = treatment_event_rate)
}

mean(ass_vec)


BPP_Func <- function(control_event_rate, treatment_event_rate){
  
  
  # Simulate future outcomes
  treatment_sim <- rbinom(1, n_t_remaining, treatment_event_rate)
  control_sim <- rbinom(1, n_c_remaining, control_event_rate)
  
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
  
  return(z>qnorm(1 - 0.025))
  
}


# Calculate the predictive power - MCMC ----------------------------------------------------------------


# Priors + data
alpha_C_prior <- 10.7
beta_C_prior <- 13.1
mu_d_prior <- 0.15
sigma_d_prior <- 1000 

stan_data <- list(
  N_C = 63,
  Y_C = 23,
  N_T = 78,
  Y_T = 31,
  alpha_C_prior = alpha_C_prior,
  beta_C_prior = beta_C_prior,
  mu_d_prior = mu_d_prior,
  sigma_d_prior = sigma_d_prior
)

# Compile Stan model
mod <- stan_model("Thesis/Chapter 5/stan_model_binary_outcome.stan")

# Sample
fit <- sampling(mod,
                data = stan_data,
                chains = 4,
                iter = 2000,
                warmup = 1000,
                seed = 1234)

# Look at results
print(fit)



print(fit, pars = c("p_C", "d"),
      probs = c(0.025, 0.5, 0.975)) # 2.5%, 50% (median), 97.5% quantiles for 95% CI


mcmc_dens(fit, pars = c("p_C", "d"),
          prob_outer = 0.95) + # Show 95% credible interval
  vline_at(0, linetype = 2, color = "red") # Add a vertical line at d=0 for the difference


# Extract posterior samples
posterior_samples <- as.data.frame(fit)

# Plotting GSD boundaries - Pocock and OBF ----------------------------------------------------------------

Pocock_Design <- rpact::getDesignGroupSequential(kMax = 4,
                                                 alpha = 0.025,
                                                 sided = 1,
                                                 typeOfDesign = "P",
                                                 informationRates = c(0.25, 0.5, 0.75, 1))

OBF_Design <- rpact::getDesignGroupSequential(kMax = 4,
                                              alpha = 0.025,
                                              sided = 1,
                                              typeOfDesign = "OF",
                                              informationRates = c(0.25, 0.5, 0.75, 1))

# Information fractions
IF <- seq(0.25, 1, by = 0.25)


png("Pocock_OBF_Boundaries.png", units="in", width=10, height=6, res=700)
# Plot Pocock boundaries without default x-axis
plot(IF, Pocock_Design$criticalValues, type = "b", pch = 16, lty = 1, col = "blue",
     ylim = c(0, 5), xlim = c(0, 1),
     ylab = "Standardized Z Statistic", xlab = "Information Fraction",
     xaxt = "n")

# Add custom x-axis ticks at 0, 0.25, 0.5, 0.75, 1
axis(1, at = seq(0, 1, by = 0.25), labels = seq(0, 1, by = 0.25))

# Add O'Brien–Fleming boundaries
lines(IF, OBF_Design$criticalValues, type = "b", pch = 17, lty = 2, col = "red")

# Add grid for readability
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")

# Add legend
legend("topright", legend = c("Pocock", "O'Brien–Fleming"),
       col = c("blue", "red"), lty = c(1, 2), pch = c(16, 17), bty = "n")
dev.off()


# Plotting GSD boundaries - Wang-Tsiatis ----------------------------------------------------------------

WT_Design_1 <- rpact::getDesignGroupSequential(kMax = 4, alpha = 0.025, sided = 1, typeOfDesign = "WT", informationRates = c(0.25, 0.5, 0.75, 1), deltaWT = 0) 
WT_Design_2 <- rpact::getDesignGroupSequential(kMax = 4, alpha = 0.025, sided = 1, typeOfDesign = "WT", informationRates = c(0.25, 0.5, 0.75, 1), deltaWT = 0.25) 
WT_Design_3 <- rpact::getDesignGroupSequential(kMax = 4, alpha = 0.025, sided = 1, typeOfDesign = "WT", informationRates = c(0.25, 0.5, 0.75, 1), deltaWT = 0.5)

# Information fractions
IF <- seq(0.25, 1, by = 0.25)

png("WT_Boundaries.png", units="in", width=10, height=6, res=700)

# Plot WT design (delta = 0)
plot(IF, WT_Design_1$criticalValues, type = "b", pch = 16, lty = 1, col = "blue",
     ylim = c(0, 5), xlim = c(0, 1),
     ylab = "Standardized Z Statistic", xlab = "Information Fraction",
     xaxt = "n")

# Add custom x-axis ticks at 0, 0.25, 0.5, 0.75, 1
axis(1, at = seq(0, 1, by = 0.25), labels = seq(0, 1, by = 0.25))

# Add other WT designs with different delta
lines(IF, WT_Design_2$criticalValues, type = "b", pch = 17, lty = 2, col = "red")
lines(IF, WT_Design_3$criticalValues, type = "b", pch = 15, lty = 3, col = "darkgreen")

# Add grid for readability
grid(nx = NA, ny = NULL, col = "gray80", lty = "dotted")

# Add legend with correct labels
legend("topright",
       legend = c(expression("WT, " * Delta * " = 0"),
                  expression("WT, " * Delta * " = 0.25"),
                  expression("WT, " * Delta * " = 0.5")),
       col = c("blue", "red", "darkgreen"),
       lty = c(1, 2, 3),
       pch = c(16, 17, 15),
       bty = "n")

dev.off()

# GSD Example ----------------------------------------------------------------


Pocock_Design <- gsDesign(
  k = 2,                     # Number of analyses
  alpha = 0.025,             # One-sided Type I error
  beta = 0.1,   # Type II error
  delta = 0.5,             # Standardized effect size
  sfu = "Pocock", 
  test.type = 1,             # One-sided test
  n.fix = 85*2,  # Fix max sample size per group
  timing = c(0.5, 1)
)


OBF_Design <- gsDesign(
  k = 2,                     # Number of analyses
  alpha = 0.025,             # One-sided Type I error
  beta = 0.1,   # Type II error
  delta = 0.5,             # Standardized effect size
  sfu = "OF", 
  test.type = 1,             # One-sided test
  n.fix = 85*2,  # Fix max sample size per group
  timing = c(0.5, 1)
)

WT_Design <- gsDesign(
  k = 2,                     # Number of analyses
  alpha = 0.025,             # One-sided Type I error
  beta = 0.1,   # Type II error
  delta = 0.5,             # Standardized effect size
  sfu = "WT",
  sfupar = 0.25,
  test.type = 1,             # One-sided test
  n.fix = 85*2,  # Fix max sample size per group
  timing = c(0.5, 1)
)

simulate_gsd_power <- function(
    gsd,
    mu_control = 120,
    mu_treat = 125,
    sigma = 10,
    n_per_group_fixed = 85,
    n_sim = 200000
) {
  k <- gsd$k
  bounds <- gsd$upper$bound
  
  # Sample sizes per group at each look
  n_looks <- c(round(n_per_group_fixed * 0.5), n_per_group_fixed)  # e.g. 43, 85
  total_n <- n_per_group_fixed
  
  success <- 0
  sample_sizes <- numeric(n_sim)
  stop_counts <- integer(k)
  
  for (i in 1:n_sim) {
    # Generate interim data
    n1 <- n_looks[1]
    x_control_1 <- rnorm(n1, mean = mu_control, sd = sigma)
    x_treat_1 <- rnorm(n1, mean = mu_treat, sd = sigma)
    
    # Interim analysis
    s2_control_1 <- var(x_control_1)
    s2_treat_1 <- var(x_treat_1)
    pooled_var_1 <- ((n1 - 1) * s2_control_1 + (n1 - 1) * s2_treat_1) / (2 * n1 - 2)
    se_1 <- sqrt(pooled_var_1 * (2 / n1))
    z_1 <- (mean(x_treat_1) - mean(x_control_1)) / se_1
    
    if (z_1 > bounds[1]) {
      success <- success + 1
      sample_sizes[i] <- 2 * n1
      stop_counts[1] <- stop_counts[1] + 1
      next
    }
    
    # Generate additional data for final analysis
    n2 <- n_looks[2] - n1
    x_control_2 <- rnorm(n2, mean = mu_control, sd = sigma)
    x_treat_2 <- rnorm(n2, mean = mu_treat, sd = sigma)
    
    # Combine data
    x_control <- c(x_control_1, x_control_2)
    x_treat <- c(x_treat_1, x_treat_2)
    
    s2_control <- var(x_control)
    s2_treat <- var(x_treat)
    pooled_var <- ((total_n - 1) * s2_control + (total_n - 1) * s2_treat) / (2 * total_n - 2)
    se <- sqrt(pooled_var * (2 / total_n))
    z <- (mean(x_treat) - mean(x_control)) / se
    
    if (z > bounds[2]) {
      success <- success + 1
    }
    sample_sizes[i] <- 2 * total_n
    stop_counts[2] <- stop_counts[2] + 1
  }
  
  list(
    power = success / n_sim,
    expected_sample_size = mean(sample_sizes),
    max_sample_size = 2 * total_n,
    stop_distribution = stop_counts / n_sim
  )
}


results <- simulate_gsd_power(gsd)
print(results)


