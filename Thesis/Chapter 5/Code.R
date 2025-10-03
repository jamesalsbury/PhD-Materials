
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
  # Output
  #cat("Conditional Power (simulated):", round(conditional_power * 100, 4), "%\n")
}


theta_f <- seq(-5, 30, by=1)/100
CP_Vec <- rep(NA, length(theta_f))

for (i in 1:length(CP_Vec)){
  CP_Vec[i] <- CP_Func(0.45, 0.45 - theta_f[i])
}

png("CP_Moxo.png", units="in", width=10, height=6, res=700)
plot(theta_f, CP_Vec, type = "l", 
     xlab = expression(theta[f]), 
     ylab= "Conditional Power", ylim = c(0,1))
dev.off()

xSeq <- seq(0, 1, by=0.01)
plot(xSeq, dbeta(xSeq, 10.6875, 10.6875*11/9), type = "l", ylab = "Density", col = "blue", lwd = 2)
lines(xSeq, dbeta(xSeq, 6, 14), col = "red", lwd = 2)
legend("topright", legend = c(expression(theta[c] %~% Be(10.7, 13)),
                              expression(theta[t] %~% Be(12, 15))), col = c("blue", "red"), lty = 1, lwd = 2)


######################################
#Do conjugate updating for theta_c
######################################

# Parameters
alpha_prior <- 10.7
beta_prior  <- 13.1

x <- 23
n <- 63

# Posterior parameters
alpha_post <- alpha_prior + x
beta_post  <- beta_prior + n - x

# Sequence of theta values
theta <- seq(0, 1, length.out = 1000)

# Prior density
prior <- dbeta(theta, alpha_prior, beta_prior)

# Likelihood (Binomial likelihood scaled for plotting)
likelihood <- dbinom(x, n, theta)
likelihood <- likelihood / max(likelihood) * max(prior)  # scale to match prior for visualization

# Posterior density
posterior <- dbeta(theta, alpha_post, beta_post)

png("Theta_c_Dens.png", units="in", width=10, height=6, res=700)

# Plot
plot(theta, prior, type="l", lwd=2, col="blue", ylim=c(0, max(prior, posterior)),
     ylab="Density", xlab=expression(theta[c]), cex.axis=1.5, cex.lab=1.5, cex.main=2)

lines(theta, likelihood, lwd=2, col="green", lty=2)
lines(theta, posterior, lwd=2, col="red")

legend("topright", legend=c("Prior", "Likelihood", "Posterior"),
       col=c("blue","green","red"), lwd=2, lty=c(1,2,1))


dev.off()


######################################
#Do conjugate updating for theta_t
######################################

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


stan_model <- stan_model("Thesis/Chapter 5/stan_model_binary_outcome.stan")


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


######################################
#Plotting GSD boundaries - Pocock and OBF
######################################

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


######################################
#Plotting GSD boundaries - Wang-Tsiatis
######################################

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

######################################
#GSD Example
######################################

Pocock_Design <- rpact::getDesignGroupSequential(kMax = 2,
                                                 informationRates = c(0.5, 1),
                                                 typeOfDesign = "P", alpha = 0.025,
                                                 sided = 1, beta = 0.1)


OF_Design <- rpact::getDesignGroupSequential(kMax = 2,
                                             informationRates = c(0.5, 1),
                                             typeOfDesign = "OF", alpha = 0.025,
                                             sided = 1, beta = 0.1)

WT_Design <- rpact::getDesignGroupSequential(kMax = 2,
                                informationRates = c(0.5, 1),
                                typeOfDesign = "WT", alpha = 0.025,
                                sided = 1, deltaWT = 0.25, beta = 0.1)

# ESS_Pocock <- rpact::getSampleSizeMeans(
#   design = Pocock_Design,
#   groups = 2,
#   alternative = 5,
#   stDev = 10
# )
# 
# 
# ESS_OF <- rpact::getSampleSizeMeans(
#   design = OF_Design,
#   groups = 2,
#   alternative = 5,
#   stDev = 10
# )
# 
# 
# ESS_WT <- rpact::getSampleSizeMeans(
#   design = WT_Design,
#   groups = 2,
#   alternative = 5,
#   stDev = 10
# )


