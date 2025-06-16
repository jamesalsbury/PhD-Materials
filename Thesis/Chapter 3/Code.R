library(survival)
library(rjags)

IPD_INTEREST_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_INTEREST_OS.csv")
IPD_ZODIAC_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_ZODIAC_OS.csv")
IPD_REVEL_OS <- read.csv(file = "Thesis/Chapter 3/Data/IPD_REVEL_OS.csv")

#png("Doce_Surv.png", units="in", width=10, height=8, res=700)
INTEREST_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_INTEREST_OS)
plot(INTEREST_kmfit, xlim = c(0, 21), col = "blue", conf.int = F, ylab = "Overall Survival", xlab = "Time (months)")
INTEREST_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_INTEREST_OS, dist = "exponential")
INTEREST_lambda <- exp(-INTEREST_myfit$coefficients)
time_seq <- seq(0, 21, by=0.01)
surv_probs <- exp(-INTEREST_lambda*time_seq)
#lines(time_seq, surv_probs, col = "blue")


ZODIAC_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_ZODIAC_OS)
lines(ZODIAC_kmfit, xlim = c(0, 21), col = "red", conf.int = F)
ZODIAC_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_ZODIAC_OS, dist = "exponential")
ZODIAC_lambda <- exp(-ZODIAC_myfit$coefficients)
surv_probs <- exp(-ZODIAC_lambda*time_seq)
#lines(time_seq, surv_probs, col = "red")

REVEL_kmfit <- survfit(Surv(Survival.time, Status)~1, data = IPD_REVEL_OS)
lines(REVEL_kmfit, xlim = c(0, 21), col = "green", conf.int = F)
REVEL_myfit <- survreg(Surv(Survival.time, Status)~1, data = IPD_REVEL_OS, dist = "exponential")
REVEL_lambda <- exp(-REVEL_myfit$coefficients)
surv_probs <- exp(-REVEL_lambda*time_seq)
#lines(time_seq, surv_probs, col= "green")

legend("topright", legend = c("INTEREST", "ZODIAC", "REVEL"), col = c("blue", "red", "green"), lty = 1)

med_INTEREST <- summary(INTEREST_kmfit)$table["median"]
med_ZODIAC <- summary(ZODIAC_kmfit)$table["median"]
med_REVEL <- summary(REVEL_kmfit)$table["median"]



lam <- c(log(2)/med_INTEREST, log(2)/med_ZODIAC, log(2)/med_REVEL) 
dc <- c(sum(IPD_INTEREST_OS$Status==1), sum(IPD_ZODIAC_OS$Status==1), sum(IPD_REVEL_OS$Status==1))
hist_data <- list(H = length(lam),                  
                dc = dc,                
                Tc = dc/lam,            
                prior_prec_tau = 1)     

# function to generate inital values for MCMC
hist_init <- function() {
  list(log_lam = rnorm(hist_data$H, mean = 0, sd = 0.5),
       mu = rnorm(1, mean = 1, sd = 0.5),
       tau = rexp(1, rate = 3)  
  )  
}

# Call WinBUGS
model <- jags.model("Thesis/Chapter 3/jags_model.txt",
                    data = hist_data,
                    inits = hist_init,
                    n.chains = 5,
                    n.adapt = 1000)

update(model, 1000) 

params <- c("lam_pred")


samples <- coda.samples(model,
                        variable.names = params,
                        n.iter = 400000,
                        thin = 20)

library(coda)
combined_samples <- do.call(rbind, samples)
lam_pred_all <- as.numeric(combined_samples[,1])


# Plot histogram
hist(lam_pred_all,
     breaks = 5000,
     main = "",
     xlim = c(0, 0.2),
     freq = F,
     ylim = c(0, 25),
     xlab = expression("Predicted "*lambda[c])
)

plot(INTEREST_kmfit, xlim = c(0, 21), col = "blue", conf.int = F, ylab = "Overall Survival", xlab = "Time (months)")
lines(ZODIAC_kmfit, xlim = c(0, 21), col = "red", conf.int = F)
lines(REVEL_kmfit, xlim = c(0, 21), col = "green", conf.int = F)
Med <- exp(-quantile(lam_pred_all, 0.5)*time_seq)
LQ <- exp(-quantile(lam_pred_all, 0.05)*time_seq)
UQ <- exp(-quantile(lam_pred_all, 0.95)*time_seq)
lines(time_seq, Med)
lines(time_seq, LQ, lty = 2)
lines(time_seq, UQ, lty = 2)

legend("topright", legend = c("INTEREST", "ZODIAC", "REVEL", "MAP Prior", "90% Interval for MAP Prior"), 
       col = c("blue", "red", "green", "black", "black"), lty = c(1, 1, 1, 1, 2))


MCMC_sample <- read.csv(file = "Thesis/Chapter 3/Data/MCMC_sample.csv")
Assurance_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                             n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                             control_dist = "Exponential",
                             control_parameters = "Distribution",
                             control_distribution = "MCMC",
                             MCMC_sample = MCMC_sample,
                             delay_time_SHELF = SHELF::fitdist(c(3, 4, 5), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                             delay_time_dist = "gamma",
                             post_delay_HR_SHELF = SHELF::fitdist(c(0.55, 0.6, 0.7), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                             post_delay_HR_dist = "gamma",  
                             P_S = 0.9,
                             P_DTE = 0.8,
                             cens_method = "Events",
                             cens_IF = 0.8,
                             rec_method = "power",
                             rec_period = 12,
                             rec_power = 1,
                             analysis_method = "LRT",
                             alternative = "one.sided",
                             alpha = 0.025,
                             nSims = 25)

Power_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              control_dist = "Exponential",
                              control_parameters = "Fixed",
                              fixed_parameters_type = "Parameter",
                              lambda_c = 0.077,
                              delay_time_SHELF = SHELF::fitdist(c(3.999, 4, 4.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                              delay_time_dist = "gamma",
                              post_delay_HR_SHELF = SHELF::fitdist(c(0.599, 0.6, 0.601), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                              post_delay_HR_dist = "gamma",  
                              P_S = 1,
                              P_DTE = 1,
                              cens_method = "Events",
                              cens_IF = 0.8,
                              rec_method = "power",
                              rec_period = 12,
                              rec_power = 1,
                              analysis_method = "LRT",
                              alternative = "one.sided",
                              alpha = 0.025,
                              nSims = 25)


Power_ND_Calc <- calc_dte_assurance(n_c = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              n_t = c(5, 60, 115, 170, 225, 280, 335, 390, 445, 500),
                              control_dist = "Exponential",
                              control_parameters = "Fixed",
                              fixed_parameters_type = "Parameter",
                              lambda_c = 0.077,
                              delay_time_SHELF = SHELF::fitdist(c(3.999, 4, 4.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                              delay_time_dist = "gamma",
                              post_delay_HR_SHELF = SHELF::fitdist(c(0.599, 0.6, 0.601), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                              post_delay_HR_dist = "gamma",  
                              P_S = 1,
                              P_DTE = 0,
                              cens_method = "Events",
                              cens_IF = 0.8,
                              rec_method = "power",
                              rec_period = 12,
                              rec_power = 1,
                              analysis_method = "LRT",
                              alternative = "one.sided",
                              alpha = 0.025,
                              nSims = 25)

outcome_DF <- data.frame(N = seq(10, 1000, length = 10),
                          Ass = sapply(Assurance_Calc, `[[`, 1),
                          Power = sapply(Power_Calc, `[[`, 1),
                          Power_ND = sapply(Power_ND_Calc, `[[`, 1))

png("PowerAss_Exp_Example.png", units="in", width=10, height=8, res=700)
plot(outcome_DF$N, outcome_DF$Ass, type= "l", ylim = c(0,1),
     xlab = "Number of Patients", ylab = "Assurance/Power")
lines(outcome_DF$N, outcome_DF$Power, type= "l", col = "blue", lty = 2)
lines(outcome_DF$N, outcome_DF$Power_ND, type= "l", col = "red", lty = 3)

legend("bottomright", legend = c("Assurance", "Power", "Power assuming no delay"),
       col = c("black", "blue", "red"), lty = 1:3)
dev.off()





