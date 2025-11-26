library(dplyr)
library(survival)
library(rjags)


source("Thesis/Chapter 6/functions.R")


# Reproduce the results from the paper -------------------------

n_c <- 340
n_t <- 340
num_events <- 512

scen_list <- list(A = list(HR1 = 0.75,
                           HR2 = 0.75,
                           delay = 0,
                           lambda_c = log(2)/12),
                  B = list(HR1 = 1,
                           HR2 = 1,
                           delay = 0,
                           lambda_c = log(2)/12),
                  C = list(HR1 = 1.3,
                           HR2 = 1.3,
                           delay = 0,
                           lambda_c = log(2)/12),
                  D = list(HR1 = 1,
                           HR2 = 0.693,
                           delay = 3,
                           lambda_c = log(2)/12),
                  E = list(HR1 = 1,
                           HR2 = 0.620,
                           delay = 6,
                           lambda_c = log(2)/12),
                  F = list(HR1 = 1.3,
                           HR2 = 0.628,
                           delay = 3,
                           lambda_c = log(2)/12),
                  G = list(HR1 = 1.3,
                           HR2 = 0.65, 
                           delay = 6,
                           lambda_c = log(2)/18),
                  H = list(HR1 = 1.3,
                           HR2 = 1.3, 
                           delay = 0,
                           lambda_c = log(2)/6),
                  I = list(HR1 = 0.7,
                           HR2 = 1.1,
                           delay = 6,
                           lambda_c = log(2)/15))



outcome_df <- data.frame(matrix(nrow = 12, ncol = 12))

for (i in 1:6){
  result <- reproduce_func(i, rec_duration = 34, n_sims = 100)
  
  outcome_df[i,1] <- mean(sapply(result$no_IA_list, function(x) x$power))
  outcome_df[i,2] <- mean(sapply(result$no_IA_list, function(x) x$duration))
  outcome_df[i,3] <- mean(sapply(result$no_IA_list, function(x) x$ss))
  
  
  outcome_df[i,4] <- mean(sapply(result$wieand_list, function(x) x$power))
  outcome_df[i,5] <- mean(sapply(result$wieand_list, function(x) x$duration))
  outcome_df[i,6] <- mean(sapply(result$wieand_list, function(x) x$ss))
  
  
  outcome_df[i,7] <- mean(sapply(result$OBF_list, function(x) x$power))
  outcome_df[i,8] <- mean(sapply(result$OBF_list, function(x) x$duration))
  outcome_df[i,9] <- mean(sapply(result$OBF_list, function(x) x$ss))
  
  
  outcome_df[i,10] <- mean(sapply(result$proposed_list, function(x) x$power))
  outcome_df[i,11] <-  mean(sapply(result$proposed_list, function(x) x$duration))
  outcome_df[i,12] <-  mean(sapply(result$proposed_list, function(x) x$ss))

  
}

for (i in 1:6){
  result <- reproduce_func(i, rec_duration = 12, n_sims = 100)
  
  outcome_df[i+6,1] <- mean(sapply(result$no_IA_list, function(x) x$power))
  outcome_df[i+6,2] <- mean(sapply(result$no_IA_list, function(x) x$duration))
  outcome_df[i+6,3] <- mean(sapply(result$no_IA_list, function(x) x$ss))
  
  outcome_df[i+6,4] <- mean(sapply(result$wieand_list, function(x) x$power))
  outcome_df[i+6,5] <- mean(sapply(result$wieand_list, function(x) x$duration))
  outcome_df[i+6,6] <- mean(sapply(result$wieand_list, function(x) x$ss))

  outcome_df[i+6,7] <- mean(sapply(result$OBF_list, function(x) x$power))
  outcome_df[i+6,8] <- mean(sapply(result$OBF_list, function(x) x$duration))
  outcome_df[i+6,9] <- mean(sapply(result$OBF_list, function(x) x$ss))
  
  outcome_df[i+6,10] <- mean(sapply(result$proposed_list, function(x) x$power))
  outcome_df[i+6,11] <-  mean(sapply(result$proposed_list, function(x) x$duration))
  outcome_df[i+6,12] <-  mean(sapply(result$proposed_list, function(x) x$ss))
  
}


colnames(outcome_df) <- rep(c("Power", "Duration", "SS"), 4) 
outcome_df[,1] <- round(outcome_df[,1], 3)
outcome_df[,2] <- round(outcome_df[,2], 1)
outcome_df[,3] <- round(outcome_df[,3], 1)
outcome_df[,4] <- round(outcome_df[,4], 3)
outcome_df[,5] <- round(outcome_df[,5], 1)
outcome_df[,6] <- round(outcome_df[,6], 1)
outcome_df[,7] <- round(outcome_df[,7], 3)
outcome_df[,8] <- round(outcome_df[,8], 1)
outcome_df[,9] <- round(outcome_df[,9], 1)
outcome_df[,10] <- round(outcome_df[,10], 3)
outcome_df[,11] <- round(outcome_df[,11], 1)
outcome_df[,12] <- round(outcome_df[,12], 1)

paper_df <- data.frame(
  Power = c(0.900, 0.0247, 0, 0.901, 0.903, 0.899, 0.899, 0.0246, 0, 0.901, 0.903, 0.899),
  Duration = c(47.4, 43.9, 41.3, 47.8, 48.3, 48.3, 34.3, 30.5, 27.7, 34.8, 35.4, 35.4),
  SS = c(680, 680, 680, 680, 680, 680, 680, 680, 680, 680, 680, 680),
  Power = c(0.898, 0.0245, 0, 0.899, 0.881, 0.872, 0.897, 0.0243, 0, 0.882, 0.804, 0.781),
  Duration = c(47.1, 34.1, 25.2, 47.3, 47.2, 47.0, 34.1, 21.5, 13.7, 33.9, 32.1, 31.6),
  SS = c(678.6, 602.6, 504.4, 676.7, 673.1, 671.7, 680, 680, 680, 680, 680, 680),
  Power = c(0.879, 0.0226, 0, 0.846, 0.786, 0.762, 0.879, 0.0223, 0, 0.794, 0.626, 0.537),
  Duration = c(46.2, 28.2, 19.9, 45.4, 43.9, 43.2, 33.3, 16.9, 10.9, 31.1, 27.1, 24.7),
  SS = c(671.7, 529.1, 399.3, 661.2, 644.3, 638.0, 679.7, 660.5, 607.5, 677.3, 670.4, 665.1),
  Power = c(0.898, 0.0245, 0, 0.895, 0.885, 0.891, 0.898, 0.0245, 0, 0.898, 0.879, 0.895),
  Duration = c(47.1, 34.3, 26.1, 47.4, 47.4, 47.8, 34.2, 23.3, 19.7, 34.5, 34.3, 35.1),
  SS = c(678.6, 605.3, 526.1, 677.4, 674.4, 677.8, 680, 680, 680, 680, 680, 680)
)

paper_df - outcome_df


# Plotting number of events v proportion of events > 3 months

#png("Events_Prop_Scen1.png", units="in", width=12, height=5, res=700)
#par(mfrow=c(1,2))
#events_v_prop_func(1, 34)
#events_v_prop_func(1, 12)
#dev.off()

# Find the situations which 'break' the rule -------------------------

break_df <- data.frame(matrix(nrow = 3, ncol = 12))

rec_vec <- c(0, 34, 12)


for (i in 1:3){
  
  result <- reproduce_func(i+6, rec_duration = rec_vec[i], n_sims = 100)

  
  break_df[i,1] <- mean(sapply(result$no_IA_list, function(x) x$power))
  break_df[i,2] <- mean(sapply(result$no_IA_list, function(x) x$duration))
  break_df[i,3] <- mean(sapply(result$no_IA_list, function(x) x$ss))
  
  
  break_df[i,4] <- mean(sapply(result$wieand_list, function(x) x$power))
  break_df[i,5] <- mean(sapply(result$wieand_list, function(x) x$duration))
  break_df[i,6] <- mean(sapply(result$wieand_list, function(x) x$ss))
  
  
  break_df[i,7] <- mean(sapply(result$OBF_list, function(x) x$power))
  break_df[i,8] <- mean(sapply(result$OBF_list, function(x) x$duration))
  break_df[i,9] <- mean(sapply(result$OBF_list, function(x) x$ss))
  
  
  break_df[i,10] <- mean(sapply(result$proposed_list, function(x) x$power))
  break_df[i,11] <-  mean(sapply(result$proposed_list, function(x) x$duration))
  break_df[i,12] <-  mean(sapply(result$proposed_list, function(x) x$ss))
  
}


colnames(break_df) <- rep(c("Power", "Duration", "SS"), 4) 
break_df[,1] <- round(break_df[,1], 3)
break_df[,2] <- round(break_df[,2], 1)
break_df[,3] <- round(break_df[,3], 1)
break_df[,4] <- round(break_df[,4], 3)
break_df[,5] <- round(break_df[,5], 1)
break_df[,6] <- round(break_df[,6], 1)
break_df[,7] <- round(break_df[,7], 3)
break_df[,8] <- round(break_df[,8], 1)
break_df[,9] <- round(break_df[,9], 1)
break_df[,10] <- round(break_df[,10], 3)
break_df[,11] <- round(break_df[,11], 1)
break_df[,12] <- round(break_df[,12], 1)

for (i in 1:3){
  
  plot_break_scenarios(i+6, rec_vec[i])

}


# Bayesian Predictive Power -------------------------


library(survival)
library(dplyr)
library(rjags)

num_patients <- 680

df <- sim_dte(n_c = num_patients/2, n_t = num_patients/2, lambda_c = log(2)/12,
              delay_time = 3, post_delay_HR = 0.75)

df$rec_time <- runif(num_patients, 0, 34)
df$pseudo_time <- df$time + df$rec_time

df <- cens_data(df, cens_method = "Events", cens_events = 256)

cens_time <- df$cens_time

df <- df$data

km_fit <- survfit(Surv(survival_time, status)~group, data = df)

png("KM_Plot_BPP_DTE.png", units="in", width=10, height=8, res=700)
plot(km_fit, col = c("blue", "red"), xlab = "Time", ylab = "Survival", cex.axis=1.5, cex.lab=1.5, cex.main=2)
legend("topright", legend = c("Control", "Treatment"), col = c("blue", "red"), lty = 1)
dev.off()    

#lambda_c_SHELF = SHELF::fitdist(c(log(2)/14, log(2)/12, log(2)/10), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10)

s1_SHELF = SHELF::fitdist(c(0.4, 0.5, 0.6), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1)
delay_SHELF = SHELF::fitdist(c(2, 3, 4), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10)
HR_SHELF = SHELF::fitdist(c(0.7, 0.75, 0.8), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 5)


posterior_df <- update_priors(data = df,
                              control_distribution = "Exponential",
                              control_model = list(s1_SHELF = s1_SHELF,
                                                   parameter_mode = "Landmark",
                                                   s1_dist = "Beta",
                                                   t_1 = 12),
                              delay_SHELF = delay_SHELF,
                              HR_SHELF = HR_SHELF,
                              delay_param_dist = "Gamma",
                              HR_param_dist = "Gamma"
)


prior_HR <- rgamma(50000, as.numeric(HR_SHELF$Gamma[1]), as.numeric(HR_SHELF$Gamma[2]))
prior_delay <- rgamma(50000, as.numeric(delay_SHELF$Gamma[1]), as.numeric(delay_SHELF$Gamma[2]))
prior_s1 <- rbeta(50000, as.numeric(s1_SHELF$Beta[1]), as.numeric(s1_SHELF$Beta[2]))
prior_lambda_c <- -log(prior_s1)/12



png("Prior_Post_BPP_DTE.png", units="in", width=15, height=8, res=700)
par(mfrow=c(1,3))

# Lambda_c
plot(density(prior_lambda_c), col = "blue", lwd = 2, xlab = expression(lambda[c]),
     main = "", ylim = c(0, max(density(prior_lambda_c)$y, density(posterior_df$lambda_c)$y)), cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(density(posterior_df$lambda_c), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Prior", "Posterior"), col = c("blue", "red"), lty = 1:2, lwd = 2)

# Delay
plot(density(prior_delay), col = "blue", lwd = 2, xlab = "Delay (months)",
     main = "", ylim = c(0, max(density(prior_delay)$y, density(posterior_df$delay_time)$y)), cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(density(posterior_df$delay_time), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Prior", "Posterior"), col = c("blue", "red"), lty = 1:2, lwd = 2)

# HR*
plot(density(prior_HR), col = "blue", lwd = 2, xlab = "HR*",
     main = "", ylim = c(0, max(density(prior_HR)$y, density(posterior_df$HR)$y)), cex.axis=1.5, cex.lab=1.5, cex.main=2)
lines(density(posterior_df$HR), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Prior", "Posterior"), col = c("blue", "red"), lty = 1:2, lwd = 2)

dev.off()





################# Calculating BPP with posteriors --------------------



max_n_c <- num_patients/2
max_n_t <- num_patients/2
max_rec_time <- 34
censoring_model <- list(method = "Events", events = 512)
analysis_model <- list(method = "LRT", alpha = 0.025, alternative_hypothesis = "one.sided")

x <- BPP_func(df, posterior_df, max_n_c, max_n_t, max_rec_time, df_cens_time = cens_time, censoring_model, analysis_model)
mean(x$BPP_df$success)
par(mfrow = c(1,1))
png("Z_Stats_BPP_DTE.png", units="in",  width=10, height=8, res=700)
hist(x$BPP_df$Z_val, xlab = "Z-Statistic", freq = F, breaks = 20, main = "Histogram of Z-Statistics")
abline(v = qnorm(0.975), lwd = 2, col = "red", lty = 2)
dev.off()


##########################################
#Code for the Examples
##########################################

#if $lambda_c$ is Gamma(14.2, 181)
#s_1(12) ~ Beta(10.2, 15.1)

#Checking
surv_beta_samples <- rbeta(50000, 10.2, 15.1)
plot(density(-log(surv_beta_samples)/12))

lambda_c_samples <- rgamma(50000, 14.2, 181)
lines(density(lambda_c_samples))



assurance_calc <- DTEAssurance::calc_dte_assurance(n_c = 400,
                                     n_t = 400,
                                     control_model = list(dist = "Exponential", 
                                                          parameter_mode = "Distribution", 
                                                          t1 = 12, 
                                                          t1_Beta_a = 10.2, 
                                                          t1_Beta_b = 15.1), 
                                     effect_model = list(delay_SHELF = SHELF::fitdist(c(3, 4, 5), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                                                         delay_dist = "gamma",
                                                         HR_SHELF = SHELF::fitdist(c(0.55, 0.6, 0.7), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                                                         HR_dist = "gamma",  
                                                         P_S = 0.9,
                                                         P_DTE = 0.8),
                                     censoring_model = list(method = "Events",
                                                            events = 650),
                                     recruitment_model = list(method = "power",
                                                              period = 12, 
                                                              power = 1),
                                     analysis_model = list(method = "LRT",
                                                           alternative_hypothesis = "one.sided",
                                                           alpha = 0.025),
                                     n_sims = 10000)
                                       



IA_Assurance_Calc <- DTEAssurance::calc_dte_assurance_interim(n_c = 400,
                                                              n_t = 400,
                                                              control_model = list(dist = "Exponential", 
                                                                                   parameter_mode = "Distribution", 
                                                                                   t1 = 12, 
                                                                                   t1_Beta_a = 10.2, 
                                                                                   t1_Beta_b = 15.1), 
                                                              effect_model = list(delay_SHELF = SHELF::fitdist(c(3, 4, 5), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                                                                                  delay_dist = "gamma",
                                                                                  HR_SHELF = SHELF::fitdist(c(0.55, 0.6, 0.7), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                                                                                  HR_dist = "gamma",  
                                                                                  P_S = 0.9,
                                                                                  P_DTE = 0.8),
                                                              recruitment_model = list(method = "power",
                                                                                       period = 12, 
                                                                                       power = 1),
                                                              GSD_model = list(events = 650,
                                                                               alpha_spending = c(rep("0.0125, 0.025", 9)),
                                                                               beta_spending = c(rep("0.05, 0.1", 9)),
                                                                               IF_vec = c("0.1, 1", "0.2, 1", "0.3, 1", "0.4, 1", 
                                                                                          "0.5, 1", "0.6, 1", "0.7, 1", "0.8, 1",
                                                                                          "0.9, 1")),
                                                              n_sims = 1000)




df2 <- IA_Assurance_Calc %>%
  mutate(
    IF_num = str_extract(IF, "^[0-9.]+") |> as.numeric(),
    success_indicator = Decision %in% c("Successful at final", "Stop for efficacy")
  )

df2 <- df2 %>%
  mutate(
    Correctness = case_when(
      # Futility
      Decision == "Stop for futility" & Final_Decision == "Unsuccessful" ~ 
        "Stopped correctly for futility",
      Decision == "Stop for futility" & Final_Decision == "Successful" ~ 
        "Stopped incorrectly for futility",
      
      # Efficacy
      Decision == "Stop for efficacy" & Final_Decision == "Successful" ~ 
        "Stopped correctly for efficacy",
      Decision == "Stop for efficacy" & Final_Decision == "Unsuccessful" ~ 
        "Stopped incorrectly for efficacy",
      
      # Final analyses
      Decision == "Successful at final"   ~ "Final success",
      Decision == "Unsuccessful at final" ~ "Final failure",
      
      # Safety fallback (should never be reached)
      TRUE ~ "Unknown"
    )
  )


summary_by_IF <- df2 %>%
  group_by(IF_num) %>%
  summarize(
    n            = n(),
    assurance    = mean(success_indicator),
    mean_duration = mean(StopTime),
    mean_sample   = mean(SampleSize),
    .groups = "drop"
  )


# Order decisions for consistent stacking
df2$Decision <- factor(
  df2$Decision,
  levels = c(
    "Unsuccessful at final",
    "Stop for futility",
    "Stop for efficacy",
    "Successful at final"
  )
)

df2$Correctness <- factor(
  df2$Correctness,
  levels = c(
    "Final failure",
    "Stopped correctly for futility",
    "Stopped incorrectly for futility",
    "Stopped incorrectly for efficacy",
    "Stopped correctly for efficacy",
    "Final success"
  )
)



assurance_df <- df2 %>%
  group_by(IF_num) %>%
  summarise(
    assurance = mean(Decision %in% c("Stop for efficacy", "Successful at final")),
    .groups = "drop"
  )



  p <- ggplot(df2, aes(x = factor(IF_num), fill = Decision)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      x = "Information Fraction",
      y = "Proportion",
      fill = "Outcome"
    ) +
    theme_minimal(base_size = 13)
  
  p +
    geom_segment(
      data = assurance_df,
      aes(
        x = as.numeric(factor(IF_num)) - 0.45,
        xend = as.numeric(factor(IF_num)) + 0.45,
        y = assurance,
        yend = assurance
      ),
      linewidth = 1.2,
      color = "black",
      inherit.aes = F
    )
  
  
  
  p <- ggplot(df2, aes(x = factor(IF_num), fill = Correctness)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      x = "Information Fraction",
      y = "Proportion",
      fill = "Outcome"
    ) +
    theme_minimal(base_size = 13)
  
  p +
    geom_segment(
      data = assurance_df,
      aes(
        x = as.numeric(factor(IF_num)) - 0.45,
        xend = as.numeric(factor(IF_num)) + 0.45,
        y = assurance,
        yend = assurance
      ),
      linewidth = 1.2,
      color = "black",
      inherit.aes = F
    )
  
  
  
  

















