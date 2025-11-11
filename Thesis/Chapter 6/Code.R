library(dplyr)
library(survival)

no_IA_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  final_time <- df$psuedo_time[num_events]
  
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  return(list(power = power,
              duration = duration,
              ss = ss))
  
}


wieand_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  interim1_time <- df$psuedo_time[num_events*0.5]
  interim2_time <- df$psuedo_time[num_events*0.75]
  final_time <- df$psuedo_time[num_events]
  
  # Evaluate at first interim
  df1 <- df[df$rec_time <= interim1_time, ]
  df1$event <- df1$psuedo_time <= interim1_time
  df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
  
  fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
  hr1 <- exp(coef(fit1))
  if (hr1 > 1){
    return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                interim = 1))
    
  }
  
  # Evaluate at second interim if not stopped

    df2 <- df[df$rec_time <= interim2_time, ]
    df2$event <- df2$psuedo_time <= interim2_time
    df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
    
    fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
    hr2 <- exp(coef(fit2))
    if (hr2 > 1){
      return(list(power = 0,  duration = interim2_time, ss = nrow(df2),
                  interim = 2))
    }

  
  
  # Evaluate at final analysis if not stopped
  
    
    df_final <- df[df$rec_time <= final_time, ]
    df_final$event <- df_final$psuedo_time <= final_time
    df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
    
    fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
    test_statistic <- summary(fit)$sctest["test"]
    observed_HR <- exp(coef(fit))
    
    power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
    duration <- final_time
    ss <- nrow(df_final)
    
  
  
  return(list(power = power,
              duration = duration,
              ss = ss, interim = 3))
  
}


OBF_func <- function(df){
  
  df <- df[order(df$psuedo_time), ]
  interim1_time <- df$psuedo_time[num_events*(1/3)]
  interim2_time <- df$psuedo_time[num_events*(2/3)]
  final_time <- df$psuedo_time[num_events]
  
  # Evaluate at first interim
  df1 <- df[df$rec_time <= interim1_time, ]
  df1$event <- df1$psuedo_time <= interim1_time
  df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
  
  fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
  hr1 <- exp(coef(fit1))
  
  
  
  if (hr1 > 0.998){
    return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                interim = 1))
    
  }
  
  # Evaluate at second interim if not stopped
  
    df2 <- df[df$rec_time <= interim2_time, ]
    df2$event <- df2$psuedo_time <= interim2_time
    df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
    
    fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
    hr2 <- exp(coef(fit2))
    if (hr2 > 0.913){
      return(list(power = 0, duration = interim2_time, ss = nrow(df2),
                  interim = 2))
    }
  
  
  
  # Evaluate at final analysis if not stopped
  
    
    df_final <- df[df$rec_time <= final_time, ]
    df_final$event <- df_final$psuedo_time <= final_time
    df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
    
    fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
    test_statistic <- summary(fit)$sctest["test"]
    observed_HR <- exp(coef(fit))
    
    power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
    duration <- final_time
    ss <- nrow(df_final)
    
  
  
  return(list(power = power,
              duration = duration,
              ss = ss, interim = 3))
}


proposed_func <- function(df) {
  df <- df[order(df$psuedo_time), ]
  
  # Event thresholds
  target1 <- as.integer(ceiling(num_events * 0.50))
  target2 <- as.integer(ceiling(num_events * 0.75))
  final_time <- df$psuedo_time[num_events]
  
  # Function to find earliest time meeting both event and delay conditions
  find_earliest_time <- function(df, target_events) {
    for (i in seq_len(nrow(df))) {
      t_i <- df$psuedo_time[i]
      df_i <- df[df$rec_time <= t_i, ]
      df_i$event <- df_i$psuedo_time <= t_i
      
      n_events <- sum(df_i$event)
      if (n_events >= target_events && n_events > 0) {
        prop_delayed <- mean((df_i$psuedo_time - df_i$rec_time)[df_i$event] >= 3)
        if (!is.na(prop_delayed) && prop_delayed >= 2/3) {
          return(list(time = t_i, prop_delayed = prop_delayed))
        }
      }
    }
    return(NULL)
  }
  
  # ---- Interim 1 ----
  ia1 <- find_earliest_time(df, target1)
  if (!is.null(ia1) && ia1$time <= final_time) {
    interim1_time <- ia1$time
    df1 <- df[df$rec_time <= interim1_time, ]
    df1$event <- df1$psuedo_time <= interim1_time
    df1$obs_time <- pmin(df1$time, interim1_time - df1$rec_time)
    
    fit1 <- coxph(Surv(obs_time, event) ~ group, data = df1)
    hr1 <- exp(coef(fit1))
    
    if (hr1 > 1) {
      return(list(power = 0, duration = interim1_time, ss = nrow(df1),
                  interim = 1, hr = hr1, prop_delayed = ia1$prop_delayed))
    }
  }
  
  # ---- Interim 2 ----
  ia2 <- find_earliest_time(df, target2)
  if (!is.null(ia2) && ia2$time <= final_time) {
    interim2_time <- ia2$time
    df2 <- df[df$rec_time <= interim2_time, ]
    df2$event <- df2$psuedo_time <= interim2_time
    df2$obs_time <- pmin(df2$time, interim2_time - df2$rec_time)
    
    fit2 <- coxph(Surv(obs_time, event) ~ group, data = df2)
    hr2 <- exp(coef(fit2))
    
    if (hr2 > 1) {
      return(list(power = 0, duration = interim2_time, ss = nrow(df2),
                  interim = 2))
    }
  }
  
  # ---- Final analysis ----
  df_final <- df[df$rec_time <= final_time, ]
  df_final$event <- df_final$psuedo_time <= final_time
  df_final$obs_time <- pmin(df_final$time, final_time - df_final$rec_time)
  
  fit <- coxph(Surv(obs_time, event) ~ group, data = df_final)
  test_statistic <- summary(fit)$sctest["test"]
  observed_HR <- exp(coef(fit))
  
  power <- (observed_HR < 1) & (test_statistic > qchisq(0.95, 1))
  duration <- final_time
  ss <- nrow(df_final)
  
  return(list(power = power, duration = duration, ss = ss, interim = 3))
}



##Reproduce the results from in the paper

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




reproduce_func <- function(scenario, rec_duration, n_sims) {
  
  wieand_list <- vector("list", length = n_sims)
  OBF_list <- vector("list", length = n_sims)
  no_IA_list <- vector("list", length = n_sims)
  proposed_list <- vector("list", length = n_sims)
  
  for (i in 1:n_sims) {
    u <- runif(n_c)
    control_times <- -log(u) / scen_list[[scenario]]$lambda_c
    
    u <- runif(n_t)
    CP <- exp(-scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1* scen_list[[scenario]]$delay)
    treatment_times <- ifelse(
      u > CP,
      -log(u) / (scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1),
      (1 / (scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR2)) * (scen_list[[scenario]]$lambda_c*scen_list[[scenario]]$HR2*scen_list[[scenario]]$delay -log(u) - scen_list[[scenario]]$lambda_c * scen_list[[scenario]]$HR1 * scen_list[[scenario]]$delay)
    )
    
    df <- data.frame(
      time = c(control_times, treatment_times),
      group = c(rep("Control", n_c), rep("Treatment", n_t))
    )
    
    df$rec_time <- runif(n_c + n_t, 0, rec_duration)
    df$psuedo_time <- df$time + df$rec_time
    
    wieand_list[[i]] <- wieand_func(df)
    OBF_list[[i]] <- OBF_func(df)
    no_IA_list[[i]] <- no_IA_func(df)
    proposed_list[[i]] <- proposed_func(df)

    

  }
  
  return(list(no_IA_list = no_IA_list,
              wieand_list = wieand_list,
              OBF_list = OBF_list,
              proposed_list = proposed_list))
  
}


outcome_df <- data.frame(matrix(nrow = 12, ncol = 12))

for (i in 1:6){
  result <- reproduce_func(i, rec_duration = 34, n_sims = 100)
  
  outcome_df[i,1] <- mean(sapply(result$no_IA_list, function(x) x$power))
  outcome_df[i,2] <- mean(sapply(result$no_IA_list, function(x) x$duration))
  outcome_df[i,3] <- mean(sapply(result$no_IA_list, function(x) x$ss))
  
  
  outcome_df[i,4] <- mean(sapply(result$wieand_list, function(x) x$power))
  outcome_df[i,5] <- mean(sapply(result$wieand_list, function(x) x$duration))
  outcome_df[i,6] <- mean(sapply(result$wieand_list, function(x) x$ss))
  sapply(result$wieand_list, function(x) x$interim)
  
  outcome_df[i,7] <- mean(sapply(result$OBF_list, function(x) x$power))
  outcome_df[i,8] <- mean(sapply(result$OBF_list, function(x) x$duration))
  outcome_df[i,9] <- mean(sapply(result$OBF_list, function(x) x$ss))
  sapply(result$OBF_list, function(x) x$interim)
  
  outcome_df[i,10] <- mean(sapply(result$proposed_list, function(x) x$power))
  outcome_df[i,11] <-  mean(sapply(result$proposed_list, function(x) x$duration))
  outcome_df[i,12] <-  mean(sapply(result$proposed_list, function(x) x$ss))
  sapply(result$proposed_list, function(x) x$interim)
  
  
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


result <- reproduce_func(9, rec_duration = 12, n_sims = 100)


