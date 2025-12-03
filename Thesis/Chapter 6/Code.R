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


# Code for the example -------------------------

Scenario1_Design1 <- readRDS("Thesis/Chapter 6/data/Scenario1_Design1.rds")
Scenario1_Design2 <- readRDS("Thesis/Chapter 6/data/Scenario1_Design2.rds")
Scenario1_Design3 <- readRDS("Thesis/Chapter 6/data/Scenario1_Design3.rds")
Scenario2_Design1 <- readRDS("Thesis/Chapter 6/data/Scenario2_Design1.rds")
Scenario2_Design2 <- readRDS("Thesis/Chapter 6/data/Scenario2_Design2.rds")
Scenario2_Design3 <- readRDS("Thesis/Chapter 6/data/Scenario2_Design3.rds")
Scenario3_Design1 <- readRDS("Thesis/Chapter 6/data/Scenario3_Design1.rds")
Scenario3_Design2 <- readRDS("Thesis/Chapter 6/data/Scenario3_Design2.rds")
Scenario3_Design3 <- readRDS("Thesis/Chapter 6/data/Scenario3_Design3.rds")
Scenario4_Design1 <- readRDS("Thesis/Chapter 6/data/Scenario4_Design1.rds")
Scenario4_Design2 <- readRDS("Thesis/Chapter 6/data/Scenario4_Design2.rds")
Scenario4_Design3 <- readRDS("Thesis/Chapter 6/data/Scenario4_Design3.rds")

results_table <- data.frame(Scenario = rep(1:4, each = 4),
                            Design = rep(c(1,2,3,4), each = 4),
                            Power = numeric(16),
                            Prob_Early_Fut = numeric(16),
                            Prob_Early_Eff = numeric(16),
                            E_SS = numeric(16),
                            E_Dur
                            ation = numeric(16))

results_table$Power[1] = Scenario1_Design1$assurance
results_table$Power[2] = mean(Scenario1_Design2$Decision %in% c("Stop for efficacy", "Successful at final"))
















