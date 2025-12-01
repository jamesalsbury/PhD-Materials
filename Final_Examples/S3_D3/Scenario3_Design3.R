library(rjags)
library(survival)
library(dplyr)
library(DTEAssurance)
library(jsonlite)
library(SHELF)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript Scenario1_Design1.R <n_sims> <seed>")
}

n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])

set.seed(seed)

n_c <- 400
n_t <- 400

fixed_params <- list(lambda_c = log(2)/9,
                     delay_time = 0,
                     post_delay_HR = 0.6
)


control_model <- list(dist = "Exponential",
                      parameter_mode = "Distribution",
                      t1 = 12,
                      t1_Beta_a = 10.2,
                      t1_Beta_b = 15.1)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(0.55, 0.6, 0.7), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                     HR_dist = "gamma",
                     P_S = 0.9,
                     P_DTE = 0.8)

recruitment_model <- list(method = "power",
                          period = 24,
                          power = 1)

GSD_model <- list(events = 650,
                  alpha_spending = c(0.0125, 0.025),
                  alpha_IF = c(0.75, 1),
                  futility_type = "BPP",
                  futility_IF = 0.5,
                  BPP_threshold = 0.2)


analysis_model <- list(method = "LRT",
                       alternative_hypothesis = "one.sided",
                       alpha = 0.025)

fixed_parameters_BPP <- function(n_c,
                                 n_t,
                                 fixed_params,
                                 control_model,
                                 effect_model,
                                 recruitment_model,
                                 GSD_model,
                                 analysis_model,
                                 n_sims = 100) {
  
  results <- future.apply::future_lapply(seq_len(n_sims), function(i) {
    
    
    
    
    trial <- sim_dte(n_c, n_t, fixed_params$lambda_c, fixed_params$delay_time, fixed_params$post_delay_HR)
    
    trial <- add_recruitment_time(trial,
                                  rec_method = recruitment_model$method,
                                  rec_period = recruitment_model$period,
                                  rec_power = recruitment_model$power,
                                  rec_rate = recruitment_model$rate,
                                  rec_duration = recruitment_model$duration)
    
    trial_data <- trial[order(trial$pseudo_time),]
    n_events   <- GSD_model$events
    t_interim  <- trial_data$pseudo_time[n_events]
    
    eligible_df <- trial_data %>%
      dplyr::filter(.data$rec_time <= t_interim)
    
    # Censoring logic
    eligible_df$status <- eligible_df$pseudo_time < t_interim
    eligible_df$survival_time <- ifelse(
      eligible_df$status,
      eligible_df$time,
      t_interim - eligible_df$rec_time
    )
    
    # Interim Cox model (for final decision)
    fit         <- survival::coxph(Surv(survival_time, status) ~ group, data = eligible_df)
    fit_summary <- summary(fit)
    z_stat      <- -fit_summary$coefficients[, "z"]
    
    
    # Build rpact design (still needed for efficacy boundaries)
    rpact_design <- DTEAssurance:::make_rpact_design_from_GSD_model(GSD_model)
    design       <- rpact_design$design
    
    outcome <- DTEAssurance:::apply_GSD_to_trial(
      n_c = n_c,
      n_t = n_t,
      trial_data        = trial,
      design            = design,
      total_events      = GSD_model$events,
      GSD_model         = GSD_model,
      control_model     = control_model,        # elicited prior
      effect_model      = effect_model,         # elicited prior
      recruitment_model = recruitment_model,
      analysis_model    = analysis_model,
      n_BPP_sims        = 50                    # or whatever you choose
    )
    
    return(data.frame(
      Trial          = i,
      Decision       = outcome$decision,
      StopTime       = outcome$stop_time,
      SampleSize     = outcome$sample_size,
      Final_Decision = ifelse(
        z_stat > stats::qnorm(1 - 0.025),
        "Successful",
        "Unsuccessful"
      )
    ))
    
  }, future.seed = TRUE)
  
  # Combine into single data frame
  results_df <- do.call(rbind, results)
  
  return(results_df)
  
}

res <- fixed_parameters_BPP(n_c = n_c,
                            n_t = n_t,
                            fixed_params = fixed_params,
                            control_model = control_model,
                            effect_model = effect_model,
                            recruitment_model = recruitment_model,
                            GSD_model = GSD_model,
                            analysis_model = analysis_model,
                            n_sims = n_sims)


outname <- sprintf("assurance_batch_%04d.rds", seed)
saveRDS(res, outname)

