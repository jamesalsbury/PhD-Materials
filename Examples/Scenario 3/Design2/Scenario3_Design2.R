#!/usr/bin/env Rscript

library(rjags)
library(survival)
library(dplyr)
library(DTEAssurance)
library(jsonlite)
library(SHELF)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript Scenario1_Design2.R <n_sims> <seed>")
}

n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])

set.seed(seed)

# ---------------------
# Mode setup (unchanged)
# ---------------------

## Design 1 -----------------
n_c = 400
n_t = 400

control_model <- list(dist = "Exponential",
                      parameter_mode = "Fixed",
                      fixed_type = "Parameters",
                      lambda = log(2)/9)


effect_model <- list(delay_SHELF = SHELF::fitdist(c(3.99, 4, 4.1), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(0.599, 0.6, 0.601), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                     HR_dist = "gamma",  
                     P_S = 1,
                     P_DTE = 0)

recruitment_model <- list(method = "power",
                          period = 12,
                          power = 1)


GSD_model <- list(events = 650,
                  alpha_spending = c(0.0125, 0.025),
                  alpha_IF = c(0.75, 1),
                  futility_type = "none")

res <- calc_dte_assurance_interim(n_c = n_c,
                                  n_t = n_t,
                                  control_model = control_model,
                                  effect_model = effect_model,
                                  recruitment_model = recruitment_model,
                                  GSD_model = GSD_model,
                                  n_sims = n_sims)

# ---------------------
# Save output
# ---------------------

outname <- sprintf("assurance_batch_%04d.rds", seed)
saveRDS(res, outname)

