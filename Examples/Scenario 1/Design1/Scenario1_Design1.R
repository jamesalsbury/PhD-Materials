#!/usr/bin/env Rscript

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


Scenario <- 1
Design <- 1






res <- calc_dte_assurance(
  n_c = n_c,
  n_t = n_t,
  control_model = control_model,
  effect_model = effect_model,
  censoring_model = censoring_model,
  recruitment_model = recruitment_model,
  analysis_model = analysis_model,
  n_sims = n_sims
)

# ---------------------
# Save output
# ---------------------

outname <- sprintf("assurance_batch_%04d.rds", seed)
saveRDS(res, outname)
