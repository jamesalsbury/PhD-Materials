#!/usr/bin/env Rscript

library(rjags)
library(survival)
library(dplyr)
library(DTEAssurance)
library(jsonlite)
library(SHELF)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript Scenario1_Design3.R <n_sims> <seed>")
}

n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])

set.seed(seed)

Scenario <- 1
Design <- 3

model_inputs <- DTEAssurance:::generate_model_inputs(Scenario, Design)

res <- calc_dte_assurance_interim(n_c = model_inputs$n_c,
                                  n_t = model_inputs$n_t,
                                  control_model = model_inputs$control_model,
                                  effect_model = model_inputs$effect_model,
                                  recruitment_model = model_inputs$recruitment_model,
                                  GSD_model = model_inputs$GSD_model,
                                  analysis_model = model_inputs$analysis_model,
                                  n_sims = n_sims)

outname <- sprintf("assurance_batch_%04d.rds", seed)
saveRDS(res, outname)









